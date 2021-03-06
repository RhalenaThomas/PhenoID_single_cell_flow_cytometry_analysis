# one the fly analysis to fine tune gating 
# gating strategy already defined
# adjust exact cut off from real time data


#1. Preprocess - save only the wanted datatype = retro transformed
#2. Cluster - with pre-estimated best conditions = say we acquire 50,000 - pretest clustering conditions - export feature maps and heatmaps
#3. Run label transfer, correlation, Random Forest  
#4. Review output - label clusters - manually - have code set up to go
#Run hypergate for predetermined groups AND get exact levels for the new data.


#### create flowset object for clustering and labeling #########
# libraries
# load libraries

######## Clustering and visualization ##########################

library(Seurat)
library(ggplot2)
library(clustree)
library(reshape2) #for plotting multiple lines (resolutions) on the same graph
library(dplyr)
library(kit) # for finding max and second max (function topn)
library(tidyr) #for the last plot in the function


# I prepared the csv 
input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/GatingPlanExperiment/Output/flowset1.csv"

input_name <- "May10_read1"
clust_method <- "Louvain" # cluster type for file name
df <- read.csv(input_path) # read in the dataframe
df2 <- df %>% dplyr::select(c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4"))
# the order of the DF is set by the order the columns are written above 
# create a matrix for later
m <- as.matrix(df2)
# create the seurat object for visualization
tm <- t(df2)
rownames(tm) <- colnames(df2)
colnames(tm) <- rownames(df2)
seu <- CreateSeuratObject(tm)
# for later to see features we define the vector here
AB <- colnames(df2)
# add the sample names into the seurat object 
seu <- AddMetaData(object=seu, metadata=df$Sample, col.name = 'Sample')
#downsample to no have too many cells
# check the cell numbers first
table(seu$Sample)

#  S1-A  S1-B  S2-A  S2-B 
# 24359 22660 22074 20312 
# I have about twice as many cells as expected
# I'll try without downsampling first

# prepare the object for clustering
seu <- ScaleData(seu)
# assume we acquire 50,000 cells 
seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)
seu <- FindNeighbors(seu, dims = 1:12, k.param =60)
seu <- RunUMAP(seu, dims = 1:12, n.neighbors = 60, spread = 2, a=0.54,b=0.84, min.dist = 0.5)
# note to me - test these resolutions and pick one here 
seu <- FindClusters(seu, resolution = c(0.8,1.2,1.5))
DimPlot(seu, group.by = 'RNA_snn_res.1.2')

########## label the clusters   ######################################
# use the 9000 labels used to create these groupings
seu.r <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/Louvkn60DifferentCellLabels220220318.Rds")

DimPlot(seu.r, group.by = 'ref.cells')
DimPlot(seu.r, group.by = 'gating.all')
Idents(seu.r) <- 'ref.cells'
seu.r <- subset(x= seu.r, downsample = 500)

# create anchors
anchors <- FindTransferAnchors(reference = seu.r, query = seu,features = AB ,reference.reduction = "pca", dim= 1:10) 
predictions <- TransferData(anchorset = anchors, refdata = seu.r$ref.cells, dims = 1:10)
# add the predicitons to the seurat object
seu <- AddMetaData(seu, metadata = predictions)
# get tables to see and plots 
t.lables <- as.data.frame(table(seu$seurat_clusters, seu$predicted.id))
pr.t.lables <- as.data.frame(prop.table(table(seu$seurat_clusters, seu$predicted.id)))
t.lables$Freq <- as.double(t.lables$Freq)
# plot to check out the data
ggplot(t.lables, aes(y = Freq, x = Var1, fill = Var2)) + geom_bar(position = "stack", stat= "identity")
DimPlot(seu, group.by = 'predicted.id')
# see tables of top labels
top.labs <- t.lables  %>% group_by(Var1)  %>% top_n(5, Freq)
top.labs
top.lab <- t.lables  %>% group_by(Var1)  %>% top_n(1, Freq)
top.lab


# see the different samples and expression
DimPlot(seu, split.by = 'Sample', group.by = 'predicted.id')
FeaturePlot(seu, split.by = 'Sample', features = "CD24", slot = "scale.data", min.cutoff = 0.25)
FeaturePlot(seu, split.by = 'Sample', features = "CD56", slot = "scale.data", min.cutoff = 0.25)
FeaturePlot(seu, split.by = 'Sample', features = "CD15", slot = "scale.data", min.cutoff = 0.25)
FeaturePlot(seu, split.by = 'Sample', features = "CD29", slot = "scale.data", min.cutoff = 0.25)
FeaturePlot(seu, split.by = 'Sample', features = "CD184", slot = "scale.data", min.cutoff = 0.25)
FeaturePlot(seu, split.by = 'Sample', features = "CD133", slot = "scale.data", min.cutoff = 0.25)
FeaturePlot(seu, split.by = 'Sample', features = "O4", slot = "scale.data", min.cutoff = 0.25)
FeaturePlot(seu, split.by = 'Sample', features = "GLAST", slot = "scale.data", min.cutoff = 0.25)
FeaturePlot(seu, split.by = 'Sample', features = "AQP4", slot = "scale.data", min.cutoff = 0.25)


DotPlot(seu, features = AB)


# now add the labels from the clusters and also relabel with only the gating groups
Idents(seu) <- "seurat_clusters"
cluster.ids <- c("Astro-AQ","Mix","Neur-CD56","Glia-CD44","Neur-CD24","Neur-Glia","Mix","Glia-CD44","RG-CD184","RG-CD133","Neur-CD15",
                  "RG-CD184","Neur-CD56","Glia-glast","Astro-CD71")
names(cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, cluster.ids)
seu$labels <- Idents(seu)
DimPlot(seu, reduction = "umap", label = TRUE, group.by = 'labels', repel = TRUE)


# then label with the main groups
seu.r <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/Louvkn60DifferentCellLabels220220318.Rds")
Idents(seu.r) <- 'gating.main'
seu.r <- subset(x= seu.r, downsample = 500)
DimPlot(seu.r)

# create anchors
anchors <- FindTransferAnchors(reference = seu.r, query = seu,features = AB ,reference.reduction = "pca", dim= 1:10) 
predictions <- TransferData(anchorset = anchors, refdata = seu.r$gating.main, dims = 1:10)
# add the predicitons to the seurat object
seu <- AddMetaData(seu, metadata = predictions)
# get tables to see and plots 
t.lables <- as.data.frame(table(seu$seurat_clusters, seu$predicted.id))
pr.t.lables <- as.data.frame(prop.table(table(seu$seurat_clusters, seu$predicted.id)))
t.lables$Freq <- as.double(t.lables$Freq)
# plot to check out the data
ggplot(t.lables, aes(y = Freq, x = Var1, fill = Var2)) + geom_bar(position = "stack", stat= "identity")
DimPlot(seu, group.by = 'predicted.id')
# see tables of top labels
top.labs <- t.lables  %>% group_by(Var1)  %>% top_n(5, Freq)
top.labs
top.lab <- t.lables  %>% group_by(Var1)  %>% top_n(1, Freq)
top.lab


##### correlation labels
# run the function or load library 
# input path above
reference_path <- "/Users/rhalenathomas/GITHUB/PhenoID_single_cell_flow_cytometry_analysis/correlation/ReferenceMatrix9celltypesOrdered.csv"
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/GatingPlanExperiment/Output/"

find_correlation(input_path, reference_path, output_path)
cor <- read.csv("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/GatingPlanExperiment/Output/correlationscorr_celltypes.csv")



sue <- AddMetaData(object=seu.q, metadata=cor$cell.label, col.name = 'cor.labels')
unique(seu.q$cor.labels)



# add the main group labels
Idents(seu) <- "seurat_clusters"
cluster.ids <- c("Neurons-CD56","RG-CD184","Neurons-mix1",
                 "RG-CD184","Neurons-mix2","RG-CD184",
                 "RG-CD184","Neurons-CD56","Neurons-mix1",
                 "Mix","RG-CD184","RG-CD184",
                 "RG-CD184","RG-CD184","RG",
                 "Neurons-CD56","RG-CD184","Mix",
                 "Neurons-CD56","Mix","Mix",
                 "Mix","RG-CD184","RG-CD184",
                  "RG-CD184","Neurons-CD24","Mix",
                 "oligo","RG-CD184")
names(cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, cluster.ids)
seu$gating.groups <- Idents(seu)
DimPlot(seu, reduction = "umap", label = TRUE, group.by = 'gating.groups', repel = TRUE)

Idents(seu) <- "seurat_clusters"
cluster.ids <- c("Neurons-CD56","RG-CD184-a","Neurons-mix1",
                 "RG-CD184-a","Neurons-CD24","RG-CD1843",
                 "RG-CD184-a","Neurons-CD56","Neurons-mix1",
                 "Mix","RG-CD184-N","RG-CD184-N",
                 "RG-CD1847","RG-CD1848","RG",
                 "Neurons-CD56","RG-CD184","Mix",
                 "Neurons-CD56","Mix","Mix",
                 "Mix","RG-CD1849","RG-CD18410",
                 "RG","Neurons-CD24","Mix",
                 "oligo","RG-CD18412")
names(cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, cluster.ids)
seu$gating.groups2 <- Idents(seu)
DimPlot(seu, reduction = "umap", label = TRUE, group.by = 'gating.groups2', repel = TRUE)


#### save the labelled object ######
saveRDS(seu,paste(output_path,"Test_May_10_seuLabels.Rds"))



################################ apply hypergates ####################################
library("hypergate")


Idents(seu) <- 'gating.groups'
i = 800
set.seed(i)
seu.down <- subset(x= seu, downsample = i)
DimPlot(seu.down)

input.xm = as.matrix(GetAssayData(seu.down, slot = 'data'))
xm.t <- t(input.xm)
cluster.labels <- as.vector(seu.down$gating.groups)

#cell types to gate

cell.types <- c("Neurons-CD56", "Neurons-CD24","Neurons-mix1","Neurons-mix2","RG-CD184","Mix")


for (cell in cell.types) {
  hg_output <- hypergate(xp = xm.t, gate_vector = cluster.labels, level = cell)
  gating_predicted = subset_matrix_hg(hg_output,xm.t)
  conf.table <- table(ifelse(gating_predicted, "Gated-in", "Gated-out"), ifelse(cluster.labels ==  cell, cell, "Others"))
  print(cell)
  print(conf.table)
  print(hgate_rule(hg_output))
 
}


# adjust and write out gating strategy


