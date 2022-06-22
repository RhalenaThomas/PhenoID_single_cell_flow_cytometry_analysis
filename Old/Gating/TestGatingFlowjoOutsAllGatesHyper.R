# I have created gates from the hypergate output of my annotations
# I gated in Flowjo and saved the fsc files concatonating across all 9MBO samples.
# 


# test the files that I gated:
# 1. Preprocess - save retrotrans only
# create the seurat object
# look at the expression across the inputs 
# 2. Count total cells and cluster with Louvain only
# 3. Use the annotation 3 methods
# 4. Annotate the clusters
# see how well they match

# there are 4 files fsc from Flowjo: Neurons1, Neurons2, Astrocytes, Radial Glia

input_path = "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/FlowDataFiles/Exported_gating_9MBO/FinalGatesApril18"

################### preprocessing ######################
# load libraries

require("flowCore") #Used for reading the data
require("ggplot2")
require("ggridges") #visualization
require("stringr") #set of functions to manipulate strings type in order to have nice titles
require("rlist") #set of functions that allows to easely manipulate lists
require("reshape2") #visualization
require("flowStats") #Alignment functions
require("scales") #scale colour intensity for visualization
require("dplyr")

library("flowCore")
library("flowStats")

# run functions
# functions Alex wrote
# will need to make these into a separate function with documentation

# function plotdensity_flow set

plotdensity_flowset <- function(flowset){ ggplot(melt(lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs)})), aes(x=value,y=L1,fill=L1)) + geom_density_ridges(alpha=.4,verbose=FALSE) +facet_wrap(~variable)+theme_light()} 

#defines a function for visualizing flowset with densityplots

# function renmame markers
# fcs files will have the marker names input during acquistion using flowjo

rename_markers<-function(flowset){#Defines a function to use marker names 
  copy_flowset=flowset[seq(along=flowset)]
  for (i in 1:length(copy_flowset)){
    marker.names=copy_flowset[[i]]@parameters@data$desc
    marker.names=lapply(marker.names,function(x){str_replace_all(x,"-","_")})
    colnames(copy_flowset[[i]]@exprs)<-unlist(lapply(marker.names, function(x){sapply(str_split(x,"_"),head,1)})) 
  }
  return(copy_flowset)
}


### user input
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/FlowDataFiles/Exported_gating_9MBO/sample_on_the_fly/flowset/"


write_fcs_files <- FALSE #Set to true to write fcs files at each step (subsampled, transformed, aligned and scaled)
# I don't want to save these files at this time 

#Create output folder if it's not already created
if(dir.exists(output_path)==FALSE){ #check
  dir.create(output_path) #create directory
}


input_path = "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/FlowDataFiles/Exported_gating_9MBO/sample_on_the_fly/flowjo"

# all the files in the input folder will be added to this data object flowset
flowset <- read.flowSet(path=input_path,transformation = FALSE ,emptyValue = FALSE,truncate_max_range = FALSE, package="flowCore") #Create a flowset object


flowset <- fsApply(flowset,function(x){x[ ,grepl("^[{'FJComp'}|{'FCS'}|{'SSC'}].*A$",colnames(flowset))]}) #Apply a function to flowSet (fsApply) that selects only the columns corresponding to areas values using "regular expression".


# see number of cells per sample
print(fsApply(flowset,function(x){dim(x@exprs)})[,1])

sampleNames(flowset) #Prints the name of each flowframe inside the flowset. User can modify it in the following chunck

# rename sample names

#sampleNames(flowset) <- c("Astro1","Astro2","Astro3","Endothelial","Epithelial","Neurons1", "Neurons2", "Oligo","RG1","RG2") #directly modify the names inside the flowset
#sampleNames(flowset) <- c("Astro","CD24++CD56++","CD44++CD184++","Glia","Neurons1-CD56","Neurons2-CD24","other_glia","RG-topgate")
# April 18 final
sampleNames(flowset) <- c("Astro-CD44+","Glia-CD184+","Neurons1-CD56+","Neurons2-CD24+")
sampleNames(flowset)

#### preprocessing - do the transformation

inversebiexponentialTransform <- function(flowset,a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0){
  copy_flowset=flowset[seq(along=flowset)] #makes a copy of the input flowset
  for (i in 1:length(copy_flowset)){ #loop though index to get each flowframe
    copy_flowset[[i]]@exprs=a*exp(b*(copy_flowset[[i]]@exprs-w))-c*exp(-d*(copy_flowset[[i]]@exprs-w))+f
  }
  return(copy_flowset)
}

biexp  <- biexponentialTransform("biexp transform",a = 0.5, b = 1, c = 0.5, d = 1, f = 0, w = 0) #creates the transformation object (the values make it equivalent to arcsinh transform)
transformed_flowset <- transform(flowset, transformList(colnames(flowset), biexp)) #Apply the transformation

if(write_fcs_files==TRUE){#Check if user set the conditional for writing several fcs files
  write.flowSet(transformed_flowset,outdir=paste0(output_path,"transformed_flowset"))#writes the flowset
  transformed_flowset=read.flowSet(path=paste0(output_path,"transformed_flowset"),phenoData = "annotation.txt")
}



# align the data
normtr=gaussNorm(transformed_flowset,colnames(transformed_flowset)[c(3,5:6,9:length(colnames(transformed_flowset)))],max.lms = 2,peak.density.thr = 0.01) #Detects and align 2 peaks on the marker 3,5,6,9...14. 
expbe_norm2=normtr$flowset
normtr=gaussNorm(expbe_norm2,colnames(expbe_norm2)[c(4,7:8)],max.lms = 1,peak.density.thr = 0.05)#Detects and align 1 peak 
aligned_transformed_flowset=normtr$flowset
retrotransformed_flowset <- inversebiexponentialTransform(aligned_transformed_flowset) #apply the function for cancelling biexp transform 

normtr <- gaussNorm(flowset = transformed_flowset,channel.names = colnames(transformed_flowset)[c(3:length(colnames(transformed_flowset)))], max.lms = 2, peak.density.thr = 0.01, peak.distance.thr = 0.5) #Detects and align 2 peaks (max.lms) for all markers except FSC and SSC (columns 3 to number of markers) the threshold are data-dependant and may have to differ from one analysis to another.

aligned_transformed_flowset <- normtr$flowset #Extract the flowset from the result of the alignment
retrotransformed_flowset <- inversebiexponentialTransform(aligned_transformed_flowset) #apply the function for cancelling biexp transform

# have a look
plotdensity_flowset(aligned_transformed_flowset)

plotdensity_flowset(transformed_flowset)

plotdensity_flowset(flowset)

# saving function
flowset_to_csv=function(flowset){  
  list_of_flowframes=fsApply(rename_markers(flowset),function(x){as.data.frame(x@exprs)})#Makes a list of dataframes
  list_names=names(list_of_flowframes)#extract names for not calling names() function at each loop
  for (index in seq_along(list_of_flowframes)){ #Iterates along the index for adding sample names
    list_of_flowframes[[index]]=list_of_flowframes[[index]] %>% #Using tidyverse package for adding features (name of batch)
      mutate(Batch=list_names[index])#Using tidyverse package for adding features (name of batch)
    #rename the columns to fix the rbind error 
    colnames(list_of_flowframes[[index]]) = colnames(list_of_flowframes[[1]])  
  }
  # this is wehre the error occurs but all the df have the same column names???
  ds=list.rbind(list_of_flowframes)#binds every fcs file in a single dataframe
  ds$cell=as.factor(unlist(lapply(as.list(c(1:length(flowset))),function(x){c(1:nrow(flowset[[x]]@exprs))})))#add cell IDs - cell count per sample 
  write.csv(ds,file=paste0(output_path,deparse(substitute(flowset)),".csv"))#save the R data for further usage
}


# use saving function flowset_to_csv to save the flowset retrotransformed file
flowset_to_csv(flowset)#apply the function
# save the subsset

# we will read in the data below


######## Clustering and visualization ##########################

library(Seurat)
library(ggplot2)
library(clustree)
library(reshape2) #for plotting multiple lines (resolutions) on the same graph
library(dplyr)


input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/Gating/FinalGatesApril18/flowset.csv"


input_name <- "Final4gated"  # processing type for file name didn't align
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

# add the cell type names
df$GatedCells <- as.factor(df$Batch) # so that when added to seurat object batch has levels
# add the meta data back in for sample groups
seu <- AddMetaData(object=seu, metadata=df$GatedCells, col.name = 'GatedCells')
# create the vector for the antibodies names for feature plotting later
AB <- colnames(df2)
# add to scale data slot
seu <- ScaleData(seu)


# make a heatmap of expression scaled
pdf(paste(output_path,input_name,clust_method,"Heatmap_gated_all_CellsInput.pdf",sep=""),width =8, height = 6)
print(DoHeatmap(seu, group.by = "GatedCells", features = AB))
dev.off()


# cluster  use the 'best' parameters: the ones for labeling  kn=60 res = 1
# this didn't work well before - maybe I should run the whole clustering functions
# I have 134,381 cells sqrt = 373, I'll use 300


#seu <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/Gating/april12-gates/gatingtree/sueSeu.downsample.Rds")

# next input - more cells squrt is 457 - maybe that is needed to find groups properly
# better option is do down sample

#seu <- subset(x = seu, downsample = 2000)# I ran this 
# visualize clusters and heatmap
# pdf("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/Gating/eachgate/UMAPs.pdf")

# save obeject in case of rerunning
#saveRDS(seu,  "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/Gating/eachgate/seuGatingallhypergates.Rds")
# saved over this object with the two neurons object

# number of cells April 16th 59393, squrt 243 

# clustering is not great 
# cluster again
# number of cells 260496 
# squrt 510 

# April 18 all cells 53715 from AIW only - sert = 231
# clusters are too overlaping I used k.par 100 not sure if it's better

Idents(seu) <- "GatedCells"

seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)
seu <- FindNeighbors(seu, dims = 1:12, k.param =100)
seu <- RunUMAP(seu, dims = 1:12, n.neighbors = 100)
seu <- FindClusters(seu, resolution = c(0.5,0.8,1,1.2))
seu <- FindClusters(seu, resolution = 1.8)
DimPlot(seu, group.by = 'GatedCells')
DimPlot(seu, group.by = 'seurat_clusters')


inputs <- "Neu1.Neur2.Astro.Glia"

# save the object without saving over it all the time
saveRDS(seu, paste(output_path,inputs,"seu.Rds",sep = ""))

AB <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")


pdf(paste(output_path,"UMAPs",inputs, ".pdf",sep=""))
  DimPlot(seu, group.by = 'GatedCells')
  DimPlot(seu, group.by = 'RNA_snn_res.0.5')
  DimPlot(seu, group.by = 'RNA_snn_res.0.8')
  DimPlot(seu, group.by = 'RNA_snn_res.1')
  DimPlot(seu, group.by = 'RNA_snn_res.1.2')
dev.off()

# use the res.1.2
test.vec <- seq(0.1,1.5,0.1)
# annotate clusters

# use the seurat object with the labels saved for gating as the reference
#seu.r <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outs/clusters/AllcellLablesMarch25.Rds")
# use the 9000 labels used to create these groupings
seu.r <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/Louvkn60DifferentCellLabels220220318.Rds")

# down sample reference
Idents(seu.r) <- 'labels5'
#Idents(seu.r) <- 'gating.all'
seu.r <- subset(x= seu.r, downsample = 500)

DimPlot(seu.r)

table(seu.r$labels5)
table(seu.r$gating.all)

# this step is computationally intensive and take a long time
# remember to down sample!!!
# predicting again with Labels6
anchors <- FindTransferAnchors(reference = seu.r, query = seu,features = AB ,reference.reduction = "pca", dim= 1:10) 


predictions <- TransferData(anchorset = anchors, refdata = seu.r$labels5, dims = 1:10)

# add the predicitons to the seurat object
seu <- AddMetaData(seu, metadata = predictions)


# check the seurat label transfer predicted lables
t.lables <- as.data.frame(table(seu$RNA_snn_res.1.2, seu$predicted.id))
pr.t.lables <- as.data.frame(prop.table(table(seu$RNA_snn_res.1.2, seu$predicted.id)))
t.lables$Freq <- as.double(t.lables$Freq)



# try bar chart
ggplot(t.lables, aes(y = Freq, x = Var1, fill = Var2)) + geom_bar(position = "stack", stat= "identity")

DimPlot(seu, group.by = 'predicted.id')
# find the top 2 in cell in the group

pdf(paste(output_path,"predict.gated",inputs, ".pdf",sep="", hieght = 5, width = 18))
DimPlot(seu, group.by = 'predicted.id', split.by = 'GatedCells')
dev.off()

top.labs <- t.lables  %>% group_by(Var1)  %>% top_n(5, Freq)
top.labs

top.lab <- t.lables  %>% group_by(Var1)  %>% top_n(1, Freq)
top.lab



# check the predictions within the gating groups
# check the seurat label transfer predicted lables
t.lables.gates <- as.data.frame(table(seu$GatedCells, seu$predicted.id))
t.lables.gates$Freq <- as.double(t.lables.gates$Freq)



# try bar chart
ggplot(t.lables.gates, aes(y = Freq, x = Var1, fill = Var2)) + geom_bar(position = "stack", stat= "identity")+ RotatedAxis()

# find the top in cell in the group


top.gates <- t.lables.gates  %>% group_by(Var1)  %>% top_n(1, Freq)
top.gates.5 <- t.lables.gates  %>% group_by(Var1)  %>% top_n(5, Freq)



# check which fsc files are in each cluster
t.lables.gates <- as.data.frame(table(seu$seurat_clusters, seu$GatedCells))
t.lables.gates$Freq <- as.double(t.lables.gates$Freq)



# try bar chart
ggplot(t.lables.gates, aes(y = Freq, x = Var1, fill = Var2)) + geom_bar(position = "stack", stat= "identity")+ RotatedAxis()

# find the top in cell in the group


top.gates <- t.lables.gates  %>% group_by(Var1)  %>% top_n(1, Freq)
top.gates.5 <- t.lables.gates  %>% group_by(Var1)  %>% top_n(5, Freq)


#seu <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/Gating/gatintreeApril16/suegatingtree.Rds")


# annotate the gated cells based on clustering

Idents(seu) <- "seurat_clusters"
cluster.ids <- c("Neurons1-CD56","Neurons2-CD24","Mix-RG-Neurons","Astro2","RG1","Mix-neurons","Astro1","RG1","Neurons1-CD56","Astro2","Astro1","Astro1","RG2","Neurons1-CD56","Astro1","Neurons1-CD56","Neurons2-CD24","Epilthelial")

names(cluster.ids) <- levels(seu)
seu <- RenameIdents(seu, cluster.ids)
seu$labels <- Idents(seu)
# save first predicted IDs
#Idents(seu) <- 'predicted.id'
#seu$predicted.id.1 <- Idents(seu)

DimPlot(seu, reduction = "umap", label = TRUE, group.by = 'labels', repel = TRUE)

# with annotation see what is in the gates

t.lables.gates <- as.data.frame(table(seu$GatedCells, seu$labels))
t.lables.gates$Freq <- as.double(t.lables.gates$Freq)

ggplot(t.lables.gates, aes(y = Freq, x = Var1, fill = Var2)) + geom_bar(position = "stack", stat= "identity")+ RotatedAxis()

top.gates <- t.lables.gates  %>% group_by(Var1)  %>% top_n(1, Freq)
top.gates.5 <- t.lables.gates  %>% group_by(Var1)  %>% top_n(5, Freq)



# save results - 

pdf(paste(output_path,input_name,"LablesinGated-fcs.pdf",sep=""))
DimPlot(seu, reduction = "umap", label = TRUE, group.by = 'labels', repel = TRUE)
ggplot(t.lables.gates, aes(y = Freq, x = Var1, fill = Var2)) + geom_bar(position = "stack", stat= "identity")+ RotatedAxis()
dev.off()


# make some heatmaps to try and understand the other populations
pdf(paste(output_path,input_name,"heatmaps.pdf",sep=""))
DoHeatmap(seu, group.by = 'predicted.id', features = AB) + theme(title = "Predicted.id")
DoHeatmap(seu, group.by = 'labels', features = AB)+ theme(title = "lables")
DoHeatmap(seu, group.by = 'seurat_clusters', features = AB)+ theme(title = "clusters")
dev.off()


# two grouping heatmap
# run function to use this
# funciton gpar in grid
library(grid)

DoMultiBarHeatmap(seu, assay = 'RNA',features = AB, group.by='GatedCells', 
                    additional.group.by = 'predicted.id.1', label = FALSE)
unique(seu$GatedCells)


DoMultiBarHeatmap(seu, features=AB, group.by='GatedCells', 
                  additional.group.by = 'predicted.id.1',additional.group.sort.by = 'predicted.id.1',label = FALSE)

