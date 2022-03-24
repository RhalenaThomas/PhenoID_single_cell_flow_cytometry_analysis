# flowsom clustering
# parameter tuning, statistic and visualizations for manual annotation


# @Shuming: make this whole script into a function.
# input option:  input_path, output_path, input_name, cluster_method, kn, resolutions

# input_name and cluster_method are for file names to keep track of the output
# kn and resolution needs to take a vector
# an option to turn off or on save plots would be helpful too.

# add plot with x axis as the number of clusters, kn by colour intensity, res by shape or different colour 
# values as dots, y axis statistic

# need to get the number of clusters at each kn/res combination


# load libraries
library(clusterSim) 
library(FlowSOM)
library(flowCore)
library(cluster)
library(fpc)
library(clv)
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)
library(reshape2) #for plotting multiple lines (resolutions) on the same graph

############# set up the data object for clustering ############################

# info to change for each comparison
# define the input pathway
# input pathway
# input pathway
input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outs_jan19allcells/prepro_outsretrotransformed_flowset.csv"
# output pathway
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outs_jan19allcells/Clusters/Louvain/"
# add input description to output files
input_name <- "retrotrans"  # this will be the different processing types

# cluster type for file name
clust_method <- "Louvain"

# read in the dataframe
df <- read.csv(input_path)
# print info to log 
print(dim(df)) # this is specific df has 73578 cells
# the preprocessing output csv needs to be cleaned - it contains live dead, FSC, SSC and the sample column
print(colnames(df))
# create a df with just the expression 
# need a way to automate this selection 
# I only want the expression values
df2 <- df %>% select(c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71"))
# the order of the DF is set by the order the columns are written above
# create a matrix for later
m <- as.matrix(df2) 

# create the seurat object for visualization

tm <- t(df2)
rownames(tm) <- colnames(df2)
colnames(tm) <- rownames(df2)
seu <- CreateSeuratObject(tm) # create a seurat object 

# add the meta data back in for sample groups
seu <- AddMetaData(object=seu, metadata=df$Batch, col.name = 'Batch')
# this doesn't work for making levels
# create the vector for the antibodies names for feature plotting later
AB <- colnames(df2)
# add to scale data slot
seu <- ScaleData(seu)

# check the data
pdf(paste(output_path,input_name,clust_method,"Heatmap_batch.pdf",sep=""),width =8, height = 6)
print(DoHeatmap(seu, group.by = "Batch", features = AB))
dev.off()

# create the UMAP
seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)


############################## explore parameters and calculate statistics ###########################



############################# define pararmater range and set up list for stats ########################################

#shuming: im getting NaN for all clusters with res = 0.01
#those clusters seem to have level 0? 
kn = c(50,100,120,140,160,200,250,300)
resolutions = c(0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5, 0.8, 1.2)



# not in the aligned transformed the number of clusters is very high at low k and higher
# more clusters are being formed in all methods


#create 3 df for stats, can be simplified later
si <- data.frame(matrix(ncol = length(resolutions), nrow = length(kn)))
colnames(si) <- resolutions
rownames(si) <- kn

ch <- data.frame(matrix(ncol = length(resolutions), nrow = length(kn)))
colnames(ch) <- resolutions
rownames(ch) <- kn

db <- data.frame(matrix(ncol = length(resolutions), nrow = length(kn)))
colnames(db) <- resolutions
rownames(db) <- kn

# track number of clusters
nc <- data.frame(matrix(ncol = length(resolutions), nrow = length(kn)))
colnames(db) <- resolutions
rownames(db) <- kn

#subsampling for silhouette score, n=1000, can make n bigger if needed
row_n <- sample(1:nrow(m), 1000)
dis <- dist(m[row_n,])


# In the loop 
# save a data object for each kn - will only keep temporarily
# the clusters will write over with each new kn


for (i in kn){
  seu <- FindNeighbors(seu, dims = 1:12, k.param = i)
  seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
  # save feature plots of this UMAP
  # file name
  UMAP_name = paste("UMAPfeatures_kn",i,".pdf",sep="")
  # save feature plots UMAP
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""))
  print(FeaturePlot(seu, features = AB,slot = 'scale.data',min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ theme(plot.title = element_text(size = 0.1)))
  dev.off()
  
  # look at batches
  UMAP_name = paste("UMAPbatches_kn",i,".pdf",sep="")
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""))
  print(DimPlot(seu,group.by = 'Batch',label.size = 1))
  dev.off()
  
  for (j in resolutions) {
    seu <- FindClusters(seu, resolution = j)
    louvainCluster <- seu@meta.data$seurat_clusters
    numb.clusters = unique(seu@meta.data$seurat_clusters)
    
    #shuming: get an error whenever i run stats on level 0 clusters (mentioned in 
    #line 94) this line skips the rest of the loop when that happens, you might  
    #want to move stats to the end of the inner loop so you don't skip the plots 
    # make plots
    # UMAP
    UMAP_name = paste("UMAPclusters_kn",i,"_res_",j,".pdf",sep="")
    print(UMAP_name) #testing 
    pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""))
    # save UMAP grouped
    print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE)) # will automatically group by active ident
    dev.off()
    # heatmap
    heatmap_name = paste("Heatmapclusters_kn",i,"_res_",j,".pdf",sep="")
    #testing 
    pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""))
    print(DoHeatmap(seu, features = AB))
    dev.off()
    
    dotplot_name = paste("Dotplotclusters_kn",i,"_res_",j,".pdf",sep="")
    #testing 
    pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""))
    print(DotPlot(seu, features = AB, group.by = 'seurat_clusters'))
    dev.off()
    
    # this should add the number of clusters into a list - it didn't seem to work, something like this should work
    nc[as.character(i), as.character(j)] <- length(unique(louvainCluster))
    
    #if (length(unique(louvainCluster))==1) next
  
    #silhouette score:
    si[as.character(i), as.character(j)] <- mean(silhouette(as.numeric(louvainCluster[row_n]),dis)[, 3])
    length(as.numeric(dis))
    silhouette(as.numeric(louvainCluster[row_n]),dis)
    # @shuming - what is the point of the next line
    as.numeric(louvainCluster[row_n])
    #Calinski-Harabasz index: 
    ch[as.character(i), as.character(j)] <- calinhara(m,louvainCluster,cn=i)
    
    # Davies–Bouldin index:
    db[as.character(i), as.character(j)] <- index.DB(df2, as.numeric(louvainCluster))$DB
    
    # save stats for each resolution
    # write.csv(stats_list, paste(output_path,list_name,sep=""))
  }
  
  # run clustree
  pdf(paste(output_path,input_name,clust_method,"kn",i,'Clustree.pdf',sep=""),width =15, height = 10)
  print(clustree(seu, prefix ='RNA_snn_res.'))
  dev.off()
  
  # save seurat object
  seu_name = paste("SeuratObject",i,".Rds",sep="")
  saveRDS(seu, paste(output_path,input_name,clust_method,seu_name,sep=""))
  # make clustree plot
  # save all stats outputs for each kn
}


# save the stats 
stats_list <- list(si,ch,db,nc)

saveRDS(stats_list,paste(output_path,input_name,clust_method,'statslist.Rds',sep=""))


############################## Plot outputs ########################################

#Shuming: this could also be simplified later by making a general plotting function

#silhouette score: ranges from -1  to 1 
#-1: bad clusters  0: neutral, indifferent  1: good clusters

pdf(paste(output_path,input_name,clust_method,"Silhouetteplot.pdf",sep=""), width = 4, height = 4)
si_new <- cbind(kn = rownames(si), si)
melted <- melt(si_new,  id.vars = 'kn', variable.name = 'resolutions')
ggplot(melted, aes(kn, value)) + 
  geom_line(aes(colour = resolutions, group = resolutions)) + 
  labs(title = "Silhouette Scores", x = "kn", y = "Average Silhouette Scores") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

write.csv(si_new,paste(output_path,input_name,clust_method,"SilhouetteStats.csv",sep=""))

# #Calinski-Harabasz index: 
# # the highest value is the optimal number of clusters

pdf(paste(output_path,input_name,clust_method,"CHIplot.pdf",sep=""), width = 4, height = 4)
ch_new <- cbind(kn = rownames(ch), ch)
melted <- melt(ch_new,  id.vars = 'kn', variable.name = 'resolutions')
ggplot(melted, aes(kn, value)) + 
  geom_line(aes(colour = resolutions, group = resolutions)) + 
  labs(title = "Calinski-Harabasz Index", x = "kn", y = "Calinski-Harabasz Index") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

write.csv(ch_new,paste(output_path,input_name,clust_method,"CHIStats.csv",sep=""))

# #Davies–Bouldin index: minimum score is zero
# #the lowest value is the optimal number of clusters

pdf(paste(output_path,input_name,clust_method,"DBplot.pdf",sep=""), width = 4, height = 4)
db_new <- cbind(kn = rownames(db), db)
melted <- melt(db_new,  id.vars = 'kn', variable.name = 'resolutions')
ggplot(melted, aes(kn, value)) + 
  geom_line(aes(colour = resolutions, group = resolutions)) + 
  labs(title = "Davies–Bouldin Index", x = "kn", y = "Davies–Bouldin Index") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


# @ shuming make a plot by the number of clusters on x axis




write.csv(db_new,paste(output_path,input_name,clust_method,"DBStats.csv",sep=""))
