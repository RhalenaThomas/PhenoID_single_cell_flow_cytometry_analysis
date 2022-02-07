# flowsom clustering
# parameter tuning, statistic and visualizations for manual annotation


# load libraries

library(FlowSOM)
library(flowCore)
library(cluster)
library(fpc)
library(clv)
library(Seurat)
library(dplyr)
library(ggplot2)

############# set up the data object for clustering ############################

# info to change for each comparison
# define the input pathway
# input pathway
input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsflowset.csv"

# output pathway
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/test/"
# add input description to ouptput files
input_name <- "Flowset"  # this will be the different processing types

# cluster type for file name
clust_method <- "Phenograph"

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
seu <- RunPCA(seu, features = AB, npcs = 25)


############################## explore parameters and calculate statistics ###########################

# stats stuff - Shuming adjust 
ch_li <-list()
dbi_li <-list()



########## function for stats if you (shuming want)





############################# loop to explore parameters ########################################


kn = c(25,50,75,100,125,150,175,200,225,250,275,300)


# save a data object for each kn - will only keep temporarily
# the clusters will write over with each new kn

kn = c(50,100)



for (i in kn){
  seu <- FindNeighbors(seu, dims = 1:12, k = i)
  seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
  # save feature plots of this UMAP
  # just for testing print
  print(FeaturePlot(seu, features = AB,slot = 'scale.data',min.cutoff = 'q1', max.cutoff ='99',label.size = 1) + theme(plot.title = element_text(size = 0.1)))
  # file name
  UMAP_name = paste("UMAPfeatures_kn",i,".pdf",sep="")
  print(UMAP_name) #testing 
  # save feature plots UMAP
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
  print(FeaturePlot(seu, features = AB,slot = 'scale.data',min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ theme(plot.title = element_text(size = 0.1)))
  dev.off()
  
  ### run phenograph clustering
  Rphenograph_out_flow <- Rphenograph(m, k = i)
  
  modularity(Rphenograph_out_flow[[2]])
  membership(Rphenograph_out_flow[[2]])
  
  # add cluster ID back into original df  - this won't work in the loop the column name needs to be the kn
  df$phenograph_cluster <- factor(membership(Rphenograph_out_flow[[2]]))
  
  clust_name = paste('Pheno.k.',i,sep="")
  # add the cluster ID into seurat object to visualize
  seu <- AddMetaData(object=seu, factor(membership(Rphenograph_out_flow[[2]])), col.name = clust_name) 

  ### make umap 
  UMAP_name = paste("UMAPclusters_kn",i,".pdf",sep="")
  print(UMAP_name) #testing 
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =8, height = 5)
  # save UMAP grouped
  print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
  dev.off()
  # heatmap
  heatmap_name = paste("Heatmapclusters_kn",i,".pdf",sep="")
  #testing 
  pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""),width =8, height = 5)
  print(DoHeatmap(seu, features = AB,group.by = clust_name))
  dev.off()
  
  #### add stats
    
    # calculate the statistics
    # silouette
    # CHI
    ch = calinhara(m,seu@meta.data$seurat_clusters,length((unique(seu@meta.data$seurat_clusters))))
    ch_li[j] <- ch
    # Davies
    
    # send stats to stats_list (or df or whatever works)
    
    # make plots
    # UMAP
    
    # save stats for each resolution
    # write.csv(stats_list, paste(output_path,list_name,sep=""))
  }
  # save df with Phenograph cluster indexes

  # make clustree plot

  # save all stats outputs for each kn




############################## Plot outputs ########################################

#silhouette score: ranges from -1  to 1 
#-1: bad clusters  0: neutral, indifferent  1: good clusters
#plot(krange, type='b', li[[1]][krange], xlab='Number of clusters', ylab='Average Silhouette Scores', frame=TRUE)

#Calinski-Harabasz index: 
# the highest value is the optimal number of clusters
pdf(paste(output_path,input_name,clust_method,"CHIplot.pdf",sep=""))
print(plot(krange, type='b',ch_li[[1]][krange], xlab='Number of clusters', ylab='Calinski-Harabasz index', frame=TRUE))
dev.off()
#Davies–Bouldin index: minimum score is zero
#the lowest value is the optimal number of clusters
pdf(paste(output_path,input_name,clust_method,"DBplot.pdf",sep=""))
print(plot(krange, type='b', dbi_li[[2]][krange], xlab='Number of clusters', ylab='Davies–Bouldin index', frame=TRUE))
dev.off()


