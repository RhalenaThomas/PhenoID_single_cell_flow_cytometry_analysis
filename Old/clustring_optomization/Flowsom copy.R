# flowsom clustering
# parameter tuning, statistic and visualizations for manual annotation


# load libraries

library(clusterSim) #new package for dbi
library(FlowSOM)
library(flowCore)
library(cluster)
library(fpc)
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)

############# set up the data object for clustering ############################

# define the input pathway
# input pathway

# input_path <- "/Users/shuming/Desktop/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsflowset.csv"
input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsflowset.csv"

# output pathway
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/test/FlowSOM/"
# add input description to ouptput files
input_name <- "Flowset"  # this will be the different processing types

# cluster type for file name
clust_method <- "FlowSOM"


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
# the order of the DF is set by the order the colunms are written above
# create a matrix for later
print(colnames(df2))
m <- as.matrix(df2) 

# create the flowframe
# if reading in a csv convert to flowset
frame <- new("flowFrame", exprs = m) #convert input to flowframe
fs <- ReadInput(frame) #convert flowframe to flowsom object
fs <- BuildSOM(fs) # build flowSOM object, no need for -1 because I cleaned the df about before making flowset 
fs <- BuildMST(fs) # build minimum spanning tree 
# BuildMST(flowSOM object generated by buildSOM)

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
seu <- RunPCA(seu, features = AB, npcs = 13)


############################## explore parameters and calculate statistics ###########################


# here k is the number of clusters
#shuming: somehow 2 doesn't work with flowsom, im not suring why
#krange = 3:30 
krange = seq(from = 3, to = 99, by = 3)

# the k will be the max k for the metaclustering clustering. 
# save a data object for each kn - will only keep temporarily
# the clusters will write over with each new kn

#stats lists
si <- list()
ch <- list()
db <- list()

#subsample for silhouette score
#shuming: here im using 1000 so it's not too slow, but 30,000 would be have better representation 
row_n <- sample(1:nrow(m), 1000)
dis <- dist(m[row_n,])



for (i in krange){
  # K is the number of clusters not the kn input
  # I'll scale the kn with the k
  kn = i*10
  print(kn)
  seu <- FindNeighbors(seu, dims = 1:12, k = kn)
  seu <- RunUMAP(seu, dims = 1:12, n.neighbors = kn)
  # save feature plots of this UMAP
  # just for testing print

  # file name
  UMAP_name = paste("UMAPfeatures_kn",kn,".pdf",sep="")
  print(UMAP_name) #testing

  # save feature plots UMAP
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
  print(FeaturePlot(seu, features = AB,slot = 'scale.data',min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ theme(plot.title = element_text(size = 0.1)))
  dev.off()
  ## run flowSOM clustering
  ## easy flowsom method : scales data nClus is the k we are forcing
  fs <- FlowSOM(
    frame,
    nClus = i,
    seed = 42
  )
  # get the clusters from FlowSom
  flowSOMcluster <- GetMetaclusters(fs, meta=metaClustering_consensus(fs$map$codes,k = i,seed=42))
  
  # name the clustering
  clust_name = paste('FlowSom.k.',i,sep="")
  # add the cluster ID into seurat object to visualize
  seu <- AddMetaData(object=seu, metadata= flowSOMcluster[fs$map$mapping[,1]], col.name = clust_name)


  ### make umap
  UMAP_name = paste("UMAPclusters_k",i,".pdf",sep="")
  print(UMAP_name) #testing
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =8, height = 5)
  # save UMAP grouped
  print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
  dev.off()
  # heatmap
  heatmap_name = paste("Heatmapclusters_k",i,".pdf",sep="")
  #testing
  pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""),width =8, height = 5)
  print(DoHeatmap(seu, features = AB,group.by = clust_name))
  dev.off()
  
  #### add stats
  
  # calculate the statistics
  
  #silhouette score:
  si[i] <- mean(silhouette(flowSOMcluster[row_n],dis)[, 3])
  
  #Calinski-Harabasz index: 
  ch[i] <- calinhara(m,flowSOMcluster,cn=i)
  
  # Davies–Bouldin index:
  db[i] <- index.DB(df2, as.numeric(flowSOMcluster))$DB
  
  # send stats to stats_list (or df or whatever works)
  
  # make plots
  # UMAP
  
  # save stats for each resolution
  # write.csv(stats_list, paste(output_path,list_name,sep=""))
}

#stats list
stats_list <- list(si, ch, db)
# write.csv(stats_list, )

#make stats plots
# shuming: commenting out because I'm not sure where you want the graphs to be saved in
# pdf(),width =, height = )

#silhouette score:
#-1: bad clusters  0: neutral, indifferent  1: good clusters
pdf(paste(output_path,input_name,clust_method,'statssilhouette.pdf',sep=""),width =4, height = 4)
print(plot(krange, type='b', stats_list[[1]][krange], xlab='k max', ylab='Average Silhouette Scores', frame=TRUE))
dev.off()

#Calinski-Harabasz index: 
# the highest value is the optimal number of clusters
pdf(paste(output_path,input_name,clust_method,'statsCalinskiHara.pdf',sep=""),width =4, height = 4)
print(plot(krange, type='b', stats_list[[2]][krange], xlab='k max', ylab='Calinski-Harabasz index', frame=TRUE))
dev.off()

#Davies–Bouldin index: minimum score is zero
#the lowest value is the optimal number of clusters
pdf(paste(output_path,input_name,clust_method,'statsDavies.pdf',sep=""),width =4, height = 4)
print(plot(krange, type='b', stats_list[[3]][krange], xlab='k max', ylab='Davies–Bouldin index', frame=TRUE))
dev.off()


# make clustree plot


pdf(paste(output_path,input_name,clust_method,'Clustree.pdf',sep=""),width =8, height = 8)
print(clustree(seu, prefix ='FlowSom.k.'))
dev.off()

# save the Seurat object
saveRDS(seu,paste(output_path,input_name,clust_method,'SeuratObject.Rds',sep=""))

# save the stats list
saveRDS(stats_list,paste(output_path,input_name,clust_method,'statslist.Rds',sep=""))


