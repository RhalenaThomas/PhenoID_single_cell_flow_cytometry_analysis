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

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("FlowSOM")

# load libraries
library(clusterSim) 
library(FlowSOM)
library(flowCore)
library(cluster)
library(fpc)
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)
library(reshape2) #for plotting multiple lines (resolutions) on the same graph


############################ function #########################################
# input_path, output_path, input_name, cluster_method, kn, resolutions

louvain_stats_plotting <- function(stats_ls, input_path, output_path, input_name, clust_method) {
  
  # ############################# 3. Plot outputs #############################
  
  ##silhouette score: ranges from -1  to 1
  ##-1: bad clusters  0: neutral, indifferent  1: good clusters
  
  pdf(paste(output_path,input_name,clust_method,"Silhouetteplot.pdf",sep=""))
  #x axis = number of cluster
  siplot1 <- ggplot(stats_ls, aes(x=nc, y=si, label=resolution)) +
    geom_line(aes(group=kn,color=factor(kn)), size=0.15) +
    geom_text(aes(label=resolution, colour=factor(kn)), 
              check_overlap = TRUE, position=position_jitter(width=0.2), size=3) +
    labs(color = "kn", title = "Silhouette Scores", 
         x = "Number of Clusters", y = "Average Silhouette Scores") +
    theme(plot.title = element_text(hjust = 0.5)) 
  print(siplot1)
  
  #x axis = kn
  siplot2 <- ggplot(stats_ls, aes(kn, si)) + 
    geom_point(aes(colour = factor(resolution), group=factor(resolution))) + 
    geom_line(aes(colour = factor(resolution), group=factor(resolution)), size=0.2) +
    labs(title = "Silhouette Scores", x = "kn", y = "Average Silhouette Scores") +
    theme(plot.title = element_text(hjust = 0.5))
  print(siplot2)
  
  dev.off()
  
  
  ##Calinski-Harabasz index:
  ## the highest value is the optimal number of clusters
  
  #x axis = number of cluster
  pdf(paste(output_path,input_name,clust_method,"CHIplot.pdf",sep=""), width = 4, height = 4)
  chplot1 <- ggplot(stats_ls, aes(x=nc, y=ch, label=resolution)) +
    geom_line(aes(group=kn,color=factor(kn)), size=0.15) +
    geom_text(aes(label=resolution, colour=factor(kn)), 
              check_overlap = TRUE, position=position_jitter(width=0.2), size=3) +
    labs(color = "kn", title = "Calinski-Harabasz Index", 
         x = "Number of Clusters", y = "Calinski-Harabasz Index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  print(chplot1)
  
  #x axis = kn
  chplot2 <- ggplot(stats_ls, aes(kn, ch)) + 
    geom_point(aes(colour = factor(resolution), group=factor(resolution))) + 
    geom_line(aes(colour = factor(resolution), group=factor(resolution)), size=0.2) +
    labs(title = "Calinski-Harabasz Index", x = "kn", y = "Calinski-Harabasz Index") +
    theme(plot.title = element_text(hjust = 0.5))
  print(chplot2)
  
  dev.off()
  
  
  
  ## #Davies–Bouldin index: minimum score is zero
  ## #the lowest value is the optimal number of clusters
  pdf(paste(output_path,input_name,clust_method,"DBplot.pdf",sep=""), width = 4, height = 4)
  #x axis = number of cluster
  dbplot1 <- ggplot(stats_ls, aes(x=nc, y=db, label=resolution)) +
    geom_line(aes(group=kn,color=factor(kn)), size=0.15) +
    geom_text(aes(label=resolution, colour=factor(kn)), 
              check_overlap = TRUE, position=position_jitter(width=0.2), size=3) +
    labs(color = "kn", title = "Davies–Bouldin index", 
         x = "Number of Clusters", y = "Davies–Bouldin index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  print(dbplot1)
  
  #x axis = kn
  dbplot2 <- ggplot(stats_ls, aes(kn, db)) + 
    geom_point(aes(colour = factor(resolution), group=factor(resolution))) + 
    geom_line(aes(colour = factor(resolution), group=factor(resolution)), size=0.2) +
    labs(title = "Davies–Bouldin index", x = "kn", y = "Davies–Bouldin index") +
    theme(plot.title = element_text(hjust = 0.5))
  print(dbplot2)
  dev.off()
}


louvain_clustering <- function(input_path, 
                               output_path, 
                               input_name, 
                               clust_method, 
                               kn = c(25,50,100,125,150,200,250,300), #test all values by default 
                               resolutions = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1.0,1.8)) {
 
  # ========================= 1. processing input data =========================
  df <- read.csv(input_path)
  
  # create a df with just the expression
  # need a way to automate this selection
  # I only want the expression values
  df2 <- df %>% select(c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29",
                         "CD56", "O4","CD140a","CD133","GLAST","CD71"))
  
  m <- as.matrix(df2)
  
  # create the seurat object for visualization
  tm <- t(df2)
  rownames(tm) <- colnames(df2)
  colnames(tm) <- rownames(df2)
  seu <- CreateSeuratObject(tm)
  
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
  
  # create PCA 
  #from Shuming: @Rhalena is this still necessary?
  seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)
  
  # not in the aligned transformed the number of clusters is very high at low k and higher
  # more clusters are being formed in all methods
  
  
  
  # ==== 2. clustering with diff combination of kn and resolution and stats ====
  
  #subsampling for silhouette score, n=1000, can make n bigger if needed
  row_n <- sample(1:nrow(m), 1000) #testing
  dis <- dist(m[row_n,])
  
  
  #create a df to store all stats
  stats_ls <- data.frame(matrix(ncol = 6, nrow = length(length(kn)*length(resolutions))))
  colnames(stats_ls) <- c("kn", "resolution", "nc","si", "ch", "db")
  count <- 0
  
  # In the loop
  # save a data object for each kn - will only keep temporarily
  # the clusters will write over with each new kn
  for (i in kn){
    seu <- FindNeighbors(seu, dims = 1:12, k.param = i)
    seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
    
    # save feature plots of this UMAP
    UMAP_name = paste("UMAPfeatures_kn",i,".pdf",sep="")  # file name
    pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
    print(FeaturePlot(seu, features = AB,slot = 'scale.data',min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ theme(plot.title = element_text(size = 0.1)))
    dev.off()
    
    # look at batches
    UMAP_name = paste("UMAPbatches_kn",i,".pdf",sep="")
    pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
    print(DimPlot(seu,group.by = 'Batch',label.size = 1))
    dev.off()
    
    for (j in resolutions) {
      seu <- FindClusters(seu, resolution = j)
      louvainCluster <- seu@meta.data$seurat_clusters
      numb.clusters = unique(seu@meta.data$seurat_clusters)
      
      if (length(unique(louvainCluster))==1) next #skip the ones with only 1 cluster
      
      count <- 1+count
      print(count)
      colnames(stats_ls) <- c("kn", "resolution", "nc","si", "ch", "db")
      stats_ls[count,"kn"] <- i
      stats_ls[count,"resolution"] <- j
      
      
      # number of clusters
      stats_ls[count,"nc"] <- length(unique(louvainCluster))
      
      # silhouette score:
      stats_ls[count,"si"] <- mean(silhouette(as.numeric(louvainCluster[row_n]),dis)[, 3])
      
      # Calinski-Harabasz index:
      stats_ls[count,"ch"] <- calinhara(m,louvainCluster,cn=i)
      
      # Davies–Bouldin index:
      stats_ls[count,"db"] <- index.DB(df2, as.numeric(louvainCluster))$DB
      
      
      # make plots
      # UMAP
      UMAP_name = paste("UMAPclusters_kn",i,"_res_",j,".pdf",sep="")
      print(UMAP_name) #testing
      pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =15, height = 10)
      # save UMAP grouped
      print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE)) # will automatically group by active ident
      dev.off()
      
      
      # heatmap
      heatmap_name = paste("Heatmapclusters_kn",i,"_res_",j,".pdf",sep="")
      print(UMAP_name) #testing
      pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""),width =15, height = 10)
      print(DoHeatmap(seu, features = AB))
      dev.off()
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
  
  #testing
  write.csv(stats_ls, paste(output_path, "stats.csv",sep=""), row.names = FALSE)
  
  # save the stats
  saveRDS(stats_ls,paste(output_path,input_name,clust_method,'statslist.Rds',sep=""))
  
  #stats plot
  louvain_stats_plotting(stats_ls, output_path, input_name, clust_method)
}



input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsaligned_transformed_flowset.csv"
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/cluster_parameters/Louvain/"
# input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
# output_path <- "/Users/shumingli/Desktop/"


input_name <- "Flowset"  # processing type for file name
clust_method <- "Louvain" # cluster type for file name

louvain_clustering(input_path, output_path, input_name, clust_method)


#next step: make saving pdf optional


