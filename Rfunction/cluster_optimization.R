
rm(list=ls()) 
library(clusterSim) #new package for dbi
library(FlowSOM)
library(flowCore)
library(cluster) #for silhouette score
library(fpc) #for calinhara
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree) #for clustree plot
library(Rphenograph)
library(reshape2) #for plotting multiple lines (resolutions) on the same graph

library(kit) # for finding max and second max (function topn)
library(tidyr) #for the last plot in the function


# ==================  I. louvain stats plot function starts  ====================

#1. function, embedded in louvain function, change this to stats plot function for all

# clust_method = one of c("louvain", "flowsom", "phenograph")
stats_plot <- function(stats_ls,
                       output_path, 
                       clust_method,
                       input_name = NULL) {
  
  # drop any rows containing NAs, they contain NAs because some kn x res 
  #give number of clusters of 1 (there is no split), and you can't run 
  #internal stats on them 
  stats_ls <- stats_ls[complete.cases(stats_ls), ]
  
  if (clust_method == "louvain") {
    ##silhouette score: ranges from -1  to 1
    ##-1: bad clusters  0: neutral, indifferent  1: good clusters
    
    #x axis = number of cluster
    siplot1 <- ggplot(stats_ls, aes(x = nc, y = si, label = resolution)) +
      geom_line(aes(group = kn, color = factor(kn)), size = 0.15) +
      geom_text(aes(label = resolution, colour = factor(kn)), 
                check_overlap = TRUE, 
                position = position_jitter(width = 0.2), 
                size = 3) +
      labs(color = "kn", title = "Silhouette Scores", 
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    #x axis = kn
    siplot2 <- ggplot(stats_ls, aes(kn, si)) + 
      geom_point(aes(colour = factor(resolution), group = (resolution))) + 
      geom_line(aes(colour = factor(resolution), group = (resolution)), size = 0.2) +
      labs(title = "Silhouette Scores", 
           x = "kn", 
           y = "Average Silhouette Scores", 
           colour = 'Resolution') +
      theme(plot.title = element_text(hjust = 0.5))
    
    ##Calinski-Harabasz index:
    ## the highest value is the optimal number of clusters
    
    #x axis = number of cluster
    chplot1 <- ggplot(stats_ls, aes(x = nc, y = ch, label = resolution)) +
      geom_line(aes(group = kn,color = factor(kn)), size = 0.15) +
      geom_text(aes(label = resolution, colour = factor(kn)), 
                check_overlap = TRUE, 
                position=position_jitter(width = 0.2), 
                size = 3) +
      labs(color = "kn", 
           title = "Calinski-Harabasz Index", 
           x = "Number of Clusters", 
           y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    #x axis = kn
    chplot2 <- ggplot(stats_ls, aes(kn, ch)) + 
      geom_point(aes(colour = factor(resolution), group = factor(resolution))) + 
      geom_line(aes(colour = factor(resolution), group = factor(resolution)), size=0.2) +
      labs(title = "Calinski-Harabasz Index", 
           x = "kn", 
           y = "Calinski-Harabasz Index",
           colour = 'Resolution') +
      theme(plot.title = element_text(hjust = 0.5))
    
    
    ## #Davies–Bouldin index: minimum score is zero
    ## #the lowest value is the optimal number of clusters
    
    #x axis = number of cluster
    dbplot1 <- ggplot(stats_ls, 
                      aes(x = nc, y = db, label = resolution)) +
      geom_line(aes(group = kn,color = factor(kn)), size = 0.15) +
      geom_text(aes(label = resolution, colour = factor(kn)), 
                check_overlap = TRUE, 
                position = position_jitter(width = 0.2), size = 3) +
      labs(color = "kn", title = "Davies-Bouldin index", 
           x = "Number of Clusters", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    #x axis = kn
    dbplot2 <- ggplot(stats_ls, aes(kn, db)) + 
      geom_point(aes(colour = factor(resolution), group = factor(resolution))) + 
      geom_line(aes(colour = factor(resolution), group = factor(resolution)), 
                size = 0.2) +
      labs(title = "Davies-Bouldin index", 
           x = "kn", 
           y = "Davies-Bouldin index", 
           colour = 'Resolution') +
      theme(plot.title = element_text(hjust = 0.5))

  } else if (clust_method == "flowsom") {
    
    siplot1 <- ggplot(stats_ls) + 
      geom_point(aes(x = krange, y = si)) +
      geom_line(aes(x = krange, y = si), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "krange", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    siplot2 <- ggplot(stats_ls) + 
      geom_point(aes(x = nc, y = si)) +
      geom_line(aes(x = nc, y = si), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    chplot1 <- ggplot(stats_ls) + 
      geom_point(aes(x = krange, y = ch)) +
      geom_line(aes(x = krange, y = ch), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "krange", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    chplot2 <- ggplot(stats_ls) + 
      geom_point(aes(x = nc, y = ch)) +
      geom_line(aes(x = nc, y = ch), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "Number of Clusters", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    dbplot1 <- ggplot(stats_ls) + 
      geom_point(aes(x = krange, y = db)) +
      geom_line(aes(x = krange, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "krange", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    dbplot2 <- ggplot(stats_ls) + 
      geom_point(aes(x = nc, y = db)) +
      geom_line(aes(x = nc, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "Number of Clusters", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
  } else if (clust_method == "phenograph") {
    siplot1 <- ggplot(stats_ls) + 
      geom_point(aes(x = kn, y = si)) +
      geom_line(aes(x = kn, y = si), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "kn", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    siplot2 <- ggplot(stats_ls) + 
      geom_point(aes(x = nc, y = si)) +
      geom_line(aes(x = nc, y = si), size = 0.1) +
      labs(title = "Silhouette Scores",
           x = "Number of Clusters", y = "Average Silhouette Scores") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    chplot1 <- ggplot(stats_ls) + 
      geom_point(aes(x = kn, y = ch)) +
      geom_line(aes(x = kn, y = ch), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "kn", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    chplot2 <- ggplot(stats_ls) + 
      geom_point(aes(x = nc, y = ch)) +
      geom_line(aes(x = nc, y = ch), size = 0.1) +
      labs(title = "Calinski-Harabasz Index",
           x = "Number of Clusters", y = "Calinski-Harabasz Index") +
      theme(plot.title = element_text(hjust = 0.5))
   
    dbplot1 <- ggplot(stats_ls) + 
      geom_point(aes(x = kn, y = db)) +
      geom_line(aes(x = kn, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "kn", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5)) 
    
    dbplot2 <- ggplot(stats_ls) + 
      geom_point(aes(x = nc, y = db)) +
      geom_line(aes(x = nc, y = db), size = 0.1) +
      labs(title = "Davies-Bouldin index",
           x = "Number of Clusters", y = "Davies-Bouldin index") +
      theme(plot.title = element_text(hjust = 0.5))  
  } else (
    warning("clust_method is incorrect. ")
  )
  
  pdf(paste(output_path, input_name, clust_method, "_stats_plot.pdf", sep = ""))
  print(siplot1)
  print(siplot2)
  print(chplot1)
  print(chplot2)
  print(dbplot1)
  print(dbplot2)
  dev.off()
  
  return(list(siplot1, siplot2, chplot1, chplot2, dbplot1, dbplot2))
}


# 2. input
stats_ls <- sdf
clust_method <- "louvain"
output_path <- "/Users/shumingli/Desktop/output_jul7/"


# 3. test

sl <- stats_plot(stats_ls, output_path, clust_method, input_name)

# ================== louvain stats plot function ends  ====================


# ========================  II. louvain function starts  ===============================

#1. louvain function: 

louvain_clustering <- function(x,
                               output_path,
                               input_name = NULL,
                               clust_method,
                               kn = c(25,50,100,125,150,200,250,300),
                               resolutions = c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1.0,1.8)) {
  
  #1. preprocessing
  dir.create(file.path(output_path), showWarnings = FALSE)
  setwd(file.path(output_path))
  
  # create a df with just the expression
  # need a way to automate this selection
  # I only want the expression values
  df2 <- df %>% dplyr::select(c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29",
                         "CD56", "O4","CD140a","CD133","GLAST","CD71"))
  
  # create the seurat object for visualization
  tm <- t(df2)
  rownames(tm) <- colnames(df2)
  colnames(tm) <- df$X
  seu <- CreateSeuratObject(tm)
  
  # add the meta data back in for sample groups
  # this doesn't work for making levels
  seu <- AddMetaData(object=seu, metadata=df$Batch, col.name = 'Batch')
  AB <- colnames(df2) # save antibody names for feature plotting later
  seu <- ScaleData(seu) # add to scale data slot
  
  
  # check the data
  pdf(paste(input_name, clust_method,"Heatmap_batch.pdf",sep=""),width =8, height = 6)
  print(DoHeatmap(seu, group.by = "Batch", features = AB))
  dev.off()
  
  # create PCA 
  seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)
  # not in the aligned transformed the number of clusters is very high at low k and higher
  # more clusters are being formed in all methods
  
  #subsampling for silhouette score, n=1000, can make n bigger if needed
  m <- as.matrix(df2)
  row_n <- sample(1:nrow(m), ifelse(nrow(m) > 10000, 10000, nrow(m))) #testing
  dis <- dist(m[row_n,])
  
  #create a list to store all stats
  statsl <- vector()
  
  #2. loop
  # In the loop
  # save a data object for each kn - will only keep temporarily
  # the clusters will write over with each new kn
  for (i in kn){
    seu <- FindNeighbors(seu, dims = 1:12, k.param = i)
    seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
    
    # save feature plots of this UMAP
    pdf(paste(input_name, 
              clust_method, 
              "UMAPfeatures_kn", i, ".pdf", sep = ""), 
        width =20, 
        height = 10)
    print(FeaturePlot(seu, 
                      features = AB,
                      slot = 'scale.data',
                      min.cutoff = 'q1', 
                      max.cutoff ='99',
                      label.size = 1)+ 
            theme(plot.title = element_text(size = 0.1)))
    dev.off()
    
    # look at batches
    pdf(paste(output_path,
              input_name,
              clust_method,
              "UMAPbatches_kn", i, ".pdf", sep=""),
        width =20, 
        height = 10)
    print(DimPlot(seu, group.by = 'Batch', label.size = 1))
    dev.off()
    
    for (j in resolutions) {
      seu <- FindClusters(seu, resolution = j)
      louvainCluster <- seu@meta.data$seurat_clusters
      numb.clusters = unique(seu@meta.data$seurat_clusters)
      
      #stats
      statsl <- c(statsl,
                  i, # kn
                  j, # resolution
                  length(unique(louvainCluster))) # number of clusters (nc
      
      if (length(unique(louvainCluster)) == 1)  {#skip the ones with only 1 cluster
        statsl <- c(statsl, rep(NA, 3))
        } else { statsl <- c(statsl,
                            mean(silhouette(as.numeric(louvainCluster[row_n]),dis)[, 3]), # silhouette score
                            calinhara(m,louvainCluster,cn=i),# Calinski-Harabasz index
                            index.DB(df2, as.numeric(louvainCluster))$DB) # Davies–Bouldin index
        
        # make UMAP grouped plots
        pdf(paste(input_name,
                  clust_method,
                  "UMAPclusters_kn", i, "_res_", j, ".pdf", sep = ""), 
            width = 15, height = 10)
        print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE)) # will automatically group by active ident
        dev.off()
        
        # heatmap
        pdf(paste(input_name,
                  clust_method,
                  "Heatmapclusters_kn", i, "_res_", j, ".pdf", sep = ""),
            width = 15, height = 10)
        print(DoHeatmap(seu, features = AB, size = 10)+
                theme(text = element_text(size = 30))) 
        dev.off()
        }
    }
      
    if (length(resolutions) > 1) {
      # run clustree
      pdf(paste(input_name,
                clust_method,
                "kn", i, 'Clustree.pdf', sep = ""),
          width = 15, height = 10)
      print(clustree(seu, prefix ='RNA_snn_res.'))
      dev.off()
    }
    
    
    # save seurat object
    saveRDS(seu, paste(input_name,
                       clust_method,
                       "SeuratObject", i, ".Rds", sep = ""))
    
    saveRDS(statsl, 'statsl.Rds') #for testing
  }
  
  stats_ls <- data.frame(matrix(statsl, ncol=6, byrow = TRUE))
  colnames(stats_ls) <- c("kn", "resolution", "nc","si", "ch", "db")
  
  # write.csv(stats_ls, paste("stats.csv", sep = ""), row.names = FALSE)
  #put it here so that if the runtime is interrupted, the most recent list is still saved
  
  # save the stats
  saveRDS(stats_ls, paste(input_name, clust_method, 'statslist.Rds', sep = ""))
  
  # stats plot
  
  sp <- stats_plot(stats_ls, output_path, clust_method, input_name)
  
  return(list(stats_ls, sp))
}

# 2. input:
input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/Old/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
output_path <- "/Users/shumingli/Desktop/output_jul7/"
df <- read.csv(input_path)
clust_method <- 'louvain'
input_name <- NULL
kn <- c(25,50,100,125,150,200,250,300)
resolutions <- c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1.0,1.8)
  
df <- df[sample(1:nrow(df), 3000),]
# kn <- c(25,50,100)
# resolutions <- c(0.05,0.1,0.2)


# 3. test:
sdf <- louvain_clustering(x = df, 
                          output_path = output_path, 
                          kn = kn,
                          resolutions = resolutions,
                          clust_method = clust_method)

# ========================  louvain function ends  ===============================


# ========================  III. flowsom function starts  ===============================

# 2. input:
input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/Old/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
output_path <- "/Users/shumingli/Desktop/output_jul7/"

df <- read.csv(input_path)
df <- df[sample(1:nrow(m), 5000),]

clust_method <- "flowsom"
krange <- 3:30




# kn <- c(25,50,100)
# resolutions <- c(0.05,0.1,0.2)

input_name <- NUL
clust_method <- 'louvain'

# 3. test:



# #1. function
# flowsom_clustering <- function(krange,
#                                input_path, 
#                                output_path, 
#                                input_name = NULL, 
#                                clust_method) {
#   #set director folder
#   dir.create(file.path(output_path), showWarnings = FALSE)
#   setwd(file.path(output_path))
#   
#   # read in the dataframe
#   df <- read.csv(input_path)
#   df2 <- df %>% select(c("AQP4", "CD24", "CD44", "CD184", "CD15", 
#                          "HepaCAM", "CD29", "CD56", "O4", "CD140a", 
#                          "CD133", "GLAST", "CD71"))
#   # the order of the DF is set by the order the colunms are written above
#   
#   m <- as.matrix(df2)  # create a matrix for later
#   
#   # create the flowframe
#   # if reading in a csv convert to flowset
#   frame <- new("flowFrame", exprs = m) #convert input to flowframe
#   fs <- ReadInput(frame) #convert flowframe to flowsom object
#   fs <- BuildSOM(fs) # build flowSOM object, no need for -1 because I cleaned the df about before making flowset 
#   fs <- BuildMST(fs) # build minimum spanning tree 
#   # BuildMST(flowSOM object generated by buildSOM)
#   
#   # create the seurat object for visualization
#   tm <- t(df2)
#   rownames(tm) <- colnames(df2)
#   colnames(tm) <- rownames(df2)
#   seu <- CreateSeuratObject(tm) # create a seurat object 
#   
#   # add the meta data back in for sample groups
#   seu <- AddMetaData(object=seu, metadata=df$Batch, col.name = 'Batch')
#   # this doesn't work for making levels
#   # create the vector for the antibodies names for feature plotting later
#   AB <- colnames(df2)
#   # add to scale data slot
#   seu <- ScaleData(seu)
#   
#   # check the data
#   pdf(paste(input_name, clust_method, "Heatmap_batch.pdf", sep = ""),
#       width = 8, height = 6)
#   print(DoHeatmap(seu, group.by = "Batch", features = AB))
#   dev.off()
#   
#   # create the UMAP
#   seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)
#   
#   # I've tried scaling the kn with the k but the values to no result in UMAP that spatial match cluster
#   # I'll just run the UMAP once with the kn = square root of the number of inputs
#   
#   kn = round(sqrt(dim(df2)[1]))
#   seu <- FindNeighbors(seu, dims = 1:12, k = kn)
#   seu <- RunUMAP(seu, dims = 1:12, n.neighbors = kn)
#   # save feature plots of this UMAP
#   # just for testing print
# 
#   # save feature plots UMAP
#   pdf(paste(input_name, 
#             clust_method, 
#             "UMAPfeatures_kn", kn, ".pdf", sep = ""), width =20, height = 10)
#   print(FeaturePlot(seu, features = AB, 
#                     slot = 'scale.data', 
#                     min.cutoff = 'q1', 
#                     max.cutoff ='99', 
#                     label.size = 1)+ 
#           theme(plot.title = element_text(size = 0.1)))
#   dev.off()
#   
#   # here k is the number of clusters
#   #shuming: somehow 2 doesn't work with flowsom, im not suring why
#   #krange = 3:30 
#   #krange = seq(from = 5, to = 100, by = 5)
#   # this cause a problem - it doesn't include the last 2 values in the loop
#   
#   # the k will be the max k for the metaclustering clustering. 
#   # save a data object for each kn - will only keep temporarily
#   # the clusters will write over with each new kn
#   
#   
#   #create a list to store all stats
#   statsl <- vector()
#   
#   #subsample for silhouette score
#   #shuming: here im using 1000 so it's not too slow, but 30,000 would have a more precise representation 
#   row_n <- sample(1:nrow(m), ifelse(nrow(m) > 10000, 10000, nrow(m))) #testing
#   dis <- dist(m[row_n, ])
#   
#   
#   for (i in krange){
#     # K max number of clusters not the kn input
#     ## run flowSOM clustering
#     ## easy flowsom method : scales data nClus is the k we are forcing
#     # fs <- FlowSOM(
#     #   frame,
#     #   nClus = i,
#     #   seed = 42
#     # )
#     # get the clusters from FlowSom
#     flowSOMcluster <- metaClustering_consensus(fs$map$codes,k = i,seed=42)
#     
#     
#     # name the clustering
#     clust_name = paste('FlowSom.k.',i,sep="")
#     # add the cluster ID into seurat object to visualize
#     seu <- AddMetaData(object=seu, metadata= flowSOMcluster[fs$map$mapping[,1]], col.name = clust_name)
#     number.clusters <- length(unique(flowSOMcluster[fs$map$mapping[,1]]))
#     
#     
#     ### make umap
#     #UMAP_name = paste("UMAPclusters_k",i,".pdf",sep="")
#     #print(UMAP_name) #testing
#     #pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =8, height = 5)
#     # save UMAP grouped
#     #print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
#     #dev.off()
#     
#     # the pdf don't work well for a quick figure
#     UMAP_name = paste("UMAPclusters_k",i,".png",sep="")
#     print(UMAP_name) #testing
#     png(paste(output_path,input_name,clust_method,UMAP_name,sep=""))
#     # save UMAP grouped
#     print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
#     dev.off()
#     
#     # heatmap
#     #heatmap_name = paste("Heatmapclusters_k",i,".pdf",sep="")
#     #testing
#     #pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""),width =8, height = 5)
#     #print(DoHeatmap(seu, features = AB,group.by = clust_name))
#     #dev.off()
#     # heatmap
#     heatmap_name = paste("Heatmapclusters_k",i,".png",sep="")
#     #testing
#     png(paste(output_path,input_name,clust_method,heatmap_name,sep=""), width = 600, height = 500)
#     print(DoHeatmap(seu, features = AB,group.by = clust_name))
#     dev.off()
#     
#     #### add stats
#     # calculate the statistics
#     
#     #number of clusters 
#     # "krange", "nc","si", "ch", "db"
#     count <- 1+count
#     
#     stats_ls[count, "krange"] <- i 
#     
#     #number of clusters:
#     stats_ls[count, "nc"] <- number.clusters # calculated above
#     
#     #silhouette score:
#     stats_ls[count, "si"] <- mean(silhouette(flowSOMcluster[fs$map$mapping[,1]][row_n],dis)[, 3])
#     
#     #Calinski-Harabasz index: 
#     stats_ls[count, "ch"] <- calinhara(m,flowSOMcluster[fs$map$mapping[,1]],cn=i)
#     
#     # Davies–Bouldin index:
#     stats_ls[count, "db"] <- index.DB(df2, as.numeric(flowSOMcluster[fs$map$mapping[,1]]))$DB
#     
#   }
#   saveRDS(stats_ls,paste(output_path,input_name,clust_method,'statslist.Rds',sep=""))
#   flowsom_stats_plot(stats_ls, output_path, input_name, clust_method)
#   
#   # make clustree plot
#   pdf(paste(output_path,input_name,clust_method,'Clustree.pdf',sep=""),width =8, height = 8)
#   print(clustree(seu, prefix ='FlowSom.k.'))
#   dev.off()
#   
#   # save the UMAP with cell types
#   
#   pdf(paste(output_path,input_name,clust_method,'UMAPcelltype.pdf',sep=""),width =8, height = 6)
#   print(DimPlot(seu,group.by = 'Batch'))
#   dev.off()
#   
#   # save the Seurat object
#   saveRDS(seu,paste(output_path,input_name,clust_method,'SeuratObject.Rds',sep=""))
#   
#   # save the stats list
#   write.csv(stats_ls, paste(output_path, "stats.csv",sep=""), row.names = FALSE)
#   
# }
# 
# 
# 
# # ========================  flowsom function ends  ===============================
# 
# 
# # ========================  IV. phenograph function starts  ===============================
# 
# phenograph_clustering <- function(kn, 
#                                   input_path, output_path, input_name, clust_method) {
#   # read in the dataframe
#   df <- read.csv(input_path)
#   # print info to log 
#   print(dim(df)) # this is specific df has 73578 cells
#   # the preprocessing output csv needs to be cleaned - it contains live dead, FSC, SSC and the sample column
#   print(colnames(df))
#   # create a df with just the expression 
#   # need a way to automate this selection 
#   # I only want the expression values
#   df2 <- df %>% select(c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71"))
#   # the order of the DF is set by the order the columns are written above
#   # create a matrix for later
#   m <- as.matrix(df2) 
#   
#   # create the seurat object for visualization
#   
#   tm <- t(df2)
#   rownames(tm) <- colnames(df2)
#   colnames(tm) <- rownames(df2)
#   seu <- CreateSeuratObject(tm) # create a seurat object 
#   
#   # add the meta data back in for sample groups
#   seu <- AddMetaData(object=seu, metadata=df$Batch, col.name = 'Batch')
#   # this doesn't work for making levels
#   # create the vector for the antibodies names for feature plotting later
#   AB <- colnames(df2)
#   # add to scale data slot
#   print("test1")
#   seu <- ScaleData(seu)
#   print("test2")
#   # check the data
#   pdf(paste(output_path,input_name,clust_method,"Heatmap_batch.pdf",sep=""),width =8, height = 6)
#   print(DoHeatmap(seu, group.by = "Batch", features = AB))
#   dev.off()
#   print("test3")
#   # create the UMAP
#   seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)
#   
#   print("test4")
#   # like FlowSOM Phenograph doesn't relate directly to the UMAP like Louvain
#   # we will make on seurat UMAP and visualize the clusters there
#   
#   
#   kn_umap = round(sqrt(dim(df2)[1]))
#   seu <- FindNeighbors(seu, dims = 1:12, k.param = kn_umap)
#   print("test5")
#   seu <- RunUMAP(seu, dims = 1:12, n.neighbors = kn_umap)
#   print("test6")
#   # save feature plots of this UMAP
#   # just for testing print
#   
#   # we also only need to plot the features once
#   # file name
#   UMAP_name = paste("UMAPfeatures_kn",kn_umap,".pdf",sep="")
#   print(UMAP_name) #testing
#   
#   # save feature plots UMAP
#   pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
#   print(FeaturePlot(seu, features = AB,slot = 'scale.data',min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ theme(plot.title = element_text(size = 0.1)))
#   dev.off()
#   print("test7")
#   # we also want to see the batch on the UMAP
#   pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =8, height = 6)
#   print(DimPlot(seu, group.by = 'Batch'))
#   dev.off()
#   print("test8")
#   
#   ############################## explore parameters and calculate statistics ###########################
#   
#   
#   #create 3 lists for stats
#   
#   #create a df to store all stats, col = nc, si, ch, db
#   stats_ls <- data.frame(matrix(ncol = 5, nrow = length(kn)))
#   colnames(stats_ls) <- c("kn", "nc","si", "ch", "db")
#   count <- 0
#   
#   #subsampling for silhouette score, n=1000, can make n bigger if needed
#   set.seed(25)
#   row_n <- sample(1:nrow(m), 1000)
#   dis <- dist(m[row_n,])
#   
#   print("test9")
#   
#   ############################# loop to explore parameters ########################################
#   # kn = c(25,50,75,100,125,150,175,200,225,250,300,350,400,450,500)
#   # kn = c(25,50,75,100,125,150,175,200,225,250,275,300)
#   # larger kn fewer clusters in general but not always
#   #kn = c(50,500)
#   # save a data object for each kn - will only keep temporarily
#   # the clusters will write over with each new kn
#   
#   
#   for (i in kn){
#     
#     ### run phenograph clustering
#     Rphenograph_out_flow <- Rphenograph(m, k = i)
#     print("test10")
#     clust_name = paste('Pheno.kn.',i,sep="")
#     # add the cluster ID into seurat object to visualize
#     seu <- AddMetaData(object=seu, factor(membership(Rphenograph_out_flow[[2]])), col.name = clust_name) 
#     print("test11")
#     number.clusters <- length(unique(factor(membership(Rphenograph_out_flow[[2]]))))
#     print("test12")
#     ### make umap 
#     
#     UMAP_name = paste("UMAPclusters_kn",i,".pdf",sep="")
#     print(UMAP_name) #testing 
#     pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
#     # save UMAP grouped
#     print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
#     dev.off()
#     # heatmap
#     heatmap_name = paste("Heatmapclusters_kn",i,".pdf",sep="")
#     #testing 
#     pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""),width =25, height = 10)
#     print(DoHeatmap(seu, features = AB,group.by = clust_name))
#     dev.off()
#     print("test13")
#     #### add stats
#     
#     
#     
#     # "kn", "nc","si", "ch", "db"
#     count <- 1+count
#     
#     # get the cluster indexes 
#     phenocluster <- factor(membership(Rphenograph_out_flow[[2]]))
#     print("test14")
#     
#     stats_ls[count, "kn"] <- i 
#     
#     #number of clusters:
#     stats_ls[count, "nc"] <- number.clusters # calculated above
#     
#     #silhouette score:
#     stats_ls[count, "si"] <- mean(silhouette(as.numeric(phenocluster[row_n]),dis)[, 3])
#     
#     #Calinski-Harabasz index: 
#     stats_ls[count, "ch"] <- calinhara(m,phenocluster,cn=i)
#     
#     # Davies–Bouldin index:
#     stats_ls[count, "db"] <- index.DB(df2, as.numeric(phenocluster))$DB
#     print("test15")
#   }
#   
#   
#   # make clustree plot
#   pdf(paste(output_path,input_name,clust_method,'Clustree.pdf',sep=""),width =15, height = 10)
#   print(clustree(seu, prefix ='Pheno.kn.'))
#   dev.off()
#   print("test16")
#   # save the Seurat object
#   saveRDS(seu,paste(output_path,input_name,clust_method,'SeuratObject.Rds',sep=""))
#   
#   
#   # save the stats list
#   
#   saveRDS(stats_ls,paste(output_path,input_name,clust_method,'statslist.Rds',sep=""))
#   
#   phenograph_stats_plot(stats_ls, output_path, input_name, clust_method)
#   
#   write.csv(stats_ls, paste(output_path, "stats.csv",sep=""), row.names = FALSE)
#   print("test17")
# }
# 
# 
# # ========================  phenograph function ends  ===============================
# 
# 
# 
