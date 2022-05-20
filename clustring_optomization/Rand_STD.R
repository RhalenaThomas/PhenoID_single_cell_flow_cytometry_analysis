# install.packages('fossil')
library(fossil) #for rand index and adjusted rand index
library(Seurat)
library(dplyr)

##run louvain 100 times under the same settings (resolution & kn)


############testing input start############
input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsretrotransformed_flowset.csv"
output_path <- "/Users/shumingli/Desktop/output_may12_rand/"
resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 1.8)
kn = c(25,50,100,125,150,200,250,300)
# resolutions <- c(0.05, 0.1)
# kn = c(25,50)
############testing input end############


calculate_rand_index <- function(input_path,
                                 output_path,
                                 resolutions,
                                 kn,
                                 n #number of iterations
                                 ) {
  df <- read.csv(input_path)
  df2 <- df %>% select(c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29",
                         "CD56", "O4","CD140a","CD133","GLAST","CD71"))
  m <- as.matrix(df2)
  tm <- t(df2)
  rownames(tm) <- colnames(df2)
  colnames(tm) <- rownames(df2)
  seu <- CreateSeuratObject(tm)
  seu <- AddMetaData(object=seu, metadata=df$Batch, col.name = 'Batch')
  AB <- colnames(df2)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE) #change npcs?
  
  #subsample to 9000. why? because of the exhausted memory, same problem as for silhouette
  row_n <- sample(1:length(seu@meta.data$orig.ident), 9000) #testing
  
  rn_ls <- round(runif(n, min=0, max=100), 0) #random number generator
  
  
  #stats df with rand index and adjusted rand index
  rand_df <- data.frame(matrix(ncol = 6))
  colnames(rand_df) <- c("kn", "resoultion", "mean", "median","sd", "ncsd") #ncsd= standard deviation of number of clusters
  count <- 0 #row number, update every iteration in the for loop
  
  for (i in kn){
    # i <- 25 #test
    seu <- FindNeighbors(seu, dims = 1:12, k.param = i)
    seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
    for (j in resolutions) {
      # j <- 0.05 #test
      
      rand_ls <- vector() #list of ri of all 100 repeats
      nc_ls <- vector()
      for (k in 1:n) {
        seu <- FindClusters(seu, random.seed = rn_ls[k], resolution = j)
        # louvainCluster <- seu@meta.data$seurat_clusters
        if (length(unique(seu@meta.data$seurat_clusters))==1) next #skip the ones with only 1 cluster
  
        seu[[paste("repeat_", k, sep = "")]] <- Idents(object = seu)
        nc_ls <- c(nc_ls, length(unique(seu@meta.data$seurat_clusters)))
        if (k==1) {next}
        for (z in 1:(k-1)) { #z is a number between 1 to k-1, for which k is one of the iterations
          # print(z)
          ri <- rand.index(as.numeric(seu@meta.data[row_n, paste("repeat_", z, sep = "")]), as.numeric(seu@meta.data[row_n, paste("repeat_", k, sep = "")]))
          # print(paste("k:", k,"z:",z,"ri:",ri))
          rand_ls <- c(rand_ls, ri)
          # print(paste(k,z))
        }
      }
      count <- count+1
      rand_df[count, 'kn'] <- i
      rand_df[count, 'resoultion'] <- j
      rand_df[count, 'mean'] <- mean(rand_ls)
      rand_df[count, 'median'] <- median(rand_ls)
      rand_df[count, 'sd'] <- sd(rand_ls)
      rand_df[count, 'ncsd'] <- sd(nc_ls)
    }
  }
  write.csv(rand_df, paste(output_path, "randIndex.csv",sep=""), row.names = FALSE)
}

#test#
calculate_rand_index(input_path,
                     output_path,
                     resolutions,
                     kn,
                     n=100)
#test#

##############################################################################
#draft from before

  
# calculate_rand_and_std <- function(resolutions, 
#                                    clustering){
#   #####calculate Rand index for clustering per resolution#####
#   res_ls <- paste("RNA_snn_res.", resolutions, sep = "") #list of resolutions
#   
#   nc_ls <- vector() #a list of number of clusters
#   for (i in res_ls) {
#     # print(i)
#     nc <- length(table(unique(clustering@meta.data[1,i]))) #number of cluster
#     nc_ls <- append(nc_ls, nc)
#     #cluster labels for each cell
#     # print(clustering@meta.data[i])
#   }
#   std <- sd(nc_ls) #standard deviation of the number of clusters
#   
#   
#   #####standard deviation of the number of clusters#####
#   
#   #subsample to 1000. why? because of the exhausted memory, same problem as for silhouette
#   row_n <- sample(1:length(clustering@meta.data$orig.ident), 1000) #testing
#   
#   
#   #stats df with rand index and adjusted rand index
#   stats_df2 <- data.frame(matrix(ncol = 4))
#   colnames(stats_df2) <- c("clustering1", "clustering2", "ri", "ari")
#   count <- 0 #row number, update every iteration in the for loop
#   
#   for (i in res_ls) {
#     for (j in res_ls) {
#       if (j<=i) next 
#       count <- count+1
#       ri <- rand.index(as.numeric(clustering@meta.data[row_n,i]), as.numeric(clustering@meta.data[row_n,j]))
#       ari <- adj.rand.index(as.numeric(clustering@meta.data[row_n,i]), as.numeric(clustering@meta.data[row_n,j]))
#       
#       stats_df2[count, "clustering1"] <- i
#       stats_df2[count, "clustering2"] <- j
#       stats_df2[count, "ri"] <- ri
#       stats_df2[count, "ari"] <- ari 
#     }
#   }
#   #return stats_df2 and std, or save 
# }
# 
# 
# #####test#####
# 
# #25 = kn here, with multiple resolutions within
# #input example:
# clustering <- readRDS("/Users/shumingli/Desktop/louvain_output/FlowsetLouvainSeuratObject25.Rds")
# resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 1.8)
# ## to call each res: 
# # kn25$RNA_snn_res.0.05
# 
# calculate_rand_and_std(clustering=clustering, resolutions=resolutions)
