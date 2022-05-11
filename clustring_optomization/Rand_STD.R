

install.packages('fossil')
library(fossil) #for rand index and adjusted rand index






calculate_rand_and_std <- function(resolutions, 
                                   clustering){
  #####calculate Rand index for clustering per resolution#####
  res_ls <- paste("RNA_snn_res.", resolutions, sep = "") #list of resolutions
  
  nc_ls <- vector() #a list of number of clusters
  for (i in res_ls) {
    # print(i)
    nc <- length(table(unique(clustering@meta.data[1,i]))) #number of cluster
    nc_ls <- append(nc_ls, nc)
    #cluster labels for each cell
    # print(clustering@meta.data[i])
  }
  std <- sd(nc_ls) #standard deviation of the number of clusters
  
  
  #####standard deviation of the number of clusters#####
  
  #subsample to 1000. why? because of the exhausted memory, same problem as for silhouette
  row_n <- sample(1:length(clustering@meta.data$orig.ident), 1000) #testing
  
  
  #stats df with rand index and adjusted rand index
  stats_df2 <- data.frame(matrix(ncol = 4))
  colnames(stats_df2) <- c("clustering1", "clustering2", "ri", "ari")
  count <- 0 #row number, update every iteration in the for loop
  
  for (i in res_ls) {
    for (j in res_ls) {
      if (j<=i) next 
      count <- count+1
      ri <- rand.index(as.numeric(clustering@meta.data[row_n,i]), as.numeric(clustering@meta.data[row_n,j]))
      ari <- adj.rand.index(as.numeric(clustering@meta.data[row_n,i]), as.numeric(clustering@meta.data[row_n,j]))
      
      stats_df2[count, "clustering1"] <- i
      stats_df2[count, "clustering2"] <- j
      stats_df2[count, "ri"] <- ri
      stats_df2[count, "ari"] <- ari 
    }
  }
  #return stats_df2 and std, or save 
}


#####test#####

#25 = kn here, with multiple resolutions within
#input example:
clustering <- readRDS("/Users/shumingli/Desktop/louvain_output/FlowsetLouvainSeuratObject25.Rds")
resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 1.8)
## to call each res: 
# kn25$RNA_snn_res.0.05

calculate_rand_and_std(clustering=clustering, resolutions=resolutions)
