# install.packages('Seurat')
library(fossil) #for rand index and adjusted rand index
library(Seurat)
library(dplyr)

##run louvain 100 times under the same settings (resolution & kn)


############testing input start############
input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsretrotransformed_flowset.csv"
# output_path <- "/Users/shumingli/Desktop/output_may12_rand/"
output_path <- "/Users/shumingli/Desktop/output_jun19/"

# #in grumio
# input_path <- "/export02/data/prepro_outsretrotransformed_flowset.csv"
# output_path <- "/export02/data/"

# resolutions <- c(0.05)
resolutions <- c(0.1, 0.3, 0.4)
kn = c(60) #just need 60  
n <- 100
# resolutions <- c(0.05, 0.1)
# kn = c(25,50)
############testing input end############


# calculate_rand_index <- function(input_path,
#                                  output_path,
#                                  resolutions,
#                                  kn,
#                                  n #number of iterations
#                                  ) {
#   df <- read.csv(input_path)
#   df2 <- df %>% select(c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29",
#                          "CD56", "O4","CD140a","CD133","GLAST","CD71"))
#   m <- as.matrix(df2)
#   tm <- t(df2)
#   rownames(tm) <- colnames(df2)
#   colnames(tm) <- rownames(df2)
#   seu <- CreateSeuratObject(tm)
#   seu <- AddMetaData(object=seu, metadata=df$Batch, col.name = 'Batch')
#   AB <- colnames(df2)
#   seu <- ScaleData(seu)
#   seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE) #change npcs?
#   
#   #subsample to 9000. why? because of the exhausted memory, same problem as for silhouette
#   row_n <- sample(1:length(seu@meta.data$orig.ident), 9000) #testing
#   
#   rn_ls <- round(runif(n, min=0, max=100), 0) #random number generator
#   
#   
#   #stats df with rand index and adjusted rand index
#   rand_df <- data.frame(matrix(ncol = 6))
#   colnames(rand_df) <- c("kn", "resoultion", "mean", "median","sd", "ncsd") #ncsd= standard deviation of number of clusters
#   count <- 0 #row number, update every iteration in the for loop
#   
#   for (i in kn){
#     # i <- 25 #test
#     seu <- FindNeighbors(seu, dims = 1:12, k.param = i)
#     seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
#     for (j in resolutions) {
#       # j <- 0.05 #test
#       
#       rand_ls <- vector() #list of ri of all 100 repeats
#       nc_ls <- vector()
#       for (k in 1:n) {
#         seu <- FindClusters(seu, random.seed = rn_ls[k], resolution = j)
#         # louvainCluster <- seu@meta.data$seurat_clusters
#         if (length(unique(seu@meta.data$seurat_clusters))==1) next #skip the ones with only 1 cluster
#   
#         seu[[paste("repeat_", k, sep = "")]] <- Idents(object = seu)
#         nc_ls <- c(nc_ls, length(unique(seu@meta.data$seurat_clusters)))
#         if (k==1) {next}
#         for (z in 1:(k-1)) { #z is a number between 1 to k-1, for which k is one of the iterations
#           # print(z)
#           ri <- rand.index(as.numeric(seu@meta.data[row_n, paste("repeat_", z, sep = "")]), as.numeric(seu@meta.data[row_n, paste("repeat_", k, sep = "")]))
#           # print(paste("k:", k,"z:",z,"ri:",ri))
#           rand_ls <- c(rand_ls, ri)
#           # print(paste(k,z))
#           
#         }
#         
#       }
#       count <- count+1
#       rand_df[count, 'kn'] <- i
#       rand_df[count, 'resoultion'] <- j
#       rand_df[count, 'mean'] <- mean(rand_ls)
#       rand_df[count, 'median'] <- median(rand_ls)
#       rand_df[count, 'sd'] <- sd(rand_ls)
#       rand_df[count, 'ncsd'] <- sd(nc_ls)
#     }
#     write.csv(rand_df, paste(output_path, "randIndex.csv",sep=""), row.names = FALSE)
#   }
# }



# calculate_rand_index <- function(input_path,
#                                  output_path,
#                                  resolutions,
#                                  kn,
#                                  n #number of iterations
# ) {

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
  colnames(rand_df) <- c("kn", "resolution", "mean", "median","sd", "ncsd") #ncsd= standard deviation of number of clusters
  count <- 0 #row number, update every iteration in the for loop
  
  for (i in kn){
    # i <- 60 #test
    # print(paste("i:", i))
    seu <- FindNeighbors(seu, dims = 1:12, k.param = i)
    seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
    for (j in resolutions) {
      # j <- 0.05 #test
      # print(paste("i:", i, "j:", j))
      rand_ls <- vector() #list of ri of all 100 repeats
      nc_ls <- vector()
      for (k in 1:n) {
        # k <- 5
        # print(paste("i:", i, "j:", j, "k:", k))
        seu <- FindClusters(seu, random.seed = rn_ls[k], resolution = j)
        # louvainCluster <- seu@meta.data$seurat_clusters
        if (length(unique(seu@meta.data$seurat_clusters))==1) next #skip the ones with only 1 cluster
        
        seu[[paste("repeat_", k, sep = "")]] <- Idents(object = seu)
        nc_ls <- c(nc_ls, length(unique(seu@meta.data$seurat_clusters)))
        if (k==1) {next}
        for (z in 1:(k-1)) { #z is a number between 1 to k-1, for which k is one of the iterations
          # z <- 4
          # print(z)
          # print(paste("i:", i, "j:", j, "k:", k, "z:", z))
          ri <- rand.index(as.numeric(seu@meta.data[row_n, paste("repeat_", z, sep = "")]), as.numeric(seu@meta.data[row_n, paste("repeat_", k, sep = "")]))
          
          # as.numeric(seu@meta.data[row_n, paste("repeat_", z, sep = "")])
          
          # print(paste("k:", k,"z:",z,"ri:",ri))
          rand_ls <- c(rand_ls, ri)
          # print(paste(k,z))
          saveRDS(rand_ls, paste(output_path, "rand_ls_", j, ".Rds",sep=""))
          saveRDS(nc_ls, paste(output_path, "nc_ls_", j, ".Rds",sep=""))
        }
        
      }
      count <- count+1
      rand_df[count, 'kn'] <- i
      rand_df[count, 'resolution'] <- j
      rand_df[count, 'mean'] <- mean(rand_ls)
      rand_df[count, 'median'] <- median(rand_ls)
      rand_df[count, 'sd'] <- sd(rand_ls)
      rand_df[count, 'ncsd'] <- sd(nc_ls)
      write.csv(rand_df, paste(output_path, "randIndex.csv",sep=""), row.names = FALSE)
    }
    write.csv(rand_df, paste(output_path, "randIndex.csv",sep=""), row.names = FALSE)
  }
# }




# #test#
# calculate_rand_index(input_path,
#                      output_path,
#                      resolutions,
#                      kn,
#                      n=5)
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
resolutions <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1, 1.8)
# ## to call each res: 
# # kn25$RNA_snn_res.0.05
# 
# calculate_rand_and_std(clustering=clustering, resolutions=resolutions)
  
  
  
  
#=============== draw graph for rand   =============== 
#x axis = resolution, y axis = cluster number, rand index 

library(ggplot2)  
  
ridf <- read.csv('/Users/shumingli/Desktop/output_jun19/randIndex0.05.csv')  


ggplot(data=ridf, aes(x=resolution)) +
  geom_line(aes(y = mean, color = 'mean rand index'))+
  geom_point(aes(y = mean, color = 'mean rand index'))+ 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.02, 
                position=position_dodge(0.05)) +
  scale_x_continuous(breaks= resolutions)



ggplot(data=ridf, aes(x=resolution)) +
  geom_line(aes(y = ncmean, color = 'sd of number of cluster'))+
  geom_point(aes(y = ncmean, color = 'sd of number of cluster'))+
  geom_errorbar(aes(ymin=ncmean-ncsd, ymax=ncmean+ncsd), width=.02,
                position=position_dodge(0.05))+
  scale_x_continuous(breaks= resolutions)



p <- ggplot(df, aes(x=Category, y=Mean, fill=Quality)) + 
  geom_point()+
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2,
                position=position_dodge(0.05))




nc_ls0.05 <-readRDS('/Users/shumingli/Desktop/output_jun19/nc_ls/nc_ls0.05.Rds')
nc_ls0.2 <- readRDS('/Users/shumingli/Desktop/output_jun19/nc_ls/nc_ls0.2.Rds')
nc_ls_0.5 <- readRDS('/Users/shumingli/Desktop/output_jun19/nc_ls/nc_ls_0.5.Rds')
nc_ls_0.6 <- readRDS('/Users/shumingli/Desktop/output_jun19/nc_ls/nc_ls_0.6.Rds')
nc_ls_0.7 <- readRDS('/Users/shumingli/Desktop/output_jun19/nc_ls/nc_ls_0.7.Rds')
nc_ls_1 <- readRDS('/Users/shumingli/Desktop/output_jun19/nc_ls/nc_ls_1.Rds')

nc <- list(nc_ls0.05, nc_ls0.2, nc_ls_0.5, nc_ls_0.6, nc_ls_0.7, nc_ls_1)
reso <- c(0.05, 0.2, 0.5, 0.6, 0.7, 1)




# rand_std_df <- data.frame(matrix(ncol = 2))
# colnames(rand_std_df) <- c("resolution", "randsd") #ncsd= standard deviation of number of clusters
# 
# nc_std_df <- data.frame(matrix(ncol = 2))
# colnames(nc_std_df) <- c("resolution", "ncsd") #ncsd= standard deviation of number of clusters

nc_std_ls1 <- c()
nc_std_ls2 <- c()

for (i in reso) {
  # count <- 1+count
  print(i)
  # ncl <- readRDS(paste('/Users/shumingli/Desktop/output_jun19/nc_ls/nc_ls_', as.character(i), '.Rds', sep=''))
  randsd <- readRDS(paste('/Users/shumingli/Desktop/output_jun19/rand_ls/rand_ls_', as.character(i), '.Rds', sep=''))
  # readRDS(paste('/Users/shumingli/Desktop/output_jun19/nc_ls/nc_ls0.05.Rds'))
  
  
  resolution <- rep(as.numeric(i), length(randsd))

  # b <- (rbind(a, randl, randl))
  # print(cbind(a,randl))
  print('test1')
  
  # cb <- cbind(resolution, randsd)
  
  print('test2')
  # nc_std_ls <- rbind(cb,nc_std_ls)
  
  nc_std_ls1 <- c(nc_std_ls1, resolution)
  nc_std_ls2 <- c(nc_std_ls2, randsd)
  
  
  print('test3')
  # print(randsd)
  # print(resolution)
  # a <- rep(i, length(ncl))
  # b <- (rbind(a, randl, ncl))
  
  # std_df <- rbind(std_df, cbind(a, ncl, randl))
}

# nc_std_ls <- c(resolution = nc_std_ls1, mean = nc_std_ls2)
# 
# c(list(ridf$resolution), nc_std_ls1)
# ridf$resolution <- nc_std_ls1[1:9]



# new_df <- data.frame(resolution=unlist(append(list(ridf$resolution), nc_std_ls1)), mean=unlist(append(list(ridf$mean), nc_std_ls2)))

# new_df <- 

# append(list(ridf$resolution), nc_std_ls1)

# ridf <- c(ridf, nc_std_ls)

# unlist(append(list(ridf$resolution), nc_std_ls1)) 
# unique(nc_std_ls[,resolution])

# nc_std_ls[,'resolution']



i <- 0.05
randl <- readRDS(paste('/Users/shumingli/Desktop/output_jun19/rand_ls/rand_ls_', as.character(i), '.Rds', sep=''))

new_df <- data.frame(resolution=unlist(append(list(ridf$resolution), nc_std_ls1)), mean=ridf$mean, rand=c(rep(NA, length(9)), nc_std_ls2))



print('coo')

for (i in 1:length(nc)) {
  print(mean(nc[[i]]))
}
  
as.numeric(nc_std_ls)

as.data.frame(nc_std_ls)
ggplot(data=ridf, aes(x=resolution)) +
  geom_line(aes(y = mean, color = 'mean rand index'))+
  geom_point(aes(y = mean, color = 'mean rand index'))+
  geom_point(data=as.data.frame(nc_std_ls), aes(x = nc_std_ls[,1], y = nc_std_ls[,2], color = 'mean rand index'))+
  geom_line(aes(y = ncsd, color = 'sd of number of cluster'))+
  geom_point(aes(y = ncsd, color = 'sd of number of cluster')) 


ggplot(df, aes(x = day)) + 
  geom_line(aes(y = sales, color = 'sales')) + 
  geom_line(aes(y = customers, color = 'customers'))


ridf

as.data.frame(nc_std_ls)

ridf$randreso <- c(rep(NA, length(9)), nc_std_ls1)
ridf$rand <- c(rep(NA, length(9)), nc_std_ls2)

ggplot(data=new_df, aes(x=resolution)) +
  geom_line(aes(y = mean, color = 'sd of number of cluster'))+
  geom_point(aes(y = mean, color = 'mean rand index'))
  
  # geom_point(data=as.data.frame(nc_std_ls), aes(x = nc_std_ls[,1], y = nc_std_ls[,2], color = 'mean rand index'))
  # geom_line(aes(y = ncsd, color = 'sd of number of cluster'))+
  # geom_point(aes(y = ncsd, color = 'sd of number of cluster')) 


new_df$mean[1:9]


Mg_DEG_up.Rds