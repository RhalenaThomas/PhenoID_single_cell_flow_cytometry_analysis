library(Seurat)
library(dplyr)
library(fossil) #for rand index and adjusted rand index
install.packages('flexclust') #another library for rand index, give it a try
library(flexclust)
#======================= Rand index starts =============================

#1. preprocessing
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

#2. function
  randindex <- function(
    df,
    resolutions,
    kn,
    n #number of iterations
    ) {


    rn_ls <- round(runif(n, min=0, max=100000), 0) #random number generator
    
  }
  

  
#===draft

  for(i in kn) {
    #shuming note: why dims=12 here? 
    seu <- FindNeighbors(seu, dims = 1:12, k.param = i) 
    seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
    
    for (j in resolutions) {
      rand_ls <- vector() #list of ri of all 100 repeats
      nc_ls <- vector()
      
      
      
      #store all repeats
      hp1 <- function(x) {
        seu <- FindClusters(seu, random.seed = rn_ls[x], resolution = j)
        seu[[paste("repeat_", x, sep = "")]] <- Idents(object = seu)
        nc_ls <- c(nc_ls, length(unique(seu@meta.data$seurat_clusters)))
      }
      
      sapply(1:n, hp1)
      
      seu <- FindClusters(seu, random.seed = rn_ls[k], resolution = j)
      
      
    }
  }
  
  #diff btw k=1, k=1.1, k=1.3

  seu <- FindClusters(seu, random.seed = 1, resolution = j)
  g1 <- unlist(seu@meta.data$seurat_clusters)
  
  seu <- FindClusters(seu, random.seed = 1.1, resolution = j)
  g2 <- unlist(seu@meta.data$seurat_clusters)
  
  seu <- FindClusters(seu, random.seed = 1.3, resolution = j)
  g3 <- unlist(seu@meta.data$seurat_clusters)
  
  seu <- FindClusters(seu, random.seed = 2, resolution = j)
  g4 <- unlist(seu@meta.data$seurat_clusters)
  
  setequal(g1,g2)
  setequal(g1,g3)
  setequal(g3,g4)
  
  
  # system.time()



  for (i in kn){
    for (j in resolutions) {
      for (k in 1:n) {
        for (z in 1:(k-1)) {
          print(paste(i,j,k,z)) 
        }
      }
    }
  }
  
 
 
  
i <- 60
j <- 0.1
seu <- FindNeighbors(seu, dims = 1:12, k.param = i)
seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
seu <- FindClusters(seu, random.seed = rn_ls[k], resolution = j)


# x <- comPart(unlist(seu[[paste("repeat_", 1, sep = "")]]), unlist(seu[[paste("repeat_", 5, sep = "")]]), type=c("ARI","RI"))


k<-5
seu[[paste("repeat_", k, sep = "")]] <- Idents(object = seu)
tab <- table(unlist(seu[[paste("repeat_", 1, sep = "")]]), unlist(seu[[paste("repeat_", 5, sep = "")]]))
x<-randIndex(tab)


rand_ls <- vector() #list of ri of all 100 repeats
nc_ls <- vector()
if (length(unique(seu@meta.data$seurat_clusters))==1) next #skip the ones 

rand.index(as.numeric(seu@meta.data[row_n, paste("repeat_", z, sep = "")]), as.numeric(seu@meta.data[row_n, paste("repeat_", k, sep = "")]))
seu[[paste("repeat_", k, sep = "")]] <- Idents(object = seu)
nc_ls <- c(nc_ls, length(unique(seu@meta.data$seurat_clusters)))

if (k==1) {next}
for (z in 1:(k-1)) { #z is a number between 1 to k-1, for which k is one of the iterations
  # z <- 4
  # print(z)
  print(paste("i:", i, "j:", j, "k:", k, "z:", z))
  ri <- rand.index(as.numeric(seu@meta.data[row_n, paste("repeat_", z, sep = "")]), as.numeric(seu@meta.data[row_n, paste("repeat_", k, sep = "")]))
  
  # as.numeric(seu@meta.data[row_n, paste("repeat_", z, sep = "")])
  
  # print(paste("k:", k,"z:",z,"ri:",ri))
  rand_ls <- c(rand_ls, ri)
  # print(paste(k,z))
  saveRDS(rand_ls, paste(output_path, "rand_ls_", j, ".Rds",sep=""))
  saveRDS(nc_ls, paste(output_path, "nc_ls_", j, ".Rds",sep=""))
}
#===draft ends
  
#3. input
  input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsretrotransformed_flowset.csv"
  output_path <- "/Users/shumingli/Desktop/output_jun19/"
  

  df <- read.csv(input_path)
  resolutions <- c(0.1)
  kn = c(60) #just need 60  
  n <- 5
  
  #4. test 
  randindex(df, resolutions, kn, n)
#======================= Rand index ends =============================