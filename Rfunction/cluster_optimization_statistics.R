library(Seurat)
library(dplyr) #for preprocessing and for plot
# install.packages('flexclust') 
library(flexclust)#for adjusted rand index
library(ggplot2) #for the plot function

#======================= Rand index starts =============================
#1. function 

  randindex <- function(
    df,
    resolutions,
    kn,
    n = 100, #number of iterations
    output_path = NULL #if null, will not save seul, ril, ncl, rdf
    ) {
    
    #preprocessing
    AB <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71",
            "CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")
    df2 <- df %>% select(AB)
    m <- as.matrix(df2)
    tm <- t(df2)
    rownames(tm) <- colnames(df2)
    colnames(tm) <- rownames(df2)
    
    #create seu object
    seu <- CreateSeuratObject(tm)
    seu <- AddMetaData(object=seu, metadata=df$Batch, col.name = 'Batch')
    seu <- ScaleData(seu)
    seu <- RunPCA(seu, features = colnames(df2), npcs = 12, approx = FALSE) 
    #change npcs?
    
    #list of random integers
    rn_ls <- round(runif(n, min = 0, max = 100000), 0)
    
    #final df with mean sd rand index and nc
    rdf <- data.frame(matrix(ncol = 8, nrow = 0))
    colnames(rdf) <- c('kn', 'resolution', 
                       'meanri', 'medianri', 'sdri', 
                       'meannc', 'mediannc', 'sdnc')
    
    #helper functions:
    hp1 <- function(x) { #store n repeats of the same ixj clustering
      seu <- FindClusters(seu, random.seed = rn_ls[x], resolution = j)
      return(Idents(object = seu))
    }
    
    #find rand index between every 2 elements in a list, no repetition 
    #(ex: 12, 21), no same number (ex: 11, 22):
    hp2 <- function(x) { 
      #1. this only caluclates ari, no subsampling needed:
      ri<-randIndex(table(as.numeric(seul[, x[1]]), 
                          as.numeric(seul[, x[2]])))
      
      ##2. this calculates ari and ri, but it's slowers:
      # ri <- comPart(unlist(seu[[paste("repeat_", x[1], sep = "")]]), 
      # unlist(seu[[paste("repeat_", x[2], sep = "")]]), type=c("ARI","RI"))
      
      ##3. this is the fossol function, toooo slow, and has to subsample:
      # ri <- rand.index(as.numeric(seu@meta.data[row_n, paste("repeat_", x[1], sep = "")]), 
      # as.numeric(seu@meta.data[row_n, paste("repeat_", x[2], sep = "")]))
      return(ri)
    }
    
    #find number of clusters
    hp3 <- function(x) { return((length(unique(x)))) }
    
    #main loops:
    for(i in kn) {
      #shuming note: why dims=12 here? 
      seu <- FindNeighbors(seu, dims = 1:12, k.param = i) 
      seu <- RunUMAP(seu, dims = 1:12, n.neighbors = i)
      for (j in resolutions) {
        seul <- sapply(1:n, hp1) #list of n repeats of clustering (Ident(seu) object)
        ncl <- apply(seul, 2, hp3) #list of number of clustering 
        ril <- apply(t(combn(1:n, 2)), 1, hp2) #list of rand index
        
        if (!is.null(output_path)) {
          saveRDS(seul,paste(output_path, "seu_ls_kn", i, "_j", j, ".Rds",sep=""))
          saveRDS(ncl,paste(output_path, "nc_ls_kn", i, "_j", j, ".Rds",sep=""))
          saveRDS(ncl,paste(output_path, "ri_ls_kn", i, "_j", j, ".Rds",sep=""))
        }
        
        rdf <-rbind(rdf, list(
                    kn = i, resolution = j,
                    meanri = mean(ril), medianri = median(ril), sdri = sd(ril), 
                    meannc = mean(ncl), mediannc = median(ncl), sdnc = sd(ncl)))
      }
    }
    return(rdf)
  }
  

  
  
#3. input
  input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsretrotransformed_flowset.csv"
  
  df <- read.csv(input_path)
  resolutions <- c(0.1)
  kn = c(60) #just need 60  
  n <- 3
  output_path <- "/Users/shumingli/Desktop/output_jul6/"
  
  #4. test 
  randindex(df, resolutions, kn, n)
#======================= Rand index ends =============================
  
#======================= Rand index plot starts =============================

#1. function
plot.randindex <- function (
    rdf,
    cp = c("orange", "violet"),
    view = c(0, 1) #zoom in x axis, this does not changes scales, just the viewable sections
) {
  
  s <- (view[2]-view[1])/max(rdf$meannc+rdf$sdnc)
  
  p <- rdf %>% 
    ggplot(aes(x = resolution)) +
    geom_line(aes(y = meanri), color = cp[1]) +
    geom_point(aes(y = meanri), color = cp[1], size=1)+ 
    geom_errorbar(aes(ymin=meanri-sdri, ymax=meanri+sdri), color = cp[1], width=.01,
                  position=position_dodge(.9))+
    geom_line(aes(y = meannc*s+view[1]), color = cp[2]) +
    geom_point(data = rdf, mapping = aes(x = resolution, y = meannc*s+view[1]), color = cp[2])+ 
    geom_errorbar(rdf, mapping = aes(x=resolution, ymin=((meannc-sdnc)*s+view[1]), ymax=((meannc+sdnc)*s+view[1])), width=.01,
                  position=position_dodge(.9), color = cp[2])+
    scale_y_continuous(limits= view, name="Mean Rand Index",
                       sec.axis = sec_axis(~ . /s - view[1]/s, 
                                           name = "Mean Number of Clusters"))+
    theme(axis.text.y  = element_text(color = cp[1]),
          axis.title.y = element_text(color=cp[1]),
          axis.text.y.right =  element_text(color = cp[2]),
          axis.title.y.right = element_text(color=cp[2]),
          plot.title = element_text(hjust = 0.5, size=10))+
    ggtitle("Plot of Mean Rand Index and \n Mean Number of Clusters")
  
  return(p)
}

# 2. input
rdf <- read.csv('/Users/shumingli/Desktop/output_jun19/randIndex0.05.csv')

#3. test 
fig <- plot.randindex(rdf, c('pink','violet'), c(0.7,1))

  

# paper figures:
# rdf %>% 
#   ggplot(aes(x = resolution)) +
#   geom_line(aes(y = meanri), color = "blue") +
#   geom_point(aes(y = meanri), color = "blue")+ 
#   geom_errorbar(aes(ymin=meanri-sdri, ymax=meanri+sdri), color = 'blue', width=.01,
#                 position=position_dodge(.9))+
#   geom_line(aes(y = meannc/22.62000*0.3+0.7), color = "red") +
#   geom_point(data = rdf, mapping = aes(x = resolution, y = meannc/22.62000*0.3+0.7), color = "red")+ 
#   geom_errorbar(rdf, mapping = aes(x=resolution, ymin=((meannc-sdnc)/22.62000*0.3+0.7), ymax=((meannc+sdnc)/22.62000*0.3+0.7)), width=.01,
#                 position=position_dodge(.9), color = 'red')+
#   scale_y_continuous(limits= c(0.7, 1.05), name="mean rand index",
#                      sec.axis = sec_axis(~ . /0.3*22.62000 - 0.7/0.3*22.62000, 
#                                          name = "mean number of clusters"))+
#   theme(axis.text.y  = element_text(color = 'blue'),
#                        axis.title.y = element_text(color='blue'),
#                        axis.text.y.right =  element_text(color = 'red'),
#                        axis.title.y.right = element_text(color='red'),
#                        plot.title = element_text(hjust = 0.5, size=10))+
#   ggtitle("Plot of mean rand index and \n mean number of clusters")


#======================= Rand index plot ends =============================
  