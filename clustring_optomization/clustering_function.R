# flowsom clustering
# parameter tuning, statistic and visualizations for manual annotation

# if(!require(devtools)){
#   install.packages("devtools") # If not already installed
# }
# devtools::install_github("JinmiaoChenLab/Rphenograph")




rm(list=ls()) 
library(clusterSim) #new package for dbi
library(FlowSOM)
library(flowCore)
library(cluster)
library(fpc)
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)
library(Rphenograph)
library(reshape2) #for plotting multiple lines (resolutions) on the same graph

library(kit) # for finding max and second max (function topn)
library(tidyr) #for the last plot in the function



#function: find_correlation
#compare pre-processed expression matrix with the expected value for each cell 
#type, return the best and second best correlated cell types for each sample in 
#the expression matrix
find_correlation <- function(test_path, reference_path, output_path, min_corr=0.1, min_diff=0.05) {
  #input type: test_path and reference_path as strings, min_corr and min_diff as numbers
  
  test <- read.csv(test_path)
  reference <- read.csv(reference_path)
  
  #replace NA in epithelial,04 with the avg expression of 04 in other cell types
  reference[9,"O4"] <- mean(reference[1:8,"O4"])
  
  #1. process and scale test and reference matrix =============================
  
  #select_col takes the markers that exist both in reference and expression matrix
  #in reference's order
  select_col <- list()
  
  for (i in colnames(reference)) {
    for (j in colnames(test)) {
      if (tolower(i) == tolower(j)) {select_col <- c(select_col, j)}
    }
  } 
  
  
  #change reference's spelling to match the test's spelling
  colnames(reference) <- select_col
  
  #select 13 markers + X in test
  test <- test %>% select(colnames(reference))
  
  #a list of markers (without X)
  markers <- unlist(select_col[-1])
  
  #z score markers expression in test (without X)
  test[,markers] <- scale(test[,markers])
  
  # z score the reference matrix 
  reference[,markers] <- scale(reference[,markers])
  
  # test <- test[sample(1:nrow(test), 5),] #testing with 5 samples
  
  #2. find best and second correlation and cell type ==========================
  
  #df will be the output df
  df <- data.frame(matrix(ncol = 6, nrow = length(test)))
  colnames(df) <- list("X", "cor.1", "best.cell.type", 
                       "cor.2", "second.cell.type","cell.label")
  
  
  #the loop that will find the best and second best correlation and cell types
  for (i in 1:nrow(test)) {
    corr_ls <- vector() #list of correlation between the reference cell types and each sample
    ct_ls <- vector() #cell type list
    for (j in 1:nrow(reference)) {
      corr <- cor(as.numeric(test[i,markers]),as.numeric(reference[j,markers])) # pearson by default and we use default
      corr_ls <- c(corr_ls, corr)
      ct_ls <- c(ct_ls, reference[j,1])
    }
    top <- topn(corr_ls, 2L, decreasing = TRUE) #return the index of the best 2
    
    df[i,"X"] <- test[i,1]
    df[i,"cor.1"] <- corr_ls[top[1]]
    df[i,"best.cell.type"] <- ct_ls[top[1]]
    df[i,"cor.2"] <- corr_ls[top[2]]
    df[i,"second.cell.type"] <- ct_ls[top[2]]
    
    # add variables for cutoffs
    # min_corr - cells are labelled unknown if best corr is less than this threshold
    # min_diff - best - second corr is less than this value then the cells are double labelled
    df[i,"cell.label"] <- ifelse(corr_ls[top[1]] < min_corr, "unknown",
                                 ifelse(corr_ls[top[1]] - corr_ls[top[2]] < min_diff, paste(ct_ls[top[1]],ct_ls[top[2]],sep = "-"), ct_ls[top[1]]))
  } 
  
  #3. saving plots and csv===================================================
  
  # save best and second best correlation and cell types as csv
  write.csv(df, paste(output_path, "corr_celltypes.csv",sep=""), row.names = FALSE)
  
  # filter to get frequency table and save as csv
  df.f <- df %>% select(cell.label)
  
  freq.table <- as.data.frame(table(df.f))
  write.csv(freq.table, paste(output_path, "Frequencytabletypes.csv",sep=""), row.names = FALSE)
  # 
  # # plot the frequencies and save as pdf
  # #plotting after filtering for cell types with more than 100 cells
  # # filter
  df.filter <- df %>% group_by(cell.label) %>% dplyr::filter(n()> 100)
  # # plot
  pdf(paste(output_path,"FreqCellTypes.pdf",sep=""),width =12, height = 6)
  plot1 <- ggplot(df.filter, aes(x=reorder(cell.label,cell.label,function(x)-length(x)), fill = cell.label))+ geom_bar()+theme_classic() +
    theme(axis.text.x=element_text(angle=90))+ xlab('Assigned cell type') + ylab('number of cell')
  print(plot1)
  dev.off()
  # 
  df.melt <- melt(df) #reformat to long df
  # 
  # # violin plot of best correlation/cell type
  pdf(paste(output_path,"vlnPlotbestcells.pdf",sep=""))
  plot2 <- ggplot(df, aes(x=best.cell.type, y=cor.1 ))+ geom_violin()+ ylim(-0.1,1)+theme_classic()+
    theme(axis.text.x=element_text(angle=90)) + ylab("correlation coefficient") + xlab("Cell type with max correlation coefficient")
  print(plot2)
  dev.off()
  # 
  # # plot the best and second best correlation together
  pdf(paste(output_path,"boxPlotdoublecelltypes.pdf",sep=""))
  plot3 <- ggplot(df.melt, aes(x=cell.label, y=value ))+ geom_boxplot()+ ylim(-0.1,1)+theme_classic()+
    theme(axis.text.x=element_text(angle=90))+ ylab("correlation coefficient") + xlab("Cell type label")
  print(plot3)
  dev.off()
  # 
  # # plot the best and second best correlation separated on the same graph
  pdf(paste(output_path,"boxPlot2corr.pdf",sep=""))
  plot4 <- ggplot(df.melt, aes(x=best.cell.type, y=value, fill= variable))+ geom_boxplot()+ ylim(-0.25,1)+theme_classic()+
    theme(axis.text.x=element_text(angle=90)) + scale_fill_manual(values = c("#4E84C4", "#52854C")) + ylab("correlation coefficient") + xlab("Cell type")
  print(plot4)
  dev.off()
  # 
  # # the second best correlation is so low it was removed from STEM with axis limit -0.1 and even -1
  # 
  # 
  # # down sample
  set.seed(64)
  df.downsample <- sample_n(df, 1000)
  df.melt.down <- melt(df.downsample)
  # 
  # # reformat the table to work with the before after plot
  # # y is the measurement in df.melt = value
  # # x is before after in df.melt = variable
  # # class another variable - in the example this is different shapes - for us this is best cell type
  # # might use facet to split the cell type - needs to be a factor
  # # id is the individual id this is the X column
  pdf(paste(output_path,"pairedPlotBestcelltype.pdf",sep=""))
  plot5 <- ggplot(df.melt.down, aes(x = variable, y = value,colour= variable, group= X)) +
    geom_line(show.legend = F, size = 0.1, color = "black") + geom_point()+ scale_color_manual(values = c("#4E84C4", "#52854C")) + ylim(-0.25,0.95) +
    facet_wrap(~(as.factor(best.cell.type))) +
    theme(legend.position = "none") +
    ylab("Correlation Coefficient") +
    xlab("")
  print(plot5)
  dev.off()
  # 
  # 
  # # this will be an excellent visualization but I need to subset only the double labels, then I can plot more cells and see more clearly.
  double.cells <- df[grep("-", df$cell.label),]
  # 
  df.melt.double <- melt(double.cells)
  # 
  pdf(paste(output_path,"pairedPlotdoubletypes.pdf",sep=""))
  plot6 <- ggplot(df.melt.double, aes(x = variable, y = value,colour= variable, group= X)) +
    geom_line(show.legend = F, size = 0.1, color = "black") + geom_point()+ scale_color_manual(values = c("#4E84C4", "#52854C")) + ylim(-0.15,0.8) +
    facet_wrap(~(as.factor(cell.label))) +
    ylab("Correlation Coefficient") +
    xlab("")
  print(plot6)
  dev.off()
  
}


#helper function: stats plotting function
flowsom_stats_plot <- function(stats_ls, output_path, input_name, clust_method) {
  #silhouette score:
  #-1: bad clusters  0: neutral, indifferent  1: good clusters
  pdf(paste(output_path,input_name,clust_method,'statssilhouette.pdf',sep=""),width =4, height = 4)
  siplot1 <- ggplot(stats_ls) + 
    geom_point(aes(x=krange, y=si)) +
    labs(title = "Silhouette Scores",
         x = "krange", y = "Average Silhouette Scores") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  siplot2 <- ggplot(stats_ls) + 
    geom_point(aes(x=nc, y=si)) +
    labs(title = "Silhouette Scores",
         x = "Number of Clusters", y = "Average Silhouette Scores") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(siplot1)
  print(siplot2)
  dev.off()
  
  #Calinski-Harabasz index: 
  # the highest value is the optimal number of clusters
  pdf(paste(output_path,input_name,clust_method,'statsCalinskiHara.pdf',sep=""),width =4, height = 4)
  chplot1 <- ggplot(stats_ls) + 
    geom_point(aes(x=krange, y=ch)) +
    labs(title = "Calinski-Harabasz Index",
         x = "krange", y = "Calinski-Harabasz Index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  
  chplot2 <- ggplot(stats_ls) + 
    geom_point(aes(x=nc, y=ch)) +
    labs(title = "Calinski-Harabasz Index",
         x = "Number of Clusters", y = "Calinski-Harabasz Index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(chplot1)
  print(chplot2)
  dev.off()
  
  #Davies–Bouldin index: minimum score is zero
  #the lowest value is the optimal number of clusters
  pdf(paste(output_path,input_name,clust_method,'statsDavies.pdf',sep=""),width =4, height = 4)
  dbplot1 <- ggplot(stats_ls) + 
    geom_point(aes(x=krange, y=db)) +
    labs(title = "Davies–Bouldin index",
         x = "krange", y = "Davies–Bouldin index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  dbplot2 <- ggplot(stats_ls) + 
    geom_point(aes(x=nc, y=db)) +
    labs(title = "Davies–Bouldin index",
         x = "Number of Clusters", y = "Davies–Bouldin index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(dbplot1)
  print(dbplot2)
  
  dev.off()
}


#helper function: stats plotting function
phenograph_stats_plot <- function(stats_ls, output_path, input_name, clust_method) {
  #silhouette score:
  #-1: bad clusters  0: neutral, indifferent  1: good clusters
  pdf(paste(output_path,input_name,clust_method,'statssilhouette.pdf',sep=""),width =4, height = 4)
  siplot1 <- ggplot(stats_ls) + 
    geom_point(aes(x=kn, y=si)) +
    labs(title = "Silhouette Scores",
         x = "kn", y = "Average Silhouette Scores") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  siplot2 <- ggplot(stats_ls) + 
    geom_point(aes(x=nc, y=si)) +
    labs(title = "Silhouette Scores",
         x = "Number of Clusters", y = "Average Silhouette Scores") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(siplot1)
  print(siplot2)
  dev.off()
  
  #Calinski-Harabasz index: 
  # the highest value is the optimal number of clusters
  pdf(paste(output_path,input_name,clust_method,'statsCalinskiHara.pdf',sep=""),width =4, height = 4)
  chplot1 <- ggplot(stats_ls) + 
    geom_point(aes(x=kn, y=ch)) +
    labs(title = "Calinski-Harabasz Index",
         x = "kn", y = "Calinski-Harabasz Index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  
  chplot2 <- ggplot(stats_ls) + 
    geom_point(aes(x=nc, y=ch)) +
    labs(title = "Calinski-Harabasz Index",
         x = "Number of Clusters", y = "Calinski-Harabasz Index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(chplot1)
  print(chplot2)
  dev.off()
  
  #Davies–Bouldin index: minimum score is zero
  #the lowest value is the optimal number of clusters
  pdf(paste(output_path,input_name,clust_method,'statsDavies.pdf',sep=""),width =4, height = 4)
  dbplot1 <- ggplot(stats_ls) + 
    geom_point(aes(x=kn, y=db)) +
    labs(title = "Davies–Bouldin index",
         x = "kn", y = "Davies–Bouldin index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  dbplot2 <- ggplot(stats_ls) + 
    geom_point(aes(x=nc, y=db)) +
    labs(title = "Davies–Bouldin index",
         x = "Number of Clusters", y = "Davies–Bouldin index") +
    theme(plot.title = element_text(hjust = 0.5)) 
  
  print(dbplot1)
  print(dbplot2)
  
  dev.off()
}



#helper function: louvain stats plot
louvain_stats_plotting <- function(stats_ls, output_path, input_name, clust_method) {
  
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

#helper function: louvain clustering
#= c(25,50,100,125,150,200,250,300)
louvain_clustering <- function(input_path, 
                               output_path, 
                               input_name, 
                               clust_method, 
                               kn,
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



#helper function: flowsom clustering
#= c(5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90)
flowsom_clustering <- function(krange,
                               input_path, output_path, input_name, clust_method) {
  
  
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
  col.names <- colnames(df2)
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
  seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)
  
  # I've tried scaling the kn with the k but the values to no result in UMAP that spatial match cluster
  # I'll just run the UMAP once with the kn = square root of the number of inputs
  
  kn = round(sqrt(dim(df2)[1]))
  seu <- FindNeighbors(seu, dims = 1:12, k = kn)
  seu <- RunUMAP(seu, dims = 1:12, n.neighbors = kn)
  # save feature plots of this UMAP
  # just for testing print
  
  # we also only need to plot the features once
  # file name
  UMAP_name = paste("UMAPfeatures_kn",kn,".pdf",sep="")
  print(UMAP_name) #testing
  
  # save feature plots UMAP
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
  print(FeaturePlot(seu, features = AB,slot = 'scale.data',min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ theme(plot.title = element_text(size = 0.1)))
  dev.off()
  
  # here k is the number of clusters
  #shuming: somehow 2 doesn't work with flowsom, im not suring why
  #krange = 3:30 
  #krange = seq(from = 5, to = 100, by = 5)
  # this cause a problem - it doesn't include the last 2 values in the loop
  
  # the k will be the max k for the metaclustering clustering. 
  # save a data object for each kn - will only keep temporarily
  # the clusters will write over with each new kn
  
  
  #create a df to store all stats, col = nc, si, ch, db
  stats_ls <- data.frame(matrix(ncol = 5, nrow = length(krange)))
  colnames(stats_ls) <- c("krange", "nc","si", "ch", "db")
  count <- 0
  
  
  #subsample for silhouette score
  #shuming: here im using 1000 so it's not too slow, but 30,000 would be have better representation 
  row_n <- sample(1:nrow(m), 1000)
  dis <- dist(m[row_n,])
  
  
  for (i in krange){
    # K max number of clusters not the kn input
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
    number.clusters <- length(unique(flowSOMcluster[fs$map$mapping[,1]]))
    
    ### make umap
    #UMAP_name = paste("UMAPclusters_k",i,".pdf",sep="")
    #print(UMAP_name) #testing
    #pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =8, height = 5)
    # save UMAP grouped
    #print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
    #dev.off()
    
    # the pdf don't work well for a quick figure
    UMAP_name = paste("UMAPclusters_k",i,".png",sep="")
    print(UMAP_name) #testing
    png(paste(output_path,input_name,clust_method,UMAP_name,sep=""))
    # save UMAP grouped
    print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
    dev.off()
    
    # heatmap
    #heatmap_name = paste("Heatmapclusters_k",i,".pdf",sep="")
    #testing
    #pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""),width =8, height = 5)
    #print(DoHeatmap(seu, features = AB,group.by = clust_name))
    #dev.off()
    # heatmap
    heatmap_name = paste("Heatmapclusters_k",i,".png",sep="")
    #testing
    png(paste(output_path,input_name,clust_method,heatmap_name,sep=""), width = 600, height = 500)
    print(DoHeatmap(seu, features = AB,group.by = clust_name, size = 10) +theme(text = element_text(size = 30)))
    dev.off()
    
    #### add stats
    # calculate the statistics
    
    #number of clusters 
    # "krange", "nc","si", "ch", "db"
    count <- 1+count
    
    stats_ls[count, "krange"] <- i 
    
    #number of clusters:
    stats_ls[count, "nc"] <- number.clusters # calculated above
    
    #silhouette score:
    stats_ls[count, "si"] <- mean(silhouette(flowSOMcluster[row_n],dis)[, 3])
    
    #Calinski-Harabasz index: 
    stats_ls[count, "ch"] <- calinhara(m,flowSOMcluster,cn=i)
    
    # Davies–Bouldin index:
    stats_ls[count, "db"] <- index.DB(df2, as.numeric(flowSOMcluster))$DB
    
  }

  flowsom_stats_plot(stats_ls, output_path, input_name, clust_method)
  
  
  
  
  # make clustree plot
  
  
  pdf(paste(output_path,input_name,clust_method,'Clustree.pdf',sep=""),width =8, height = 8)
  print(clustree(seu, prefix ='FlowSom.k.'))
  dev.off()
  
  # save the UMAP with cell types
  
  pdf(paste(output_path,input_name,clust_method,'UMAPcelltype.pdf',sep=""),width =8, height = 6)
  print(DimPlot(seu,group.by = 'Batch'))
  dev.off()
  
  # save the Seurat object
  saveRDS(seu,paste(output_path,input_name,clust_method,'SeuratObject.Rds',sep=""))
  
  # save the stats list
  saveRDS(stats_ls,paste(output_path,input_name,clust_method,'statslist.Rds',sep=""))
  write.csv(stats_ls, paste(output_path, "stats.csv",sep=""), row.names = FALSE)
  
}



#helper function: phenograph clustering
#= c(25,50,75,100,125,150,175,200,225,250,300,350,400,450,500)
phenograph_clustering <- function(kn, 
                                  input_path, output_path, input_name, clust_method) {
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
  print("test1")
  seu <- ScaleData(seu)
  print("test2")
  # check the data
  pdf(paste(output_path,input_name,clust_method,"Heatmap_batch.pdf",sep=""),width =8, height = 6)
  print(DoHeatmap(seu, group.by = "Batch", features = AB))
  dev.off()
  print("test3")
  # create the UMAP
  seu <- RunPCA(seu, features = AB, npcs = 12, approx = FALSE)
  
  print("test4")
  # like FlowSOM Phenograph doesn't relate directly to the UMAP like Louvain
  # we will make on seurat UMAP and visualize the clusters there
  
  
  kn_umap = round(sqrt(dim(df2)[1]))
  seu <- FindNeighbors(seu, dims = 1:12, k.param = kn_umap)
  print("test5")
  seu <- RunUMAP(seu, dims = 1:12, n.neighbors = kn_umap)
  print("test6")
  # save feature plots of this UMAP
  # just for testing print

  # we also only need to plot the features once
  # file name
  UMAP_name = paste("UMAPfeatures_kn",kn_umap,".pdf",sep="")
  print(UMAP_name) #testing
  
  # save feature plots UMAP
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
  print(FeaturePlot(seu, features = AB,slot = 'scale.data',min.cutoff = 'q1', max.cutoff ='99',label.size = 1)+ theme(plot.title = element_text(size = 0.1)))
  dev.off()
  print("test7")
  # we also want to see the batch on the UMAP
  pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =8, height = 6)
  print(DimPlot(seu, group.by = 'Batch'))
  dev.off()
  print("test8")
  
  ############################## explore parameters and calculate statistics ###########################
  
  
  #create 3 lists for stats
  
  #create a df to store all stats, col = nc, si, ch, db
  stats_ls <- data.frame(matrix(ncol = 5, nrow = length(kn)))
  colnames(stats_ls) <- c("kn", "nc","si", "ch", "db")
  count <- 0
  
  #subsampling for silhouette score, n=1000, can make n bigger if needed
  set.seed(25)
  row_n <- sample(1:nrow(m), 1000)
  dis <- dist(m[row_n,])

  print("test9")
  
  ############################# loop to explore parameters ########################################
  # kn = c(25,50,75,100,125,150,175,200,225,250,300,350,400,450,500)
  # kn = c(25,50,75,100,125,150,175,200,225,250,275,300)
  # larger kn fewer clusters in general but not always
  #kn = c(50,500)
  # save a data object for each kn - will only keep temporarily
  # the clusters will write over with each new kn


  for (i in kn){
    
    ### run phenograph clustering
    Rphenograph_out_flow <- Rphenograph(m, k = i)
    print("test10")
    clust_name = paste('Pheno.kn.',i,sep="")
    # add the cluster ID into seurat object to visualize
    seu <- AddMetaData(object=seu, factor(membership(Rphenograph_out_flow[[2]])), col.name = clust_name) 
    print("test11")
    number.clusters <- length(unique(factor(membership(Rphenograph_out_flow[[2]]))))
    print("test12")
    ### make umap 
    
    UMAP_name = paste("UMAPclusters_kn",i,".pdf",sep="")
    print(UMAP_name) #testing 
    pdf(paste(output_path,input_name,clust_method,UMAP_name,sep=""),width =20, height = 10)
    # save UMAP grouped
    print(DimPlot(seu,reduction = "umap", repel = TRUE, label = TRUE, group.by = clust_name)) # will automatically group by active ident
    dev.off()
    # heatmap
    heatmap_name = paste("Heatmapclusters_kn",i,".pdf",sep="")
    #testing 
    pdf(paste(output_path,input_name,clust_method,heatmap_name,sep=""),width =25, height = 10)
    print(DoHeatmap(seu, features = AB,group.by = clust_name))
    dev.off()
    print("test13")
    #### add stats
    
  
  
    # "kn", "nc","si", "ch", "db"
    count <- 1+count
    
    # get the cluster indexes 
    phenocluster <- factor(membership(Rphenograph_out_flow[[2]]))
    print("test14")
    
    stats_ls[count, "kn"] <- i 
    
    #number of clusters:
    stats_ls[count, "nc"] <- number.clusters # calculated above
    
    #silhouette score:
    stats_ls[count, "si"] <- mean(silhouette(as.numeric(phenocluster[row_n]),dis)[, 3])
    
    #Calinski-Harabasz index: 
    stats_ls[count, "ch"] <- calinhara(m,phenocluster,cn=i)
    
    # Davies–Bouldin index:
    stats_ls[count, "db"] <- index.DB(df2, as.numeric(phenocluster))$DB
    print("test15")
  }


  # make clustree plot
  pdf(paste(output_path,input_name,clust_method,'Clustree.pdf',sep=""),width =15, height = 10)
  print(clustree(seu, prefix ='Pheno.kn.'))
  dev.off()
  print("test16")
  # save the Seurat object
  saveRDS(seu,paste(output_path,input_name,clust_method,'SeuratObject.Rds',sep=""))
  
  
  # save the stats list

  saveRDS(stats_ls,paste(output_path,input_name,clust_method,'statslist.Rds',sep=""))
  
  phenograph_stats_plot(stats_ls, output_path, input_name, clust_method)
  
  write.csv(stats_ls, paste(output_path, "stats.csv",sep=""), row.names = FALSE)
  print("test17")
}


#MAIN FUNCTION function with all clustering methods
clustering <- function(resolution=NULL, kn, input_path, output_path, input_name, clust_method) {
  case_when(clust_method=="Louvain" ~ louvain_clustering(input_pat=input_path, output_path=output_path, input_name=input_name, clust_method=clust_method, kn=kn),
            clust_method=="FlowSOM" ~ flowsom_clustering(input_path=input_path, output_path=output_path, input_name=input_name, clust_method=clust_method, krange=kn, resolution=resolution),
            clust_method=="Pheno" ~ phenograph_clustering(input_path=input_path, output_path=output_path, input_name=input_name, clust_method=clust_method,kn=kn),
  )
}


#plot stats for all
#takes 3 stats list (louvain, flowsom, phenograph) and compare each stat 
#for all three clustering types
#input is df, so should read in rds or csv as df beforehand
plot_comparison <- function(input_name, output_path, louvain, flowsom, phenograph) {
  
  #krange in flowsom=kn? 
  colnames(flowsom)[colnames(flowsom) == "krange"] <- "kn"
  flowsom['resolution'] <- NA
  phenograph['resolution'] <- NA
 
  #add a column for clustering methodd to each df
  flowsom <- cbind(method="flowsom",flowsom)
  phenograph <- cbind(method="phenograph",phenograph)
  louvain <- cbind(method="louvain",louvain)
  
  #rename the clustering method for louvain with kn  
  for (i in row(louvain)) {
    louvain[i, 'method'] <- paste("louvain ",louvain[i,'kn'],sep="")
  } 
  

  #merge 3 df in one file (louvain method + kn)
  stats_ls <-rbind(flowsom, phenograph, louvain)
  
  #rename the clustering method for louvain with resolution
  for (i in row(louvain)) {
    louvain[i, 'method'] <- paste("louvain ",louvain[i,'resolution'],sep="")
  } 
  
  #create a df that contains louvain method + resolution
  stats_ls2 <-rbind(flowsom, phenograph, louvain)
  
  #reorder the method so that the legend is displayed in order
  stats_ls$method <- factor(stats_ls$method,
                            levels = unique(stats_ls$method))
  
  stats_ls2$method <- factor(stats_ls2$method,
                            levels = unique(stats_ls2$method))

  #color palette, can add more color to it later
  palette <- c("#d62828", "#fd9e02", "#a9d6e5", "#89c2d9", "#61a5c2", "#468faf", "#2c7da0", "#2a6f97", "#014f86","#01497c", "#013a63", "#012a4a")
  
  #plot comparison for silhouette score
  pdf(paste(output_path,input_name,clust_method,"Silhouette_plot_comparison.pdf",sep=""))
  
  siplot1 <- ggplot(stats_ls, aes(x=nc, y=si, label=method)) +
    geom_point(aes(group=method,color=method), size=1) +
    geom_line(aes(group=method,color=method), size=0.1) +
    scale_colour_manual(values=palette) +
    labs(color = "method", title = "Silhouette Scores",
         x = "Number of Clusters", y = "Average Silhouette Scores") +
    theme(plot.title = element_text(hjust = 0.5))
  print(siplot1)
  
  siplot2 <- ggplot(stats_ls2, aes(x=kn, y=si, label=method)) +
    geom_point(aes(group=method,color=method), size=1) +
    geom_line(aes(group=method,color=method), size=0.1) +
    scale_colour_manual(values=palette) +
    labs(color = "method", title = "Silhouette Scores",
         x = "kn", y = "Average Silhouette Scores") +
    theme(plot.title = element_text(hjust = 0.5))
  print(siplot2)

  dev.off()
  
  
  
  #plot comparison for Calinski-Harabasz index
  pdf(paste(output_path,input_name,clust_method,"calinski_harabasz_plot_comparison.pdf",sep=""))
  
  chplot1 <- ggplot(stats_ls, aes(x=nc, y=ch, label=method)) +
    geom_point(aes(group=method,color=method), size=1) +
    geom_line(aes(group=method,color=method), size=0.1) +
    scale_colour_manual(values=palette) +
    labs(color = "method", title = "Calinski-Harabasz index",
         x = "Number of Clusters", y = "Calinski-Harabasz index") +
    theme(plot.title = element_text(hjust = 0.5))
  print(chplot1)
  
  chplot2 <- ggplot(stats_ls2, aes(x=kn, y=ch, label=method)) +
    geom_point(aes(group=method,color=method), size=1) +
    geom_line(aes(group=method,color=method), size=0.1) +
    scale_colour_manual(values=palette) +
    labs(color = "method", title = "Calinski-Harabasz index",
         x = "kn", y = "Calinski-Harabasz index") +
    theme(plot.title = element_text(hjust = 0.5))
  print(chplot2)
  
  dev.off()
  
  
  #plot comparison for Davies–Bouldin index
  pdf(paste(output_path,input_name,clust_method,"davies_bouldin_plot_comparison.pdf",sep=""))
  
  dbplot1 <- ggplot(stats_ls, aes(x=nc, y=db, label=method)) +
    geom_point(aes(group=method,color=method), size=1) +
    geom_line(aes(group=method,color=method), size=0.1) +
    scale_colour_manual(values=palette) +
    labs(color = "method", title = "Davies–Bouldin index",
         x = "Number of Clusters", y = "Davies–Bouldin index") +
    theme(plot.title = element_text(hjust = 0.5))
  print(dbplot1)
  
  dbplot2 <- ggplot(stats_ls2, aes(x=kn, y=db, label=method)) +
    geom_point(aes(group=method,color=method), size=1) +
    geom_line(aes(group=method,color=method), size=0.1) +
    scale_colour_manual(values=palette) +
    labs(color = "method", title = "Davies–Bouldin index",
         x = "kn", y = "Davies–Bouldin index") +
    theme(plot.title = element_text(hjust = 0.5))
  print(dbplot2)
  
  dev.off()
}


#testing for plot comparison

#shuming's input path
louvain <- read.csv("/Users/shumingli/Desktop/louvain_output/stats.csv")
flowsom <- readRDS("/Users/shumingli/Desktop/flowsom_output/FlowsetFlowSOMstatslist.Rds")
phenograph <- readRDS("/Users/shumingli/Desktop/phenograph_output/FlowAlignTransPhenostatslist.Rds")

input_name <- "AlignTrans"
output_path <- "/Users/shumingli/Desktop/"

plot_comparison(input_name=input_name, output_path=output_path, louvain=louvain, flowsom=flowsom, phenograph=phenograph)






# #Testing for phenogrpah
# 
# # input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsaligned_transformed_flowset.csv"
# input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
# 
# # output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/cluster_parameters/Pheno/"
# output_path <- "/Users/shumingli/Documents/"
# 
# # add input description to output files
# input_name <- "FlowAlignTrans"  # this will be the different processing types
# 
# # cluster type for file name
# clust_method <- "Pheno"
# 
# kn = c(25, 50)
# 
# clustering(kn=kn, input_path=input_path, output_path=output_path, input_name=input_name, clust_method=clust_method)


# 
# ############# testing for louvain ############################
# 
# 
# input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsaligned_transformed_flowset.csv"
# output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/cluster_parameters/Louvain/"
# # input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
# # output_path <- "/Users/shumingli/Desktop/"
# 
# 
# input_name <- "Flowset"  # processing type for file name
# clust_method <- "Louvain" # cluster type for file name
# 
# louvain_clustering(input_path, output_path, input_name, clust_method)
# 
# 
# 
# 
# 
# 
# ############# testing for flowsom  ############################
# 
# define the input pathway
# input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsaligned_transformed_flowset.csv"
# input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
# 
# # output pathway
# # output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/cluster_parameters/FlowSom-cat/"
# output_path <- "/Users/shumingli/Desktop/"
# 
# 
# 
# # add input description to output files
# input_name <- "Flowset"  # this will be the different processing types
# 
# # cluster type for file name
# clust_method <- "FlowSOM"
# 
# flowsom_clustering(input_path=input_path, output_path=output_path, input_name=input_name, clust_method=clust_method)
# 
# 
# 
# 
# 
# ############# set up the data object for phenograph clustering ############################

# # info to change for each comparison
# # define the input pathway
# # input pathway
# # input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsaligned_transformed_flowset.csv"
# input_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
# 
# # output pathway
# # output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/cluster_parameters/Pheno/"
# output_path <- "/Users/shumingli/Desktop/"
# 
# # add input description to output files
# input_name <- "FlowAlignTrans"  # this will be the different processing types
# 
# # cluster type for file name
# clust_method <- "Pheno"
# 
# kn = c(25,50,75,100,125,150,175,200,225,250,300,350,400,450,500)
# 
# phenograph_clustering(kn=kn, input_path=input_path, output_path=output_path, input_name=input_name, clust_method=clust_method)


# stats_ls <- readRDS("/Users/shumingli/Desktop/FlowAlignTransPhenostatslist.Rds")
# 
# phenograph_stats_plot(stats_ls, output_path, input_name, clust_method)





