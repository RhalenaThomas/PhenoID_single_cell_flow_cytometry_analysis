# make scatter plots of expression - pair by pair


############################ load libraries ###################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(gridExtra)

# input path
input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/test/Louvain/"
filename <- "FlowsetSeuratSeuratObject100.Rds"

# read in a seurat object 
seu <- readRDS(paste(input_path,filename,sep=""))


##################################### make scatter plots in a loop and save #############################################
Idents(seu) <- 'orig.ident'  # instead of group by set the grouping ident

print(FeatureScatter(seu, feature1 = "AQP4", feature2 = "CD44", slot = 'scale.data'))

# make a loop to compare each pair - need to make some kind of grid plot to save.

Idents(seu) <- 'Batch'  # instead of group by set the grouping ident

print(FeatureScatter(seu, feature1 = "AQP4", feature2 = "CD44", slot = 'scale.data'))

# vector of the antibodies (features to compare)
antibodies <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")

p <- list()

for (AB in antibodies){
  for (i in 1:13){
  p[[i]]<- FeatureScatter(seu, feature1 = AB, feature2 = antibodies[i], slot = 'scale.data')
  } 
  do.call(grid.arrange, p)
}


