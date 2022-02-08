# make scatter plots of expression - pair by pair


############################ load libraries ###################################

library(Seurat)
library(dplyr)
library(ggplot2)

# input path
input_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/test/Louvain/"
filename <- "FlowsetSeuratSeuratObject100.Rds"

# read in a seurat object 
seu <- readRDS(paste(input_path,filename,sep=""))


##################################### make scatter plots in a loop and save #############################################

print(FeatureScatter(seu, feature1 = "AQP4", feature2 = "CD44"))


