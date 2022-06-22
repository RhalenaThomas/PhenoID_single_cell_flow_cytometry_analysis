# explore clustering results

# load libraries
library(clusterSim) 
library(FlowSOM)
library(flowCore)
library(cluster)
library(fpc)
library(clv)
library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)
library(reshape2)

# read in the seurat objects

seu <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/cluster_parameters/retro-louv-moreparm/retrotransLouvainSeuratObject20.Rds")

# Feature list
# "AQP4"    "CD24"    "CD44"    "CD184"   "CD15"    "HepaCAM" "CD29"    "CD56"    "O4"      "CD140a" 
# "CD133"   "GLAST"   "CD71"  

# make dotplot
DotPlot(seu, features = AB, group.by = 'seurat_clusters') + +RotatedAxis()+coord_flip()+theme(axis.text.x=element_text(size=7),axis.text.y=element_text(size=7))