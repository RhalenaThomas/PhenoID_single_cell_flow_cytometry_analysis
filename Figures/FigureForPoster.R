# analyse the 9 MBO full dataset - make figures for poster

# set up the environment 

library(Seurat)
library(dplyr)
library(ggplot2)
library(clustree)
library(reshape2)


# read in the annotated data ojbect
seu <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outs/clusters/AllcellLablesMarch25.Rds")

outpath <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PosterFigures/"

DimPlot(seu, group.by = 'cell.types')
DimPlot(seu, group.by = 'labels.groups')
# go back to the many labels of cell type groups

DimPlot(seu, group.by = 'labels')
# I'm going to rename some of the clusters. I plotted the ordered Dot plot and look to see what is really the same groups.
# current names in lables
# c("Mix","NPC","Astro1","RG1","RG2","N-E-NPC-Astro","Astro2","Astro3","RG3","N-OPC-Mix","Astro4","Endothelial","Neurons1","EarlyNeurons","Neurons2","Astro5","Neurons3","NPC-Neuro","Oligo","Epithelial","Neurons4")

Idents(seu) <- 'labels'

cluster.ids <- c("Mix","Neural Precursors","Astrocytes 1","Radial Glia 1","Radial Glia 2","Glia","Astrocytes 1","Astrocytes 1","Radial Glia 1","OPC","Astrocytes 2","Endothelial","Neurons 1","Neurons 2","Neurons 1","Astrocytes 2","Neurons 1","Neural Precursors","Oligodendrocytes","Epithelial","Mix")
names(cluster.ids) <- levels(seu)    # get the levels
seu <- RenameIdents(seu, cluster.ids) # rename  
seu$cluster.ids <- Idents(seu)   # add a new dataslot

# make a vector of the desired plotting order
cell.type.order <- c("Neurons 1","Neurons 2","Neural Precursors","Astrocytes 1","Astrocytes 2","Glia","Radial Glia 1","Radial Glia 2","Oligodendrocytes","Epithelial","Endothelial","OPC","Mix")
cell.type.order <- rev(cell.type.order)

colours <- c("mediumorchid2","purple3","mediumpurple1","darkorange1","lightsalmon2","orange1","cyan",
             "cyan3","darkblue","aquamarine","aquamarine2","darkgrey","azure3")

DimPlot(seu, group.by = 'cluster.ids', label = TRUE, repel = TRUE,order = cell.type.order, cols = colours, shuffle = TRUE, label.size = 4)

Idents(seu) <- 'cluster.ids'
png(paste(outpath,"UMAPlabelledClusters.png"),width = 1100, height = 600)
 DimPlot(seu, order = cell.type.order, cols = colours, shuffle = TRUE, raster=FALSE,pt.size = 0.02) +
 theme(legend.text = element_text(size=30), axis.title.y = element_text(size=25), 
       axis.title.x = element_text(size=25), axis.text.y = element_text(size =25),axis.text.x = element_text(size =25))
dev.off()

# get the mean expression per cluster to do the dot plots and heatmaps

express.by.cluster <- as.data.frame(AverageExpression(seu, features = AB, group.by = 'cluster.ids', slot = 'scale.data'))
express.by.cluster <- as.data.frame(scale(express.by.cluster))

col.names.RNA <- colnames(express.by.cluster)

# have to find a better renaming methods - if I just get levels the order is not the same
# I'm taking cluster id names 
re.col.name <- c("Mix","Neural Precursors","Astrocytes 1","Radial Glia 1","Radial Glia 2","Glia","OPC","Astrocytes 2","Endothelial","Neurons 1","Neurons 2","Oligodendrocytes","Epithelial")


names(express.by.cluster) <- re.col.name
AB2 <- row.names(express.by.cluster)
ex.data <- cbind(AB2,express.by.cluster)


# DOTPLOT
# can I get the proportion of cells expressing a marker?
a <- DotPlot(seu, features = AB, group.by = 'cluster.ids')
a$data
pct.exp <- as.data.frame(a$data) %>% select(features.plot, id, pct.exp)

# add the mean expression and the percent cells expressing together
# rename the AB and cell labels 
colnames(longData) <- c("features.plot","id","expression")

df.exp.pct <- merge(longData, pct.exp, by = c("features.plot", "id"))


# reorder the data so that it is easier to read

new.order <- c("Neurons 1","Neurons 2","Neural Precursors","Astrocytes 1","Astrocytes 2","Glia","Radial Glia 1","Radial Glia 2","Oligodendrocytes","Epithelial","Endothelial","OPC","Mix")
new.order <- rev(new.order)
AB.order <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")


data <- df.exp.pct %>% mutate(Cell.type = factor(id, levels = new.order))
data <- data %>% mutate(Marker = factor(features.plot, levels = AB.order))

ggplot(data = data, aes(x = Marker, y = Cell.type, color = expression, size = pct.exp)) +
  geom_point() +
  scale_color_gradient(low = "grey", high = "blue") + ylab("Cell Phenotypes") + xlab("Antibodies") + RotatedAxis() +theme_classic()

# make the heatmap
# want to reorder again

new.order <- c("Neurons 1","Neurons 2","Neural Precursors","Astrocytes 1","Astrocytes 2","Glia","Radial Glia 1","Radial Glia 2","Oligodendrocytes","Epithelial","Endothelial","OPC","Mix")

AB.order <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")
AB.order <- rev(AB.order)

data <- df.exp.pct %>% mutate(Cell.type = factor(id, levels = new.order))
data <- data %>% mutate(Marker = factor(features.plot, levels = AB.order))


png(paste(outpath,"HeatmapCelltypes.png"),height = 550, width = 500)
ggplot(data, aes(x = Cell.type, y = Marker)) + 
  geom_raster(aes(fill=expression)) + 
  scale_fill_gradient(name = "Relative\nExpression",low="darkorchid4", high="cyan1", na.value = "grey", 
                      breaks = c(min(data$expression),mean(data$expression),max(data$expression)),labels = c("Min","Mean","Max")) +
  labs(x="Cell types", y="Antibodies") +
  theme_bw() + theme(axis.text.x=element_text(size=20, colour = "black", angle=90, hjust=1),
                     axis.text.y=element_text(size=18, colour = "black"),
                     axis.title.x = element_text(size=24, colour = "black"), 
                     axis.title.y = element_text(size=24, colour = "black"),legend.text = element_text(size=15),
                     legend.title=element_text(size=18))
dev.off()



#### make a split plot between the different MBO
# I'll split by genotype
# Genotype
Idents(seu) <- 'Batch'
cluster.ids <- c("3450","3450","3450","AIW002","AIW002","AIW002","AJG001C","AJG001C","AJG001C")
# vector with the new names - you need this vector from me

names(cluster.ids) <- levels(seu)    # get the levels
seu <- RenameIdents(seu, cluster.ids) # rename  
seu$Genotype <- Idents(seu)   # add a new dataslot

seu$Genotype <- factor(x = seu$Genotype, levels = c("AIW002","AJG001C","3450")) # to reorder the split.by
Idents(seu) <- 'cluster.ids'
# save the plot Split DimPlot
cell.type.order <- c("Neurons 1","Neurons 2","Neural Precursors","Astrocytes 1","Astrocytes 2","Glia","Radial Glia 1","Radial Glia 2","Oligodendrocytes","Epithelial","Endothelial","OPC","Mix")
cell.type.order <- rev(cell.type.order)

# add everything to match the UMAP above
png(paste(outpath,"UMAPsplitbyGenotype.png"), width= 1500, height = 350)
DimPlot(seu, reduction = "umap", label = FALSE, split.by = 'Genotype', order = cell.type.order, cols = colours, 
        shuffle = TRUE, raster=FALSE,pt.size = 0.02) + theme(legend.text = element_text(size=20), axis.title.y = element_text(size=22), 
        axis.title.x = element_text(size=22), axis.text.y = element_text(size =22),axis.text.x = element_text(size =22))
dev.off()


# make bar charts to show proportion of cells 

sample.lables <- as.data.frame(table(seu$Genotype, seu$cluster.ids))

# bar chart of with percent 
ggplot(sample.lables, aes(x = Var1,y=Freq ,fill = Var2)) + 
  geom_bar(position= "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + theme_classic() + theme(text = element_text(size=15),
  axis.text.x = element_text(angle=90, hjust=1))  + xlab('Genotype') + ylab('Percent of Cell type') 

# now I need to reorder, recolour and resize the fonts


new.order <- c("Neurons 1","Neurons 2","Neural Precursors","Astrocytes 1","Astrocytes 2","Glia","Radial Glia 1","Radial Glia 2","Oligodendrocytes","Epithelial","Endothelial","OPC","Mix")

# data <- df.OG %>% mutate(NewColumnName = factor(current_column_name, levels = new.order))
sample.lables <- sample.lables %>% mutate(Cell.Types = factor(Var2, levels = new.order))

# save bar chart with proportions

png(paste(outpath,"BarChartGenotypePrpCells.png"), width = 550, height = 400)
ggplot(sample.lables, aes(x = Var1, y=Freq, fill = Cell.Types)) + 
  geom_bar(position= "fill", stat = "identity") + 
  scale_y_continuous(labels = scales::percent_format()) + theme_classic() + theme(legend.text = element_text(size=18),
                      text = element_text(size=21, colour = "black"),
                      axis.text.x = element_text(angle=90, hjust=1, colour = "black", size = 18), 
                      axis.text.y = element_text(colour = "black", size= 16))+ 
                      xlab('Genotype') + 
                      ylab('Percent of Cell type') +
                      scale_fill_manual(values= colours, name = " ")
dev.off()


#### now I need some why to compare expression between groups 
##### get the expression values 

exp.by.cell.geneotype <- as.data.frame(AverageExpression(seu, features = AB, group.by = c('cluster.ids','Genotype'), slot = 'scale.data'))
exp.by.cell.geneotype <- as.data.frame(scale(exp.by.cell.geneotype))
# add the rownames as a column

AB2 <- row.names(exp.by.cell.geneotype)
ex.data.c <- cbind(AB2,exp.by.cell.geneotype)


# reformat
longData<- melt(ex.data.c)

#now we need to split the variable Cell type and Genotype to have less crazy labels - maybe use facet for each cell type.
# maybe filter out some cell types 


# heatmap

ggplot(longData, aes(x = variable, y = AB2)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="blue", high="pink", na.value = "grey") +
  labs(x="Cell types", y="Antibodies") +
  theme_bw() + theme(axis.text.x=element_text(size=12, angle=90, vjust=0.3),
                     axis.text.y=element_text(size=12),
                     plot.title=element_text(size=12))



### for now I'll save as a pdf and relabel in illustrator 

# I can remove some of the cell types 
colnames(ex.data.c)
ex.data.c.2 <- ex.data.c %>% select(-c(RNA.Mix_AIW002,RNA.Mix_AJG001C,RNA.Mix_3450,RNA.OPC_3450,RNA.OPC_AIW002,RNA.OPC_AJG001C,
                                       RNA.Oligodendrocytes_AJG001C,RNA.Oligodendrocytes_AIW002,RNA.Oligodendrocytes_AJG001C,
                                       RNA.Endothelial_3450,RNA.Endothelial_AIW002,RNA.Endothelial_AJG001C,RNA.Epithelial_3450,
                                       RNA.Epithelial_AIW002,RNA.Epithelial_AJG001C,RNA.Oligodendrocytes_3450))
#Rename so I can reorder
colnames(ex.data.c.2)
re.col.name <- c("Marker","NeuralPrecursors AIW002",  "NeuralPrecursors AJG001C",
                 "NeuralPrecursors 3450","Astrocytes1 AIW002","Astrocytes1 AJG001C",    
                 "Astrocytes1 3450", "RadialGlia1 AIW002", "RadialGlia1 AJG001C",    
                 "RadialGlia1 3450", "RadialGlia2 AIW002", "RadialGlia2 AJG001C",    
                 "RadialGlia2 3450", "Glia AIW002", "Glia AJG001C",              
                 "Glia_3450", "Astrocytes2 AIW002","Astrocytes2 AJG001C",
                  "Astrocytes2 3450", "Neurons1 AIW002", "Neurons1 AJG001C",
                  "Neurons1 3450", "Neurons2 AIW002",         
                  "Neurons2 AJG001C", "Neurons2 3450")           
names(ex.data.c.2) <- re.col.name
longdata.2 <- melt(ex.data.c.2)


new.order <- c("Neurons1 AIW002", "Neurons1 AJG001C",
               "Neurons1 3450","Neurons2 AIW002",         
               "Neurons2 AJG001C", "Neurons2 3450","NeuralPrecursors AIW002",  "NeuralPrecursors AJG001C",
               "NeuralPrecursors 3450","Astrocytes1 AIW002","Astrocytes1 AJG001C",    
               "Astrocytes1 3450","Astrocytes2 AIW002","Astrocytes2 AJG001C",
               "Astrocytes2 3450","Glia AIW002", "Glia AJG001C",              
               "Glia_3450","RadialGlia1 AIW002", "RadialGlia1 AJG001C",    
               "RadialGlia1 3450", "RadialGlia2 AIW002", "RadialGlia2 AJG001C",    
               "RadialGlia2 3450")
#new.order <- rev(new.order)
AB.order <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")
AB.order <- rev(AB.order)

data.2 <- longdata.2 %>% mutate(Cell.Type = factor(variable, levels = new.order))
data.2 <- data.2 %>% mutate(Marker = factor(Marker, levels = AB.order))

pdf(paste(outpath,"Heatmap.cell.genotype.pdf"),width = 6, height = 3.5)
ggplot(data.2, aes(x = Cell.Type, y = Marker)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="blue", high="pink", na.value = "grey") +
  labs(x="Cell types", y="Antibodies") +
  theme_bw() + theme(axis.text.x=element_text(size=8, angle=90, hjust=1),
                     axis.text.y=element_text(size=8),
                     plot.title=element_text(size=8))
dev.off()



###### to make a screen shot of heatmap #####
DoHeatmap(seu, group.by = 'seurat_clusters', features = AB)


