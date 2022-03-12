# install.packages("kit")
# install.packages("tidyr")

# clear the environment
rm(list=ls()) 
library(kit) # for finding max and second max (function topn)
library(ggplot2)  # for plotting
library(reshape2) # for plotting (function melt)
library(dplyr) # for df formating (function select)
library(tidyr) #for the last plot in the function

#function: find_correlation
#compare pre-processed expression matrix with the expected value for each cell 
#type, return the best and second best correlated cell types for each sample in 
#the expression matrix
find_correlation <- function(test_path, reference_path, output_path, min_corr=0.1, min_diff=0.05) {
  #input type: test_path and reference_path as strings, min_corr and min_diff as numbers
  test <- read.csv(test_path)
  reference <- read.csv(reference_path)
  
  
  #1. process and scale test and reference matrix =============================
  
  # remove pericytes #what is it again, check after testing
  reference <- reference[-c(9),]
  
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
  
  #replace NA in epithelial,04 with the avg expression of 04 in other cell types
  reference[10,"O4"] <- mean(reference[1:9,"O4"])
  
  #z score markers expression in test (without X)
  test[,markers] <- scale(test[,markers])
  
  # z score the reference matrix 
  reference[,markers] <- scale(reference[,markers])
  
  # test <- test[sample(1:nrow(test), 5),] #testing with 5 samples
  
  #2. find best and second correlation and cell type ==========================
  
  #corr_df will be the output df
  corr_df <- data.frame(matrix(ncol = 6, nrow = length(test)))
  colnames(corr_df) <- list("X", "cor.1", "best.cell.type", 
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
    
    corr_df[i,"X"] <- test[i,1] 
    corr_df[i,"cor.1"] <- corr_ls[top[1]]
    corr_df[i,"best.cell.type"] <- ct_ls[top[1]]
    corr_df[i,"cor.2"] <- corr_ls[top[2]]
    corr_df[i,"second.cell.type"] <- ct_ls[top[2]]
    
    # add variables for cutoffs
    # min_corr - cells are labelled unknown if best corr is less than this threshold
    # min_diff - best - second corr is less than this value then the cells are double labelled
    corr_df[i,"cell.label"] <- ifelse(corr_ls[top[1]] < min_corr, "unknown",
                                      ifelse(corr_ls[top[1]] - corr_ls[top[2]] < min_diff, paste(ct_ls[top[1]],ct_ls[top[2]],sep = "-"), ct_ls[top[1]]))
  }
  
  #3. saving plots and csv===================================================
  
  # save best and second best correlation and cell types as csv 
  write.csv(df, paste(output_path, "corr_celltypes.csv",sep=""), row.names = FALSE)
  
  # filter to get frequency table and save as csv
  df.f <- df %>% select(cell.label)
  freq.table <- as.data.frame(table(df.f))
  write.csv(freq.table, paste(output_path, "Frequencytabletypes.csv",sep=""), row.names = FALSE)
  
  # plot the frequencies and save as pdf
  #plotting after filtering for cell types with more than 100 cells
  # filter
  df.filter <- df %>% group_by(cell.label) %>% dplyr::filter(n()> 100)
  # plot
  pdf(paste(output_path,"FreqCellTypes.pdf",sep=""),width =12, height = 6)
  plot1 <- ggplot(df.filter, aes(x=reorder(cell.label,cell.label,function(x)-length(x)), fill = cell.label))+ geom_bar()+theme_classic() +
    theme(axis.text.x=element_text(angle=90))+ xlab('Assigned cell type') + ylab('number of cell')
  print(plot1)
  dev.off()
  
  df.melt <- melt(df) #reformat to long df
  
  # violin plot of best correlation/cell type
  pdf(paste(output_path,"vlnPlotbestcells.pdf",sep=""))
  plot2 <- ggplot(df, aes(x=best.cell.type, y=cor.1 ))+ geom_violin()+ ylim(-0.1,1)+theme_classic()+
    theme(axis.text.x=element_text(angle=90)) + ylab("correlation coefficient") + xlab("Cell type with max correlation coefficient")
  print(plot2)
  dev.off()
  
  # plot the best and second best correlation together
  pdf(paste(output_path,"boxPlotdoublecelltypes.pdf",sep=""))
  plot3 <- ggplot(df.melt, aes(x=cell.label, y=value ))+ geom_boxplot()+ ylim(-0.1,1)+theme_classic()+
    theme(axis.text.x=element_text(angle=90))+ ylab("correlation coefficient") + xlab("Cell type label")
  print(plot3)
  dev.off()
  
  # plot the best and second best correlation separated on the same graph
  pdf(paste(output_path,"boxPlot2corr.pdf",sep=""))
  plot4 <- ggplot(df.melt, aes(x=best.cell.type, y=value, fill= variable))+ geom_boxplot()+ ylim(-0.25,1)+theme_classic()+
    theme(axis.text.x=element_text(angle=90)) + scale_fill_manual(values = c("#4E84C4", "#52854C")) + ylab("correlation coefficient") + xlab("Cell type")
  print(plot4)
  dev.off()
  
  # the second best correlation is so low it was removed from STEM with axis limit -0.1 and even -1
  
  
  # down sample
  set.seed(64)
  df.downsample <- sample_n(df, 1000)
  df.melt.down <- melt(df.downsample)
  
  # reformat the table to work with the before after plot
  # y is the measurement in df.melt = value
  # x is before after in df.melt = variable
  # class another variable - in the example this is different shapes - for us this is best cell type
  # might use facet to split the cell type - needs to be a factor
  # id is the individual id this is the X column
  pdf(paste(output_path,"pairedPlotBestcelltype.pdf",sep=""))
  plot5 <- ggplot(df.melt.down, aes(x = variable, y = value,colour= variable, group= X)) +
    geom_line(show.legend = F, size = 0.1, color = "black") + geom_point()+ scale_color_manual(values = c("#4E84C4", "#52854C")) + ylim(-0.25,0.95) +
    facet_wrap(~(as.factor(best.cell.type))) +
    #theme_few() +
    #theme(legend.position = "none") +
    ylab("Correlation Coefficient") +
    xlab("")
  print(plot5)
  dev.off()
  
  
  # this will be an excellent visualization but I need to subset only the double labels, then I can plot more cells and see more clearly.
  double.cells <- df[grep("-", df$cell.label),]
  
  df.melt.double <- melt(double.cells)
  
  pdf(paste(output_path,"pairedPlotdoubletypes.pdf",sep=""))
  plot6 <- ggplot(df.melt.double, aes(x = variable, y = value,colour= variable, group= X)) +
    geom_line(show.legend = F, size = 0.1, color = "black") + geom_point()+ scale_color_manual(values = c("#4E84C4", "#52854C")) + ylim(-0.15,0.8) +
    facet_wrap(~(as.factor(cell.label))) +
    ylab("Correlation Coefficient") +
    xlab("")
  print(plot6)
  dev.off()
  
  # compare original cell type input
  # X column has cell type but also the index
  df.split <- df %>% separate(X, ".")
  df.split <- rename(df.split, c("culture.type" = "."))
  pdf(paste(output_path,"PointplotCelllabelvcellculture.pdf",sep=""))
  plot7 <- ggplot(df.split, aes(x =cell.label, y = cor.1, colour = factor(culture.type))) + geom_point(size = 0.05)+
    geom_jitter(width = 0.5, height = 0, size = 0.1) +  theme(axis.text.x=element_text(angle=90))
  print(plot7)
  dev.off()
}

## input 
# test_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv"
# reference_path <- "/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/correlation/ReferenceMatrix10celltypes.csv"
# output_path <- "/Users/shumingli/Desktop/"
test_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsretrotransformed_flowset.csv"
reference_path <- "/Users/rhalenathomas/GITHUB/PhenoID_single_cell_flow_cytometry_analysis/correlation/ReferenceMatrix10celltypes.csv"
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/correlation/retro/8celltypes/"

find_correlation(test_path, reference_path, output_path)