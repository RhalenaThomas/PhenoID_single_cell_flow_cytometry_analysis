

############################## correlation ########################################

# clear the environment
rm(list=ls()) 

library(ggplot2)  # for plotting

library(reshape2) # for melt

library(dplyr) # for df formating

# prepossessed matrix

# 9000 cells per 9 MBO
# use aligned sample transformed
#preprocessed <- read.csv("/Users/shuming/Desktop/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsaligned_transformed_flowset.csv")
preprocessed <- read.csv("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsaligned_transformed_flowset.csv")
# 2D cultures
#preprocessed <- read.csv("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/2Dcells_surface/preprocessing/select/2DcellsSelectflowset.csv")


# read in the reference matrix
expected_val=read.csv("/Users/rhalenathomas/GITHUB/PhenoID_single_cell_flow_cytometry_analysis/correlation/ReferenceMatrix10celltypes.csv")

# expected list: X AQP4	CD24	CD44	CD184	CD15	HepaCam	CD29	CD56	O4	CD140a	CD133	Glast	CD71
# changed to: X AQP4 CD24  CD44  CD184  CD15  HepaCAM  CD29  CD56 O4 CD140a  CD133  GLAST  CD71

#there's probably a more elegant way to do it
#make a list of markers with preprocessed's spelling but in expected's order
select_col <- list()
for (ev in colnames(expected_val)) {
  for (pp in colnames(preprocessed)) {
    if (tolower(ev) == tolower(pp)) {select_col <- c(select_col, pp)}
  }
} 

#change expected_val's spelling to match the proprocessed's spelling
colnames(expected_val) <- select_col

#select 13 markers + X in preprocessed
preprocessed <- preprocessed %>% select(colnames(expected_val))

#a list of markers (without X)
markers <- unlist(select_col[-1])

#z score markers expression (without X)
preprocessed[,markers] <- scale(preprocessed[,markers])

# z score the reference matrix expected value
expected_val[,markers] <- scale(expected_val[,markers])

# subsample <- preprocessed[sample(1:nrow(m), 5),] #test
subsample <- preprocessed


corr_df <- data.frame(matrix(ncol = 6, nrow = length(subsample)))
colnames(corr_df) <- list("X", "best.correlation", "best.cell.type", 
                          "second.correlation", "second.cell.type","cell.label")


for (i in 1:nrow(subsample)) {
  #initial values
  best_cor <- -1 #maybe replace it with -inf in R? @shuming perfect correlation coefficient is 1
  second_cor <- -1
  best_ct <- ""
  
  for (j in 1:nrow(expected_val)) {
    corr <- cor(as.numeric(subsample[i,markers]),as.numeric(expected_val[j,markers]))
    if ((is.na(corr) == FALSE) & (corr > best_cor)) {
      #update the second best correlation & cell type
      second_cor <- best_cor 
      second_ct <- best_ct
      
      #update the best correlation & cell type
      best_cor <- corr
      best_ct <- expected_val[j,1]
    }
  }
  corr_df[i,"X"] <- subsample[i,1]
  corr_df[i,"best.correlation"] <- best_cor
  corr_df[i,"best.cell.type"] <- best_ct
  corr_df[i,"second.correlation"] <- second_cor
  corr_df[i,"second.cell.type"] <- second_ct
  corr_df[i,"cell.label"] <- ifelse(best_cor < 0.1, "unknown", 
                                    ifelse(best_cor - second_cor < 0.05, paste(best_ct,second_ct,sep = "-"),best_ct))
}

# now add a column for cell labels
df <- corr_df

#df$cell.lable <- ifelse(df$`best.correlation`<0.1, "unknown", ifelse(df$`best.correlation` - df$`second.correlation` < 0.05, paste(df$`best.cell.type`,df$`second.cell.type`,sep = "-"),df$`best.cell.type`))



# filter to get frequency table
df.f <- df %>% select(cell.label)

freq.table <- as.data.frame(table(df.f))




#2:26:!3
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/correlation/"

#write.csv(corr_df,"/Users/shuming/Desktop/PhenoID_single_cell_flow_cytometry_analysis/correlation/correlation_feb14.csv", row.names = FALSE)
write.csv(df, paste(output_path, "corr_celltypes.csv",sep=""))

write.csv(freq.table, paste(output_path, "Frequencytabletypes.csv",sep=""))


# plot the frequencies
#plotting after filtering for cell types with more than 10 or 100 cells
# filter 
df.filter <- df %>% group_by(cell.label) %>% dplyr::filter(n()> 100)

pdf(paste(output_path,"FreqCellTypes.pdf",sep=""),width =8, height = 6)
ggplot(df.filter, aes(x=cell.label))+ geom_bar()+theme_classic()+
  theme(axis.text.x=element_text(angle=90))
dev.off()




# need to sample before melt to keep matching data.
df <- rename(df, c("best.correlation"= "cor.1", "second.correlation" = "cor.2"))

#### reformat to long df

df.melt <- melt(df)

# second corr values are sometime so low it messes up the plot

# vln plot of the 

pdf(paste(output_path,"vlnPlotbestcells.pdf",sep=""))
ggplot(df, aes(x=best.cell.type, y=cor.1 ))+ geom_violin()+ ylim(-0.1,1)+theme_classic()+
  theme(axis.text.x=element_text(angle=90)) + ylab("correlation coefficient") + xlab("Cell type with max correlation coefficient")
dev.off()

# this plot the best and second best corr are together 
pdf(paste(output_path,"boxPlotdoublecelltypes.pdf",sep=""))
ggplot(df.melt, aes(x=cell.label, y=value ))+ geom_boxplot()+ ylim(-0.1,1)+theme_classic()+
  theme(axis.text.x=element_text(angle=90))+ ylab("correlation coefficient") + xlab("Cell type label")
dev.off()

# splitting the max and 2nmax corr
pdf(paste(output_path,"boxPlot2corr.pdf",sep=""))
ggplot(df.melt, aes(x=best.cell.type, y=value, fill= variable))+ geom_boxplot()+ ylim(-0.25,1)+theme_classic()+
  theme(axis.text.x=element_text(angle=90)) + scale_fill_manual(values = c("#4E84C4", "#52854C")) + ylab("correlation coefficient") + xlab("Cell type")

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
ggplot(df.melt.down, aes(x = variable, y = value,colour= variable, group= X)) +
  geom_line(show.legend = F, size = 0.1, color = "black") + geom_point()+ scale_color_manual(values = c("#4E84C4", "#52854C")) + ylim(-0.25,0.95) +
  facet_wrap(~(as.factor(best.cell.type))) +
  #theme_few() +
  #theme(legend.position = "none") +
  ylab("Correlation Coefficient") +
  xlab("")
dev.off()


# this will be an excellent visualization but I need to subset only the double labels, then I can plot more cells and see more clearly. 

double.cells <- df[grep("-", df$cell.label),]

df.melt.double <- melt(double.cells)

pdf(paste(output_path,"pairedPlotdoubletypes.pdf",sep=""))
ggplot(df.melt.double, aes(x = variable, y = value,colour= variable, group= X)) +
  geom_line(show.legend = F, size = 0.1, color = "black") + geom_point()+ scale_color_manual(values = c("#4E84C4", "#52854C")) + ylim(-0.15,0.8) +
  facet_wrap(~(as.factor(cell.label))) +
  ylab("Correlation Coefficient") +
  xlab("")
dev.off()

# compare original cell type input
# X column has cell type but also the index
library(tidyr)
df.split <- df %>% separate(X, ".")
df.split <- rename(df.split, c("." = "culture.type"))
pdf(paste(output_path,"PointplotCelllabelvcellculture.pdf",sep=""))
ggplot(df.split, aes(x =cell.label, y = cor.1, colour = factor(culture.type))) + geom_point(size = 0.05)+
  geom_jitter(width = 0.5, height = 0, size = 0.1) +  theme(axis.text.x=element_text(angle=90))
dev.off()
