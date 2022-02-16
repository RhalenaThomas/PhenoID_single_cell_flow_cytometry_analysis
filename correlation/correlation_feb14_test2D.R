############################## correlation ########################################
#we probably don't need these many packages for correlation, just coping it from preprocessed

require("flowCore") #Used for reading the data
require("ggplot2")
require("ggridges") #visualization
require("stringr") #set of functions to manipulate strings type in order to have nice titles
require("rlist") #set of functions that allows to easily manipulate lists
require("reshape2") #visualization
require("flowStats") #Alignment functions
require("scales") #scale colour intensity for visualization
require("dplyr")

# prepossessed matrix
#preprocessed <- read.csv("/Users/shuming/Desktop/PhenoID_single_cell_flow_cytometry_analysis/preprocessing/outputs/prepro_outsflowset.csv")
preprocessed <- read.csv("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/prepro_outsflowset.csv")

# read in the reference matrix
expected_val=read.csv("/Users/rhalenathomas/GITHUB/PhenoID_single_cell_flow_cytometry_analysis/correlation/ReferenceMatrix.csv")

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



# subsample <- preprocessed[sample(1:nrow(m), 5),] #test
subsample <- preprocessed


corr_df <- data.frame(matrix(ncol = 5, nrow = length(subsample)))
colnames(corr_df) <- list("X", "best correlation", "best cell type", 
                          "second correlation", "second cell type")


for (i in 1:nrow(subsample)) {
  #initial values
  best_cor <- -1000000 #maybe replace it with -inf in R? @shuming perfect correlation coefficient is 1
  second_cor <- -1000000
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
  corr_df[i,"best correlation"] <- best_cor
  corr_df[i,"best cell type"] <- best_ct
  corr_df[i,"second correlation"] <- second_cor
  corr_df[i,"second cell type"] <- second_ct
}

head(corr_df)

#2:26:!3
output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/correlations"

#write.csv(corr_df,"/Users/shuming/Desktop/PhenoID_single_cell_flow_cytometry_analysis/correlation/correlation_feb14.csv", row.names = FALSE)
write.csv(corr_df, paste(output_path, "corr_df_r2andcelltype.csv",sep=""))


