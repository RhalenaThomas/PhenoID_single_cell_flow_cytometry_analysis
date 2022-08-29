# statistics on full annotation
rm(list=ls())

# loading the appropriate libraries

### for anova test: 
library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table) #for transpose()
library(reshape2)#used to rename melted df

library(readr)
library(multcompView)
# are these 2 used?

library(scProportionTest) ### for prop test 
# devtools::install_github("rpolicastro/scProportionTest")



########start of the preprcessing#########

## read in the annotated Seurat object
# outpath <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/"
# input_path <- paste("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/","All9MOannaote.12072022.Rds")

input_path <- "/Users/shumingli/Downloads/AllCellsFinal.Rds"
output_path <- "/Users/shumingli/Desktop/output_jul7/"

input <- readRDS(input_path)
df <- transpose(as.data.frame(GetAssayData(input, slot = 'scale.data')))
AB <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")

# I think the data should be in the way it appears in the input df, which 
# AB <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")
colnames(df) <- AB

# add the cluster labels and the sample IDs
Sample.input <- as.matrix(input@meta.data$Batch)
cell.types <- input@meta.data$cell.types
df.2 <- cbind(df, label = cell.types, Sample = Sample.input)

#add new columns batch, genotype, experiment day and age based on Batch (the old column)
# this could also be added from the meta data in the seurat object

#batch
df.2[which(grepl("0317A", df.2$Sample)), "batch"] <- "A" #batch A (0317A = A)
df.2[which(grepl("0317B|0306", df.2$Sample)), "batch"] <- "B" #batch B (0317B ir 0306 = B)

#genotype
df.2[which(grepl("3450", df.2$Sample)), "genotype"] <- "3450" # genotype 3450
df.2[which(grepl("AIW002", df.2$Sample)), "genotype"] <- "AIW002" # genotype AIW002
df.2[which(grepl("AJG001C", df.2$Sample)), "genotype"] <- "AJG001C" # genotype AJG001C

#experiment date
df.2[which(grepl("0306", df.2$Sample)), "exdate"] <- "0306" # experiment date 0306
df.2[which(grepl("0317", df.2$Sample)), "exdate"] <- "0317" # experiment date 0317

#age
df.2[which(grepl("0306", df.2$Sample)), "age"] <- "273" # age 273, same as ex date 0306
df.2[which(grepl("0317B", df.2$Sample)), "age"] <- "284" # age 284
df.2[which(grepl("0317A", df.2$Sample)), "age"] <- "263" # age 263

#make a df with 9 entries (9 samples)
#take the mean of each sample (for categorical, there should only be 1 value) 
df.3 <- c()
for (i in unique(df.2$Sample)) {
  for (j in AB) {df.3 <- c(df.3, mean(df.2[which(df.2$Sample == i), j]))}
  df.3 <- c(df.3, i)
  for (j in iv) {df.3 <- c(df.3, unique(df.2[which(df.2$Sample == i), j]))}
}

df.3 <- matrix(lapply(as.list(df.3), type.convert, as.is=TRUE), ncol = length(c(AB,'Sample', iv)), byrow = TRUE)
colnames(df.3) <- c(AB,'Sample', iv)

melt_df.3 <- melt(df.3, measure.vars = AB, variable.name = 'antibody', value.name = 'AB.mean.expression') # variable = AB name, value = AB mean expression/sample

#### create a df with the mean for each cell type for each of the 9 samples
melt.df2 <- melt(df.2, measure.vars = AB, variable.name = 'antibody')

### I do want the other variables to be added back
df.group <- melt.df2 %>% group_by(label, Sample, antibody, batch, genotype, exdate, age
) %>% dplyr::summarize(Mean = mean(value, na.rm=TRUE))

########end of the preprcessing#########

######### run ANOVAs ####################

#function
twanova <- function(A,  #one of independent variables, i.e. one of c('batch', 'genotype', 'exdate', 'age')
                    B, #antibody or label 
                    dv, #one of dependent variable, antibody mean value here
                    df, #input dataframe
                    group, #group to loop through and to use as sample, antibody or label. 
                    output_path = NULL) { #if output_path is given, save stats lists
                    # If B=antibody, group=label; if B=label, group=antibody
  
  #lists to store stats 
  anv.l <- c() #store anova 
  tk.l <- vector( "list" , 4) #store tukey 
  tk.l.useful <- vector( "list" , 4) #store useful comparison of tukey in interaction 
  tk.l.sig <- c() #store significant outcome (p<0.05) within useful tukey comparison in interaction
  
  # i <- "Glia Lineage"
  
  for (i in unique(unlist(df[, group]))){
    # run 2 way anova with interaction effect
    # subset the cell type
    df.group.sub <- df[which(df[, group]== i), ]
    df.group.sub <- as.data.frame(df.group.sub)
    res.aov2 <- aov(df.group.sub[, dv] ~ df.group.sub[, B] * df.group.sub[, A]) #diff
    

    anv.l <- c(anv.l, 
               i, B, summary(res.aov2)[[1]][["Pr(>F)"]][1], 
               i, A, summary(res.aov2)[[1]][["Pr(>F)"]][2], 
               i, paste(B,':', A, sep=''), summary(res.aov2)[[1]][["Pr(>F)"]][3]
               )

    # now the posthoc tests
    tukey.results <- TukeyHSD(res.aov2)
    

    for (j in 1:length(tukey.results)) {
      itrt = switch(j, B, A, paste(B,':', A, sep=''))
      tk.l[[1]] <- c(tk.l[[1]], rep(i, nrow(tukey.results[[j]])))
      tk.l[[2]] <- c(tk.l[[2]], rep(itrt, nrow(tukey.results[[j]])))
      tk.l[[3]] <- c(tk.l[[3]], rownames(tukey.results[[j]]))
      tk.l[[4]] <- c(tk.l[[4]], unname(tukey.results[[j]][, 'p adj']))
    }
    
    length(tukey.results)
    nrow(tukey.results[[]])
    
    #variable to the left of ":" is always B
    rowl <- c() #row numbers of comparisons that are useful, only same B is useful
    for (j in 1:nrow(tukey.results[[3]])) {
      split_rowname <- strsplit(strsplit(rownames(tukey.results[[3]])[j], split = "-")[[1]], split = ":")
      if (split_rowname[[1]][1] == split_rowname[[2]][1]) {
        # print(rownames(tukey.results[[3]])[j])
        rowl <- c(rowl, j)
      }
    }
    
    for (j in 1:2) {
      itrt = switch(j, B, A)
      tk.l.useful[[1]] <- c(tk.l.useful[[1]], rep(i, nrow(tukey.results[[j]])))
      tk.l.useful[[2]] <- c(tk.l.useful[[2]], rep(itrt, nrow(tukey.results[[j]])))
      tk.l.useful[[3]] <- c(tk.l.useful[[3]], rownames(tukey.results[[j]]))
      tk.l.useful[[4]] <- c(tk.l.useful[[4]], unname(tukey.results[[j]][, 'p adj']))
    }
    tk.l.useful[[1]] <- c(tk.l.useful[[1]], rep(i, length(rowl)))
    tk.l.useful[[2]] <- c(tk.l.useful[[2]], rep(paste(B,':', A, sep=''), length(rowl)))
    tk.l.useful[[3]] <- c(tk.l.useful[[3]], rownames(tukey.results[[3]][rowl, ]))
    tk.l.useful[[4]] <- c(tk.l.useful[[4]], unname(tukey.results[[3]][rowl, 'p adj']))
  }
  
  # # filter for interactions that are significant
  # tk.l.sig[[i]] <- tukey.df[rowl, ] %>% filter(`p adj`<= 0.05)
  
  anv <- as.data.frame(matrix(anv.l, ncol = 3, byrow=TRUE))
  colnames(anv) <- c('group', 'subgroup', 'p.val')
  
  tk <- data.frame(group = tk.l[[1]],
                          subgroup = tk.l[[2]],
                          comparison = tk.l[[3]],
                          p.val = tk.l[[4]])
  
  tk.useful <- data.frame(group = tk.l.useful[[1]],
                     subgroup = tk.l.useful[[2]],
             comparison = tk.l.useful[[3]],
             p.val = tk.l.useful[[4]])
  
  tk.sig <- tk.useful %>% filter(`p.val`<= 0.05)

  
  # save results
  if(!is.null(output_path)) {
    saveRDS(anv, paste(output_path,"anova.", A, ".", B, ".RDS", sep=''))
    saveRDS(tk, paste(output_path,"tukey.", A, ".", B, ".RDS", sep=''))
    saveRDS(tk.useful, paste(output_path,"tukey.sig.useful.", A, ".", B, ".RDS", sep='')) 
    saveRDS(tk.sig, paste(output_path,"tukey.sig.", A, ".", B, ".RDS", sep=''))
    
    write.csv(anv, paste(output_path,"anova.", A, ".", B, ".csv", sep=''), row.names = FALSE)
    write.csv(tk, paste(output_path,"tukey.", A, ".", B, ".csv", sep=''), row.names = FALSE)
    write.csv(tk.sig, paste(output_path,"tukey.sig.useful.", A, ".", B, ".csv", sep=''), row.names = FALSE) 
    write.csv(tk.useful, paste(output_path,"tukey.sig.", A, ".", B, ".csv", sep=''), row.names = FALSE)
  }
  return(list(anv = anv, tk = tk, tk.useful = tk.useful, tk.sig = tk.sig))
}

#parameters
A <- 'genotype'
B <- 'antibody' 
dv <- 'Mean'
df <- df.group 
group <- 'label' #label or antibody, if B is antibody, then group is label, if B is label, then group is antibody

#test
test <- twanova(A, B, dv, df, group, output_path)




### prop test ####

#function:
prop_test <- function(input_path, ci, s1, s2, si, output_path = NULL) {
  seu <- readRDS(input_path)
  prop_test <- sc_utils(seu)
  # set up the comparison
  # need to compare separately
  prop_test <- permutation_test(
    prop_test, cluster_identity = ci,
    sample_1 = s1, sample_2 = s2,
    sample_identity = si
  )
  p <- permutation_plot(prop_test)
  # save the plot 
  if(!is.null(output_path)) {
    png(paste(output_path,"prp.test.AIWvAJG.png"))
    p
    dev.off()
  }
  return(p)
}

#parameters:
input_path <- "/Users/shumingli/Downloads/AllCellsFinal.Rds"
output_path <- "/Users/shumingli/Desktop/output_jul7/"
ci <- "cell.types"
s1 <- "AIW002"
s2 <- "AJG001" # one of c("AJG001", "3450", "AIW002")
si <- "Genotype"

#test:
test <- prop_test(input_path, ci, s1, s2, si, output_path)




# # make a group of AIW vs other
# Idents(seu) <- "Batch"
# cluster.ids <- c("AJG-3450","AJG-3450","AJG-3450","AIW002","AIW002","AIW002",
#                  "AJG-3450","AJG-3450","AJG-3450")
# names(cluster.ids) <- levels(seu)
# seu <- RenameIdents(seu, cluster.ids)
# seu$ipsc <- Idents(seu)
# DimPlot(seu, group.by = 'ipsc')
# 
# # join together AJG and 3450
# prop_test <- permutation_test(
#   prop_test, cluster_identity = "cell.types",
#   sample_1 = "AIW002", sample_2 = "AJG-3450",
#   sample_identity = "ipsc")
# # make the plot
# permutation_plot(prop_test)
# # save the plot 
# pdf(paste(output_path,"prp.test.AIWvs3450-AJG.pdf"),width = 9, height = 4)
# permutation_plot(prop_test) + theme_bw() + 
#   theme(axis.text.x=element_text(size=15),
#         axis.text.y=element_text(size=15))
# dev.off()