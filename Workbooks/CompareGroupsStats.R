# statistics on full annotation


library(Seurat)
library(dplyr)
library(ggplot2)
library(data.table) #for transpose()
library(reshape2)#used to rename melted df

# read in the annotated Seurat object


outpath <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/"


input_path <- paste("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/","All9MOannaote.12072022.Rds")
input <- readRDS(input_path)

df <- transpose(as.data.frame(GetAssayData(input,slot = 'scale.data')))

 
AB <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")

# I think the data should be in the way it appears in the input df, which 
# AB <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")
colnames(df) <- AB

# add the cluster labels and the sample IDs
Sample.input <- as.data.frame(input@meta.data$Batch)
cell.types <- input@meta.data$cell.types


df.2 <- cbind(df, label = cell.types)
df.2 <- cbind(df.2, Sample=Sample.input)
names(df.2)[names(df.2) == "input@meta.data$Batch"] <- "Sample"


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


######### run ANOVAs ####################
iv <- c('batch', 'genotype', 'exdate', 'age')

df.3 <- data.frame(matrix(ncol = 18))
colnames(df.3) <- c(AB,'Sample', 'batch', 'genotype', 'exdate', 'age')
count <- 0

#make a df with 9 entries (9 samples)
#take the mean of each sample (for categorical, there should only be 1 value) 
for (i in unique(df.2$Sample)) {
  count <- count+1
  for (j in AB) {df.3[count, j] <- mean(df.2[which(df.2$Sample == i), j])}
  df.3[count, 'Sample'] <- i
  for (j in iv){df.3[count, j] <- unique(df.2[which(df.2$Sample == i), j])}
}




melt_df.3 <- melt(df.3, measure.vars = AB, variable.name = 'antibody', value.name = 'AB.mean.expression') # variable = AB name, value = AB mean expression/sample


#### create a df with the mean for each cell type for each of the 9 samples
melt.df2 <- melt(df.2, measure.vars = AB, variable.name = 'antibody')

df.group <- melt.df2 %>% group_by(label, Sample, antibody) %>% 
  dplyr::summarize(Mean = mean(value, na.rm=TRUE))

### I do want the other variables to be added back
#batch
df.group[which(grepl("0317A", df.group$Sample)), "batch"] <- "A" #batch A (0317A = A)
df.group[which(grepl("0317B|0306", df.group$Sample)), "batch"] <- "B" #batch B (0317B ir 0306 = B)

#genotype
df.group[which(grepl("3450", df.group$Sample)), "genotype"] <- "3450" # genotype 3450
df.group[which(grepl("AIW002", df.group$Sample)), "genotype"] <- "AIW002" # genotype AIW002
df.group[which(grepl("AJG001C", df.group$Sample)), "genotype"] <- "AJG001C" # genotype AJG001C

#experiment date
df.group[which(grepl("0306", df.group$Sample)), "exdate"] <- "0306" # experiment date 0306
df.group[which(grepl("0317", df.group$Sample)), "exdate"] <- "0317" # experiment date 0317

#age
df.group[which(grepl("0306", df.group$Sample)), "age"] <- "273" # age 273, same as ex date 0306
df.group[which(grepl("0317B",df.group$Sample)), "age"] <- "284" # age 284
df.group[which(grepl("0317A", df.group$Sample)), "age"] <- "263" # age 263


head(df.group)
#### run 2 way anovas ###############

# loading the appropriate libraries
library(readr)
library(multcompView)


##### see the data 

df.cd24 <- df.group %>% filter(antibody == 'CD24')

qplot(x = as.factor(label), y = Mean, geom = "point", data = df.cd24) +
  facet_grid(.~genotype, labeller = label_both)

### This is still for all the cells

# you need to add * to get the interaction
res.aov2 <- aov(Mean ~ genotype * antibody, data = df.group)
summary(res.aov2)

# input which variable to do the invividual test on 
TukeyHSD(res.aov2, which = "antibody")
### this compares the AB with each other and is not helpful

# now this plots the interaction each AB @ a specific genotype
TukeyHSD(res.aov2)

# I'll need to group these results
# need to split 'genotype:antibody'

tukey.results <- TukeyHSD(res.aov2)

tukey.df <- tukey.results$`genotype:antibody`
head(tukey.df)



# great now I need to run this for each cell type 
A <- 'genotype'
B <- 'antibody'
# compare antibody expression between different variables
anv.l <- list()
tk.l <- list()
tk.l.sig <- list()
for (i in unique(df.group$label)){
  # subset the cell type
  df.group.sub <- df.group %>% filter(label == i)
  # run 2 way anova with interaction effect
  res.aov2 <- aov(Mean ~ antibody * genotype, data = df.group.sub)
  # print results and add results into a list
  print(paste("Cell type: ",i))
  print(summary(res.aov2))
  anv.l[[i]] <- summary(res.aov2)
  # now the posthoc tests
  tukey.results <- TukeyHSD(res.aov2)
  tukey.df <- as.data.frame(tukey.results$`antibody:genotype`)
  # save all results
  tk.l[[i]] <- tukey.results
  # filter for interactions that are significant
  tk.l.sig[[i]] <- tukey.df %>% filter(`p adj`<= 0.05)
  # see the significant comparisons
  print(tukey.df %>% filter(`p adj`<= 0.05))
}

# save results

saveRDS(anv.l, paste(output_path,"anova.geno.AB.RDS"))
saveRDS(tk.l, paste(output_path,"tukey.geno.AB.RDS"))
saveRDS(tk.l.sig, paste(output_path,"tukey.sig.geno.AB.RDS"))

# age and antibody for each cell type

A <- 'age'
B <- 'antibody'

anv.l <- list()
tk.l <- list()
tk.l.sig <- list()
for (i in unique(df.group$label)){
  # run 2 way anova with interaction effect
  # subset the cell type
  df.group.sub <- df.group %>% filter(label == i)
  res.aov2 <- aov(Mean ~ antibody * age, data = df.group.sub)
  # print results and add results into a list
  print(paste("Cell type: ",i))
  print(summary(res.aov2))
  anv.l[[i]] <- summary(res.aov2)
  # now the posthoc tests
  tukey.results <- TukeyHSD(res.aov2)
  tukey.df <- as.data.frame(tukey.results$`antibody:age`)
  # save all results
  tk.l[[i]] <- tukey.results
  # filter for interactions that are significant
  tk.l.sig[[i]] <- tukey.df %>% filter(`p adj`<= 0.05)
  # see the significant comparisons
  print(tukey.df %>% filter(`p adj`<= 0.05))
}


# save results

saveRDS(anv.l, paste(output_path,"anova.age.AB.RDS"))
saveRDS(tk.l, paste(output_path,"tukey.age.AB.RDS"))
saveRDS(tk.l.sig, paste(output_path,"tukey.sig.age.AB.RDS"))


##### compare batch

A <- 'batch'
B <- 'antibody'

anv.l <- list()
tk.l <- list()
tk.l.sig <- list()
for (i in unique(df.group$label)){
  # run 2 way anova with interaction effect
  # subset the cell type
  df.group.sub <- df.group %>% filter(label == i)
  res.aov2 <- aov(Mean ~ antibody * batch, data = df.group.sub)
  # print results and add results into a list
  print(paste("Cell type: ",i))
  print(summary(res.aov2))
  anv.l[[i]] <- summary(res.aov2)
  # now the posthoc tests
  tukey.results <- TukeyHSD(res.aov2)
  tukey.df <- as.data.frame(tukey.results$`antibody:batch`)
  # save all results
  tk.l[[i]] <- tukey.results
  # filter for interactions that are significant
  tk.l.sig[[i]] <- tukey.df %>% filter(`p adj`<= 0.05)
  # see the significant comparisons
  print(tukey.df %>% filter(`p adj`<= 0.05))
}


# save results


saveRDS(anv.l, paste(output_path,"anova.batch.AB.RDS"))
saveRDS(tk.l, paste(output_path,"tukey.batch.AB.RDS"))
saveRDS(tk.l.sig, paste(output_path,"tukey.sig.batch.AB.RDS"))


# experiment date

A <- 'batch'
B <- 'antibody'

anv.l <- list()
tk.l <- list()
tk.l.sig <- list()
for (i in unique(df.group$label)){
  # run 2 way anova with interaction effect
  # subset the cell type
  df.group.sub <- df.group %>% filter(label == i)
  res.aov2 <- aov(Mean ~ antibody * exdate, data = df.group.sub)
  # print results and add results into a list
  print(paste("Cell type: ",i))
  print(summary(res.aov2))
  anv.l[[i]] <- summary(res.aov2)
  # now the posthoc tests
  tukey.results <- TukeyHSD(res.aov2)
  tukey.df <- as.data.frame(tukey.results$`antibody:exdate`)
  # save all results
  tk.l[[i]] <- tukey.results
  # filter for interactions that are significant
  tk.l.sig[[i]] <- tukey.df %>% filter(`p adj`<= 0.05)
  # see the significant comparisons
  print(tukey.df %>% filter(`p adj`<= 0.05))
}


# save results and see the results again
anv.l
tk.l.sig
unique(df.group$label)

saveRDS(anv.l, paste(output_path,"anova.ex.AB.RDS"))
saveRDS(tk.l, paste(output_path,"tukey.ex.AB.RDS"))
saveRDS(tk.l.sig, paste(output_path,"tukey.sig.ex.AB.RDS"))


# make some plots to see what might be interesting
outpath <- "/Users/rhalenathomas/Documents/Projects_Papers/PhenoID/ForFigures/"

pdf(paste(output,"stats.boxplot.geno.ab.cell.pdf"),width = 11, height = 10)
ggplot(df.group, aes(x = label, y = Mean, colour = as.factor(genotype))) +
  geom_boxplot() +
  facet_wrap(~ antibody) + theme(axis.text.x = element_text(angle = 90))
dev.off()


# genotype
pdf(paste(outpath,"stats.boxplot.cell.geno.ab.pdf"),width = 11, height = 10)
ggplot(df.group, aes(x = antibody, y = Mean, colour = as.factor(genotype))) +
  geom_boxplot() +
  facet_wrap(~ label, scales="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# age
pdf(paste(outpath,"stats.boxplot.cell.age.ab.pdf"),width = 11, height = 10)
ggplot(df.group, aes(x = antibody, y = Mean, colour = as.factor(age))) +
  geom_boxplot() +
  facet_wrap(~ label, scales="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# batch
pdf(paste(outpath,"stats.boxplot.cell.geno.ab.pdf"),width = 11, height = 10)
ggplot(df.group, aes(x = antibody, y = Mean, colour = as.factor(batch))) +
  geom_boxplot() +
  facet_wrap(~ label, scales="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# exdate

pdf(paste(outpath,"stats.boxplot.cell.exdate.ab.pdf"),width = 11, height = 10)
ggplot(df.group, aes(x = antibody, y = Mean, colour = as.factor(exdate))) +
  geom_boxplot() +
  facet_wrap(~ label, scales="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# batch

pdf(paste(outpath,"stats.boxplot.cell.batch.ab.pdf"),width = 11, height = 10)
ggplot(df.group, aes(x = antibody, y = Mean, colour = as.factor(batch))) +
  geom_boxplot() +
  facet_wrap(~ label, scales="free") + theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#### check stats in a different way --- split first by AB then compare genotypes per cell type

anv.l <- list()
tk.l <- list()
tk.l.sig <- list()
for (i in unique(df.group$antibody)){
  # run 2 way anova with interaction effect
  # subset the cell type
  df.group.sub <- df.group %>% filter(antibody == i)
  res.aov2 <- aov(Mean ~ label * exdate, data = df.group.sub)
  # print results and add results into a list
  print(paste("AB: ",i))
  #print(summary(res.aov2))
  anv.l[[i]] <- summary(res.aov2)
  # now the posthoc tests
  tukey.results <- TukeyHSD(res.aov2)
  tukey.df <- as.data.frame(tukey.results$`label:exdate`)
  # save all results
  tk.l[[i]] <- tukey.results
  # filter for interactions that are significant
  tk.l.sig[[i]] <- tukey.df %>% filter(`p adj`<= 0.05)
  # see the significant comparisons
  #print(tukey.df %>% filter(`p adj`<= 0.05))
}


# save results and see the results again
anv.l
tk.l.sig


saveRDS(anv.l, paste(output_path,"anova.ex.label.RDS"))
saveRDS(tk.l, paste(output_path,"tukey.ex.label.RDS"))
saveRDS(tk.l.sig, paste(output_path,"tukey.sig.ex.label.RDS"))


### genotype

anv.l <- list()
tk.l <- list()
tk.l.sig <- list()
for (i in unique(df.group$antibody)){
  # run 2 way anova with interaction effect
  # subset the cell type
  df.group.sub <- df.group %>% filter(antibody == i)
  res.aov2 <- aov(Mean ~ label * genotype, data = df.group.sub)
  # print results and add results into a list
  print(paste("AB: ",i))
  #print(summary(res.aov2))
  anv.l[[i]] <- summary(res.aov2)
  # now the posthoc tests
  tukey.results <- TukeyHSD(res.aov2)
  tukey.df <- as.data.frame(tukey.results$`label:genotype`)
  # save all results
  tk.l[[i]] <- tukey.results
  # filter for interactions that are significant
  tk.l.sig[[i]] <- tukey.df %>% filter(`p adj`<= 0.05)
  # see the significant comparisons
  #print(tukey.df %>% filter(`p adj`<= 0.05))
}


# save results and see the results again
anv.l
tk.l.sig

unique(melt.df2$antibody)

saveRDS(anv.l, paste(output_path,"anova.gene.label.ab.RDS"))
saveRDS(tk.l, paste(output_path,"tukey.gene.label.ab.RDS"))
saveRDS(tk.l.sig, paste(output_path,"tukey.gene.gene.label.ab.RDS"))

### batch

anv.l <- list()
tk.l <- list()
tk.l.sig <- list()
for (i in unique(df.group$antibody)){
  # run 2 way anova with interaction effect
  # subset the cell type
  df.group.sub <- df.group %>% filter(antibody == i)
  res.aov2 <- aov(Mean ~ label * batch, data = df.group.sub)
  # print results and add results into a list
  print(paste("AB: ",i))
  #print(summary(res.aov2))
  anv.l[[i]] <- summary(res.aov2)
  # now the posthoc tests
  tukey.results <- TukeyHSD(res.aov2)
  tukey.df <- as.data.frame(tukey.results$`label:batch`)
  # save all results
  tk.l[[i]] <- tukey.results
  # filter for interactions that are significant
  tk.l.sig[[i]] <- tukey.df %>% filter(`p adj`<= 0.05)
  # see the significant comparisons
  #print(tukey.df %>% filter(`p adj`<= 0.05))
}


# save results and see the results again
anv.l
tk.l.sig

unique(melt.df2$antibody)

saveRDS(anv.l, paste(output_path,"anova.batch.label.ab.RDS"))
saveRDS(tk.l, paste(output_path,"tukey.batch.label.ab.RDS"))
saveRDS(tk.l.sig, paste(output_path,"tukey.batch.gene.label.ab.RDS"))

###age

anv.l <- list()
tk.l <- list()
tk.l.sig <- list()
for (i in unique(df.group$antibody)){
  # run 2 way anova with interaction effect
  # subset the cell type
  df.group.sub <- df.group %>% filter(antibody == i)
  res.aov2 <- aov(Mean ~ label * age, data = df.group.sub)
  # print results and add results into a list
  print(paste("AB: ",i))
  #print(summary(res.aov2))
  anv.l[[i]] <- summary(res.aov2)
  # now the posthoc tests
  tukey.results <- TukeyHSD(res.aov2)
  tukey.df <- as.data.frame(tukey.results$`label:age`)
  # save all results
  tk.l[[i]] <- tukey.results
  # filter for interactions that are significant
  tk.l.sig[[i]] <- tukey.df %>% filter(`p adj`<= 0.05)
  # see the significant comparisons
  #print(tukey.df %>% filter(`p adj`<= 0.05))
}


# save results and see the results again
anv.l
tk.l.sig


saveRDS(anv.l, paste(output_path,"anova.age.label.ab.RDS"))
saveRDS(tk.l, paste(output_path,"tukey.age.label.ab.RDS"))
saveRDS(tk.l.sig, paste(output_path,"tukey.age.label.ab.RDS"))


###### note - maybe pair wise comparisons with the BH or other correction would work better for the post hoc test.






















#preprocessing for 2 way anova ends




#2. 2 way anova function starts
twanova <- function(df, A, B, dv) { #B is a list
  #result dfs: anova df and post-hoc df
  aov_table <- data.frame(matrix(ncol = 5, nrow = 0))
  colnames(aov_table) <- c("A", "B", "pval.A", "pval.B", "pval.interaction")
  ph_table <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(ph_table) <- c("type", "comparison", "pval")
  
  for (i in B) {
    aov_result <- aov(df[,dv] ~
                        df[,i] + df[,A] + df[,i]:df[,A])
    aov_result <-  summary(aov_result)[[1]][["Pr(>F)"]] 
    #anova test: dv (AB mean expression level) ~ iv B (Antibody) +
    #iv A (one of the age, batch...) + iv B: iv A (interaction effect)
    
    aov_table <-rbind(aov_table, list(
      A = A, B = i, pval.A = aov_result[1], 
      pval.B = aov_result[1], pval.interaction = aov_result[3]))
    
    #post hoc for A or B if significant: interaction effect is excluded from ph
    if ((aov_result[1] < 0.05) & 
        (aov_result[2] < 0.05)) {
      tkwhich <- c('df[, A]','df[, i]') #if A & B are significant
    } else if (aov_result[1] < 0.05) {
      tkwhich <- 'df[, A]' #if only A is significant
    } else if (aov_result[2] < 0.05) {
      tkwhich <- 'df[, i]' #if only B is significant
    } else { next }
    
    tk <- TukeyHSD(aov_result, which = tkwhich, conf.level=.95) #Tukey test
    
    for(j in names(tk)) {
      #temporary list, type = A or B main effect
      tkl <- list(type = rep(if (j=="df[, A]") A else if (j=="df[, i]") i, 
                             length(rownames(tk[[j]]))), 
                  comparison = rownames(tk[[j]]),
                  pval = as.numeric((tk[[j]])[,'p adj'])
      )
      ph_table <- rbind(ph_table, tkl)
      #shuming note: probably not good to use rbind in a for loop, change later
    }
  }
  return(list(aov = aov_table, ph = ph_table))
}

#3. input for the function
df <- melt_df.3
A <- 'antibody'
B <- c('batch', 'genotype', 'exdate', 'age')
dv <- 'AB.mean.expression'


#4. test the function
test <- twanova(df, A, B, dv)
test$aov
test$ph


#### this appears to test only the mean values from all the cells not individual cell types



# I didn't both to save the results

output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/ANOVA_CHI/"
write.csv(test$aov, paste(output_path, "aov_test.csv",sep=""), row.names = FALSE)
saveRDS(aov, paste(output_path, "aov_and_ph.Rds",sep=""))

# write.csv(ph_table, paste(output_path, "ph_table.csv",sep=""), row.names = FALSE)

##################### proportion tests #####################################################

#1. function
proptest <- function(test=c("prop1", "prop2", "chisq1"),
                     df,
                     c, #c = matrix column, ex: cell types, the name of the column in og df
                     r){ #r = matrix row, ex: genotype, the name of the column in og df
  
  #preprocessing:
  testl <- vector() #a list of cell counts 
  for (i in unique(df[, c])) { 
    for (j in unique(df[, r])) {
      testl <- c(testl, length(which((df[, c] == i) & (df[, r] == j))))
    }
  }
  
  #create a matrix: col = cell type, r = batch in this case
  testm <- matrix(testl,ncol = length(unique(df[, r])), byrow = TRUE)
  
  #proportion test #1, compare 1 batch x 1 cell type at a time
  #if input contains prop1, this test will run
  if ("prop1" %in% test) { 
    pl1 <- vector()
    for (i in 1:nrow(testm)) {  
      for (j in 1:ncol(testm)) {
        prop <- prop.test(x = testm[i, j], # x=the cell count tested
                          n = sum(testm[i,]), # n=sum of the row
                          p = 1/ncol(testm)) #p=probability in each group (genotype) 
        pl1 <- c(pl1, prop$p.value) #p value list
      }
    }
    pd1 <- data.frame(c = rep(unique(df[,c]), length(unique(df[,r]))),
                      r = rep(unique(df[,r]), length(unique(df[,c]))),
                      pval = pl1)
  } else {pd1 <- NULL}
  
  #proportion test #2, compare 3 batchs x 1 cell type at a time
  if ("prop2" %in% test) {
    hp2 <- function(x) { #helper function
      #x=each row in the matrix, n=sum of the row 
      prop <- prop.test(x, rep(sum(x), length(x)), p=rep(1/length(x),length(x)))
      return(prop$p.value)
    }
    pd2 <- data.frame(c = unique(df[,c]), pval = apply(testm,1,hp2)) 
  } else {pd2 <- NULL}
  
  #chi square test, compare 3 bacths x 1 cell type at a time
  if ("chisq1" %in% test) {
    hp3 <- function(x) {
      prop <- chisq.test(x)
      return(prop$p.value)
    } 
    pd3 <- data.frame(c = unique(df[, c]), pval = apply(testm, 1, hp3))
  } else {pd3 <- NULL}
  return(list(prop1 = pd1, prop2 = pd2, chisq1 = pd3))
}


### run the test  ####
### inputs

test = c("prop1", "prop2", "chisq1") #run which tests? options: c("prop1", "prop2", "chisq1")
df = df.2 #see preprocessing at the very beginning of this R script
c = "label" #cell type
r = "genotype" #we want to test the difference among batches

# 3. test

dfl <- proptest(test, df, c, r)

prop1.df <- dfl$prop1
prop2.df <- dfl$prop2
chi.df <- dfl$chisq1


dfex <- proptest(test, df, 'label', 'age')
dfex$prop1
