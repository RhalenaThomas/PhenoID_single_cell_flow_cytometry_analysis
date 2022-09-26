# library(dplyr) #for
library(data.table) #for transpose()
library(reshape2)#used to rename melted df


# === preprocessing for all functions in this file starts ===

#this step is to create columns of features such as genotype, experiment date, 
#age and batch. no need to run this if users have these categories in input df

#in order to run two way anova, users should first test the assumptions (on their own)

input_path <- "/Users/shumingli/Downloads/AllcellLablesMarch25.Rds"
input <- readRDS(input_path)
df <- transpose(as.data.frame(GetAssayData(input,slot = 'scale.data')))
AB <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")
colnames(df) <- AB 
df.2 <- cbind(df, label = input@meta.data$cluster.ids, Batch = input$Batch)

#add new columns batch, genotype, experiment day and age based on Batch (the old column)

#batch
df.2[which(grepl("0317A", df.2$Batch)), "batch"] <- "A" #batch A (0317A = A)
df.2[which(grepl("0317B|0306", df.2$Batch)), "batch"] <- "B" #batch B (0317B ir 0306 = B)

#genotype
df.2[which(grepl("3450", df.2$Batch)), "genotype"] <- "3450" # genotype 3450
df.2[which(grepl("AIW002", df.2$Batch)), "genotype"] <- "AIW002" # genotype AIW002
df.2[which(grepl("AJG001C", df.2$Batch)), "genotype"] <- "AJG001C" # genotype AJG001C

#experiment date
df.2[which(grepl("0306", df.2$Batch)), "exdate"] <- "0306" # experiment date 0306
df.2[which(grepl("0317", df.2$Batch)), "exdate"] <- "0317" # experiment date 0317

#age
df.2[which(grepl("0306", df.2$Batch)), "age"] <- "273" # age 273, same as ex date 0306
df.2[which(grepl("0317B", df.2$Batch)), "age"] <- "284" # age 284
df.2[which(grepl("0317A", df.2$Batch)), "age"] <- "263" # age 263

# === preprocessing for all functions ends ===
#df.2 now has 4 extra columns: batch, age...


#  ========================= I. ANOVA starts =========================

#1. preprocessing for 2 way anova starts (can be used for any anova)
#this preprocessing stacks 13 columns of antibody expression into 1 column 

AB <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")
iv <- c('batch', 'genotype', 'exdate', 'age')



df.3 <- data.frame(matrix(ncol = 18))
colnames(df.3) <- c(AB,'Batch', 'batch', 'genotype', 'exdate', 'age')
count <- 0

#make a df with 9 entries (9 samples)
#take the mean of each sample (for categorical, there should only be 1 value) 

for (i in unique(df.2$Batch)) {
  count <- count+1
  for (j in AB) {df.3[count, j] <- mean(df.2[which(df.2$Batch == i), j])}
  df.3[count, 'Batch'] <- i
  for (j in iv){df.3[count, j] <- unique(df.2[which(df.2$Batch == i), j])}
}




melt_df.3 <- melt(df.3, measure.vars = AB, variable.name = 'antibody', value.name = 'AB.mean.expression') # variable = AB name, value = AB mean expression/sample


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
aov <- twanova(df, A, B, dv)


output_path <- "/Users/shumingli/Desktop/output_jul6/"
write.csv(test$aov, paste(output_path, "aov_test.csv",sep=""), row.names = FALSE)
saveRDS(aov, paste(output_path, "aov_and_ph.Rds",sep=""))

# write.csv(ph_table, paste(output_path, "ph_table.csv",sep=""), row.names = FALSE)

#  ========================= ANOVA ends =========================



#  ================= II. Chi square & proportion test starts ===================

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

# 2. input
test = c("prop1", "prop2", "chisq1") #run which tests? options: c("prop1", "prop2", "chisq1")
df = df.2 #see preprocessing at the very beginning of this R script
c = "label" #cell type
r = "genotype" #we want to test the difference among batches

# 3. test

dfl <- proptest(test, df, c, r)


output_path <- "/Users/shumingli/Desktop/output_jul6/"
saveRDS(dfl, paste(output_path, "prop_and_chisq.Rds",sep=""))

# ====================== chi square & proportion test ends =====================

saveRDS(melt_df.3, paste(output_path, "aov_input_data.Rds",sep=""))
saveRDS(df2, paste(output_path, "prop_input_data.Rds",sep=""))



