library(dplyr) #used to rename melted df

# === preprocessing for all functions in this file starts ===

#this step is to create columns of features such as genotype, experiment date, 
#age and batch. no need to run this if users have these categories in input df

#in order to run two way anova, users should first test the assumptions (on their own)

input_path <- "/Users/shumingli/Downloads/AllcellLablesMarch25.Rds"
input <- readRDS(input_path)
df <- transpose(as.data.frame(GetAssayData(input,slot = 'scale.data')))
col.names <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")
colnames(df) <- col.names 
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


#  ========================= ANOVA starts =========================

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
} #shuming notes: maybe replace with apply later

melt_df.3 <- melt(df.3) # variable = AB name, value = AB mean expression/sample
melt_df.3 <- rename(melt_df.3, antibody = variable)

#preprocessing for 2 way anova ends


#2. 2 way anova function starts
twanova <- function(df, A, B, dv) { #B is a list
  
  #result dfs: anova df and post-hoc df
  aov_table <- data.frame(matrix(ncol = 5))
  colnames(aov_table) <- c("A", "B", "pval.A", "pval.B", "pval.interaction")
  count <- 0 #row number, update every iteration in the for loop
  ph_table <- data.frame(matrix(ncol = 3))
  colnames(ph_table) <- c("type", "comparison", "pval")
  ctkl <- ls()
  
  #shuming note: will replace for loop with apply later
  for (i in B) {
    aov_result <- aov(df[,dv] ~
                        df[,i] + df[,A] + df[,i]:df[,A])
  #anova test: dv (AB mean expression level) ~ iv B (Antibody) +
    #iv A (one of the age, batch...) + iv B: iv A (interaction effect)

    count <- count + 1
    aov_table[count, "A"] <- 'antibody'
    aov_table[count, "B"] <- i
    aov_table[count, "pval.A"] <- summary(aov_result)[[1]][["Pr(>F)"]][1]
    aov_table[count, "pval.B"] <- summary(aov_result)[[1]][["Pr(>F)"]][2]
    aov_table[count, "pval.interaction"] <- summary(aov_result)[[1]][["Pr(>F)"]][3]
    
    #post hoc for A or B if significant: interaction effect is excluded from ph
    if ((summary(aov_result)[[1]][["Pr(>F)"]][1] < 0.05) & 
        (summary(aov_result)[[1]][["Pr(>F)"]][2] < 0.05)) {
      tkwhich <- c('df[, A]','df[, i]') #if A & B are significant
    } else if (summary(aov_result)[[1]][["Pr(>F)"]][1] < 0.05) {
      tkwhich <- 'df[, A]' #if only A is significant
    } else if (summary(aov_result)[[1]][["Pr(>F)"]][2] < 0.05) {
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
  #remove a row if the whole row is NAs
  ph_table <- filter(ph_table, rowSums(is.na(ph_table)) != ncol(ph_table)) 
  
  return(list(aov = aov_table, ph = ph_table))
}
#2 way anova function ends

#3. input for the function
df <- melt_df.3
A <- 'antibody'
B <- c('batch', 'genotype', 'exdate', 'age')
dv <- 'value'

#4. test the function
test <- twanova(df, A, B, dv)

# output_path <- "/Users/shumingli/Desktop/output_jun19/"
# write.csv(aov_table, paste(output_path, "aov_table.csv",sep=""), row.names = FALSE)
# write.csv(ph_table, paste(output_path, "ph_table.csv",sep=""), row.names = FALSE)

#  ========================= ANOVA ends =========================



#  ========================= Chi square starts =========================



#  ========================= Chi square ends =========================

