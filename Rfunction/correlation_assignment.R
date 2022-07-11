

library(dplyr) #rename
library(kit) # for finding max and second max (function topn) 
library(reshape2) #for melt?

# ========= preprocessing for all functions in this file starts =========

# 1. corr function
find_correlation <- function(test, 
                             reference, 
                             output_path = NULL, #if null, don't save 
                             min_corr = 0.1, 
                             min_diff = 0.05) {

  #1. process and scale test and reference matrix
  
  #find intersect between reference and test, ignore case, 
  #change reference spelling, order to test's
  testc <- vector()
  refc <- vector()
  for (i in colnames(test)) {
    for (j in colnames(reference)) {
      if (tolower(i) == tolower(j)) {
        testc <- c(testc, i)
        refc <- c(refc, j)
      }
    }
  } 
  
  reference <- reference[, refc] #select markers + X in reference
  colnames(reference) <- testc #change reference's spelling to test's
  test <- test[, testc] #select markers + X in test
  markers <- colnames(select_if(test, is.numeric)) #a list of markers (without X)
  test[, markers] <- scale(test[, markers]) #z score test (without X) 
  reference[, markers] <- scale(reference[, markers])  # z score the reference matrix 
  
  # test <- test[sample(1:nrow(test), 1000),] #testing with 5 samples
  
  #2. find best and second correlation and cell type
  result <- vector()

  #the loop that will find the best and second best correlation and cell types
  for (i in 1:nrow(test)) {
    corr_ls <- vector() #list of correlation between the reference cell types and each sample
    for (j in 1:nrow(reference)) {
      corr <- cor(as.numeric(test[i,markers]), as.numeric(reference[j,markers])) # pearson by default and we use default
      corr_ls <- c(corr_ls, corr)
    }
    
    top <- topn(corr_ls, 2L, decreasing = TRUE) #return the index of the best 2
    result <- c(result, 
                test[i, 'X'], #col 1: cell sample
                corr_ls[top[1]], #col 2: 1st correlation
                reference[top[1], 'X'], #col 3: 1st best cell type
                corr_ls[top[2]], #col 4: 2nd correlation 
                reference[top[2], 'X'], #col 5: 2nd best cell type
                ifelse(corr_ls[top[1]] < min_corr, "unknown",
                       ifelse(corr_ls[top[1]] - corr_ls[top[2]] < min_diff, 
                              paste(reference[top[1], 1], 
                                    reference[top[2], 1], sep = "-"), 
                              reference[top[1], 1]))) #col 6: assigned cell type
    # if best corr < min_corr, assign unknown cell type
    # if best corr - second best corr < min diff, assign combined cell type 
    # else, assign best cell type
  
  #convert the result list to a df  
  cdf <- data.frame(matrix(result, ncol=6, byrow = TRUE))
  colnames(cdf) <- c("X", "cor.1", "best.cell.type",
                       "cor.2", "second.cell.type", "cell.label")
  cdf$cor.1 <- as.numeric(cdf$cor.1)
  cdf$cor.2 <- as.numeric(cdf$cor.2)
  return(cdf)
  }
}

# 2. input
test <- read.csv("/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/Old/preprocessing/outputs/prepro_outstransformed_flowset.csv")
reference <- read.csv("/Users/shumingli/Documents/GitHub/PhenoID_single_cell_flow_cytometry_analysis/Old/correlation/ReferenceMatrix9celltypesOrdered.csv")

#fix reference: replace NA in [epithelial, 04] with the mean expression of 04 in other cell types
reference[9,"O4"] <- mean(reference[1:8, "O4"])

# 3. test
x <- find_correlation(test, reference)


# ========================= plot starts ============================


plot_corr <- function(df) {
  # filter to get frequency table and save as csv
  df.f <- df %>% select(cell.label)
  freq.table <- as.data.frame(table(df.f))
  df.filter <- df %>% group_by(cell.label) %>% dplyr::filter(n()> 100)
  
  # plot the frequencies
  plot1 <- 
    ggplot(df.filter, aes(x = reorder(
      cell.label, cell.label, function(x) - length(x)),
      fill = cell.label)) + geom_bar() + theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) + 
    xlab('Assigned cell type') + 
    ylab('number of cell') + 
    labs(fill='Cell Types')
  # print(plot1)
  
  # violin plot of best correlation/cell type
  plot2 <- 
    ggplot(df, aes(x = best.cell.type, y = cor.1)) +
    geom_violin()+ 
    ylim(-0.1, 1)+
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90)) + 
    ylab("correlation coefficient") + 
    xlab("Cell type with max correlation coefficient")
  # print(plot2)
  
  df.melt <- melt(df) #reformat to long df
  
  # plot the best and second best correlation together
  plot3 <- 
    ggplot(df.melt, aes(x = cell.label, y = value ))+ 
    geom_boxplot()+ ylim(-0.1, 1)+theme_classic()+
    theme(axis.text.x = element_text(angle = 90))+ 
    ylab("correlation coefficient") + 
    xlab("Cell type label")
  # print(plot3)
  
  # plot the best and second best correlation separated on the same graph
  plot4 <- 
    ggplot(df.melt, aes(x = best.cell.type, y = value, fill = variable))+ 
    geom_boxplot()+ 
    ylim(-0.25, 1)+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90)) + 
    scale_fill_manual(values = c("#4E84C4", "#52854C")) + 
    ylab("correlation coefficient") + xlab("Cell type")
  # the second best correlation is so low it was removed from STEM with axis limit -0.1 and even -1
 
  # print(plot4)
  
  # down sample
  set.seed(64)
  df.downsample <- sample_n(df, 1000)
  df.melt.down <- melt(df.downsample)
  
  # # reformat the table to work with the before after plot
  # # y is the measurement in df.melt = value
  # # x is before after in df.melt = variable
  # # class another variable - in the example this is different shapes - for us this is best cell type
  # # might use facet to split the cell type - needs to be a factor
  # # id is the individual id this is the X column
  plot5 <- 
    ggplot(df.melt.down, aes(x = variable, y = value,colour = variable, group = X)) +
    geom_line(show.legend = F, size = 0.1, color = "black") + 
    geom_point()+ 
    scale_color_manual(values = c("#4E84C4", "#52854C")) + 
    ylim(-0.25, 0.95) +
    facet_wrap(~(as.factor(best.cell.type))) +
    theme(legend.position = "none") +
    ylab("Correlation Coefficient") +
    xlab("")
  # print(plot5)
  
  double.cells <- df[grep("-", df$cell.label),]
  df.melt.double <- melt(double.cells)
  
  # # this will be an excellent visualization but I need to subset only the double labels, then I can plot more cells and see more clearly.
  plot6 <- 
    ggplot(df.melt.double, aes(x = variable, y = value,colour= variable, group= X)) +
    geom_line(show.legend = F, size = 0.1, color = "black") + 
    geom_point()+ 
    scale_color_manual(values = c("#4E84C4", "#52854C")) + 
    ylim(-0.15,0.8) +
    facet_wrap(~(as.factor(cell.label))) +
    ylab("Correlation Coefficient") +
    xlab("")
  # print(plot6)
  
  return(list(freq.table, plot1, plot2, plot3, plot4, plot5, plot6))
}

#2. test 
y <- plot_corr(x)

# ========================= plot ends ============================

