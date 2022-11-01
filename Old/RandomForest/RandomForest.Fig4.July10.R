library(randomForest)
library(caret)
library(Seurat)
library(data.table)



# query data
seu.r<- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/Seu9000annot.08072021.RDS")
AB <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")


df <- transpose(as.data.frame(GetAssayData(seu.r,slot = 'scale.data')))
dim(df)

# the data frame should be in the order I made the df but the column names are not the AB names
colnames(df)

colnames(df) <- AB
print(colnames(df))
# add in the cell labels

# I will use the smallest group labels 
annotations <- seu.r$subgroups

df.l <- cbind(df, lables = annotations)
print(table(df.l$lables))

# split the data
set.seed(222)
ind <- sample(2, nrow(df.l), replace = TRUE, prob = c(0.5, 0.5)) # prop is the proportions
# split the data into train and test
train <- df.l[ind==1,]
test <- df.l[ind==2,]


# from the https://www.guru99.com/r-random-forest-tutorial.html#3 
# there are paramaters to tune:
# ntree: number of trees in the forest
# mtry: number of canditates to draw - defualt is the square of the number of columns
# maxnodes: max turminal nodes

# ML optomization : Random search and Grid search
#Grid Search definition
#The grid search method is simple, the model will be evaluated over all the combination you pass in the function, using cross-validation.
# k fold cross validation:

#trainControl(method = "cv", number = n, search ="grid")
# method of resampling
# number - of k fold
# search : grid or random


# run default model
rf <- randomForest(lables~., data=test, proximity=TRUE, ntrees = 50)  

trControl <- trainControl(method = "cv", number = 10, search ="grid")

#You will use caret library to evaluate your model. The library has one function called train() to evaluate almost all machine learning algorithm.


#train(formula, df, method = "rf", metric= "Accuracy", trControl = trainControl(), tuneGrid = NULL)


# run default model
rf.default <- randomForest(lables~., data=train, method = "rf", metric = "Accuracy", trControl = trControl) 
# print results
print(rf.default)

# now try tuning
set.seed(48)

tuneGrid <- expand.grid(.mtry = c(1: 10))
rf_mtry <- train(lables~., 
                 data=train,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 14,
                 ntree = 300)
print(rf_mtry)

# the mtry with the best accuracy was mtry = 6
rf_mtry$bestTune$mtry
max(rf_mtry$results$Accuracy)

best_mtry <- rf_mtry$bestTune$mtry 
best_mtry

# search the best node size
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(12: 25)) {
  set.seed(1234)
  rf_maxnode <- train(lables~.,
                      data = train,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  current_iteration <- toString(maxnodes)
  store_maxnode[[current_iteration]] <- rf_maxnode
}
results_mtry <- resamples(store_maxnode)
summary(results_mtry)

# node size 25 is best
# higher values might give better accuracy

store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(20: 30)) {
  set.seed(1234)
  rf_maxnode <- train(lables~.,
                      data = train,
                      method = "rf",
                      metric = "Accuracy",
                      tuneGrid = tuneGrid,
                      trControl = trControl,
                      importance = TRUE,
                      nodesize = 14,
                      maxnodes = maxnodes,
                      ntree = 300)
  key <- toString(maxnodes)
  store_maxnode[[key]] <- rf_maxnode
}
results_node <- resamples(store_maxnode)
summary(results_node)

# max node = 30


# now search of the best number of trees
store_maxtrees <- list()
for (ntree in c(300, 400, 500, 600, 800, 1000, 2000)) {
  set.seed(5678)
  rf_maxtrees <- train(lables~.,
                       data = train,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 25,
                       maxnodes = 30,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
# best number of trees is different between media and max 
# 2000 best for max 300 best for mean and median
# go with 1000 close to highest max and better mean and median

# now fit the model with the best conditions
#  mtry = 6
# maxnodes = 30
# node size = 25
# ntree = 1000

rf <- randomForest(lables~.,
                   train,
                   mtry = 6,
                   nodesize = 25,
                   ntree = 1000,
                   maxnodes = 30)

# save outputs
saveRDS(rf, output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/RFM_trained.11072022.Rds")



##### main function finishes here ########################
##### I think we want some kind of print out of the model accuracy


# helper functions:  
# 1. maker confusion matrix
# 2. print out accuracy in the training text

# function two - predict cell types unlabelled from a given cell type and keep the results.

# check the accuracy of the model 

prediction.train <-predict(rf, train)
prediction.test <-predict(rf, test)

### check the model accuracy 
# to do so we need to run the data 

print("predict training data")
confusionMatrix(prediction.train, train$lables)

# this makes a df of the frequency of predictions that can make the confusion matrix
output.df <- as.data.frame(table(prediction.train, train$lables))


print("predict test data")
confusionMatrix(prediction.test, test$lables)

# there are also overall statistics to print out that it would be better to change into a different format
# maybe we just print to log
# better to plot the confusion matrix than print it out



#Overall Statistics

#Accuracy : 0.6613          
#95% CI : (0.6564, 0.6661)
#No Information Rate : 0.1369          
#P-Value [Acc > NIR] : < 2.2e-16       
#Kappa : 0.629    



# now the predictions
p1 <- predict(fit_rf, train)
confusionMatrix(p1, train$lables)

p2 <- predict(fit_rf, test)
c2 <- confusionMatrix(p2, test$lables)
c2.table <- as.data.frame(c2$table)

# try to plot the results

library(ggplot2)
ggplot(data =  c2.table, mapping = aes(x = Prediction, y = Reference, fill= Freq)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f",Freq)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") + theme(axis.text.x = element_text(angle = 90))



output_path <- "/Users/rhalenathomas/Documents/Projects_Papers/PhenoID/ForFigures/"
pdf(paste(output_path, "RFMconfustion.Test.results"))
ggplot(data =  c2.table, mapping = aes(x = Prediction, y = Reference, fill= Freq)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f",Freq)), vjust = 1) +
  scale_fill_gradient(low = "grey92", high = "mediumpurple4") + theme(axis.text.x = element_text(angle = 90))
dev.off()

# save the model for later

saveRDS(fit_rf, output_path <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/RFM_trained.11072022.Rds")
write.csv(c2.table, "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/PaperFigures/RFMtest.results.11072022.csv")


##### another blog that might be better for visualizing tuning parameters
# https://www.r-bloggers.com/2016/07/random-forests-in-r-2/
  
########################## function/s to see how the model is performing

# vignette
# https://cran.rstudio.com/web/packages/randomForestExplainer/vignettes/randomForestExplainer.html

devtools::install_github("MI2DataLab/randomForestExplainer")
library(randomForestExplainer)

min_depth_frame <- min_depth_distribution(rf)
head(min_depth_frame, n = 10)

plot_min_depth_distribution(min_depth_frame)

# not sure how to interpret the above

# more importance measures
importance_frame <- measure_importance(rf)

# plot_multi_way_importance(forest, size_measure = "no_of_nodes") # gives the same result as below but takes longer
plot_multi_way_importance(importance_frame, size_measure = "no_of_nodes")
# low mean_min_depth and high times_a_root are the most important factors

plot_multi_way_importance(importance_frame, 
                          x_measure = "mean_min_depth", 
                          y_measure = "no_of_nodes",
                          size_measure = "p_value", no_of_labels = 6)

plot_multi_way_importance(importance_frame, 
                          x_measure = "mean_min_depth", 
                          y_measure = "gini_decrease",
                          size_measure = "p_value", no_of_labels = 13)


#### try to plot error rate
oob.err.data1 <- data.frame(
  Trees = rep(1:nrow(rf$err.rate), 20), 
  Type = rep(c("OOB","Unknown","Mixed","Neurons 1","Neurons 2","Neurons 3",
               "NPC","Neural stem","Radial Glia 1","Radial Glia 2","Radial Glia 3", 
               "Epithelial","Endothelial","Oligodendrocytes",
               "Astrocytes 1","Astrocytes 2","Astrocytes 3","Astrocytes mature",
               "Stem-like 1","Stem-like 2"), each = nrow(rf$err.rate)),
  Error = c(rf$err.rate[,"OOB"], rf$err.rate[,"Unknown"], rf$err.rate[,"Mixed"], 
            rf$err.rate[,"Neurons 1"],rf$err.rate[,"Neurons 2"],rf$err.rate[,"Neurons 3"],
            rf$err.rate[,"NPC"],rf$err.rate[,"Neural stem"],rf$err.rate[,"Radial Glia 1"],
            rf$err.rate[,"Radial Glia 2"],rf$err.rate[,"Radial Glia 3"],
            rf$err.rate[,"Epithelial"],rf$err.rate[,"Endothelial"],rf$err.rate[,"Oligodendrocytes"],
            rf$err.rate[,"Astrocytes 1"],rf$err.rate[,"Astrocytes 2"],rf$err.rate[,"Astrocytes 3"],
            rf$err.rate[,"Astrocytes mature"],rf$err.rate[,"Stem-like 1"],rf$err.rate[,"Stem-like 2"]
            ))

# this works but we need a way to not have to type in all the cell type names

ggplot(data = oob.err.data1, aes(x = Trees, y= Error)) + geom_line(aes(color = Type))


############################# run the random forest on the unlabelled data ############
seu.q <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outs/clusters/Louvain/Allcellsretros_LouvainSeuratObject60.Rds")

df <- transpose(as.data.frame(GetAssayData(seu.q,slot = 'scale.data')))
dim(df)
colnames(df) <- AB
rfm.pred <- predict(rf,df)
head(rfm.pred)
results.df <- as.data.frame(rfm.pred)
head(results.df)

# save the results
write.csv(results.df, paste(outpath, "RFMpredictionsAllcells.12072022.csv"))
