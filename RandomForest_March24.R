library(randomForest)
library(caret)
library(Seurat)
library(data.table)



# add some data
seu <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/Louvkn60DifferentCellLabels220220318.Rds")
df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))

dim(df)


# I need to know the names of the columns now
# from the input df
col.names <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")

colnames(df) <- col.names
print(col.names)
# add in the cell lables

annotations <- seu$labels6

df.2 <- cbind(df, lables = annotations)
print(table(df.2$lables))

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

# not try tuning
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

# the mtry with the best accuracy was mtry = 5

rf_mtry$bestTune$mtry
max(rf_mtry$results$Accuracy)

best_mtry <- rf_mtry$bestTune$mtry 
best_mtry

# search the best maxnodes
store_maxnode <- list()
tuneGrid <- expand.grid(.mtry = best_mtry)
for (maxnodes in c(5: 15)) {
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

# max node best accuracy = 14 

# I won't run the higher node values

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

# best max node 29
# node size 15 was best

# now search of the best number of trees
store_maxtrees <- list()
for (ntree in c(250, 300, 350, 400, 450, 500, 550, 600, 800, 1000, 2000)) {
  set.seed(5678)
  rf_maxtrees <- train(lables~.,
                       data = train,
                       method = "rf",
                       metric = "Accuracy",
                       tuneGrid = tuneGrid,
                       trControl = trControl,
                       importance = TRUE,
                       nodesize = 15,
                       maxnodes = 30,
                       ntree = ntree)
  key <- toString(ntree)
  store_maxtrees[[key]] <- rf_maxtrees
}
results_tree <- resamples(store_maxtrees)
summary(results_tree)
# best number of trees is 2000

# now fit the model with the best conditions
#  mtry = 5
# maxnodes = 14
# node size is min nodes lets us
# ntree = TBD

# results from test data
#Accuracy : 0.7632          
#95% CI : (0.7588, 0.7676)
#No Information Rate : 0.1859          
#P-Value [Acc > NIR] : < 2.2e-16       
#Kappa : 0.7262  

# varImpPlot(fit_rf) # can't run because the model is not class randomForst but train.formula

# train with RandomForest and save the model
rf <- randomForest(lables~.,
                   train,
                   mtry = 5,
                   nodesize = 10,
                   ntree = 2000,
                   maxnodes = 15)

# now the predictions
p1 <- predict(rf, train)
confusionMatrix(p1, train$lables)

p2 <- predict(rf, test)
c2 <- confusionMatrix(p2, test$lables)
c2.table <- as.data.frame(c2$table)

# try to plot the results

library(ggplot2)
ggplot(data =  c2.table, mapping = aes(x = Prediction, y = Reference, fill= Freq)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f",Freq)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") + theme(axis.text.x = element_text(angle = 90))


# save the model for later - model with the firts annotation

saveRDS(rf, "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/rf_trained_march15.Rds")
write.csv(c2.table, "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/confusionm_values_march15.csv")



# model with new annotation

saveRDS(rf, "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/rf_trained_march25.Rds")
write.csv(c2.table, "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/confusionm_values_march25.csv")


rf.old <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/rf_trained_march15.Rds")


# make a predict and plot function 
# input arguments: the rf model, the data to test
# outputs : predictions and heatmap of results

#not test data must have column labels 
# this will only work on the input dataset spllit into training and test

rf_confusion <- function(rf_model, test_data) {
  p2 <- predict(rf_model, test_data)
  c2 <- confusionMatrix(p2, test_data$lables)
  p.table <- as.data.frame(c2$table)
  print(p.table)
  
  p <- ggplot(data =  p.table, mapping = aes(x = Prediction, y = Reference, fill= Freq)) +
    geom_tile(aes(fill = Freq), colour = "white") +
    geom_text(aes(label = sprintf("%1.0f",Freq)), vjust = 1) +
    scale_fill_gradient(low = "white", high = "steelblue") + theme(axis.text.x = element_text(angle = 90))
  print(p)
  
  
}





rf_confusion(rf.old, test.old)
# rf_confusion(rf.old, test) # doesn't work - the data levels in RF and test have to match ... try to fix
# later


# try to make predctions on unlableled data 

rf_model = rf
test_data = test
test_data <- test_data

pred <- predict(rf_model, test_data)
print(pred)


# this make preditions but only how many cells in each 
# try to get prediction by cells


c3 <- as.data.frame(pred) # this give a vector of predictions 

dim(test_data)
dim(c3)
labels <- as.data.frame(test_data$lables)
dim(labels)

compares <- cbind(labels,c3)



rf_predictions <- function(rf_model, test_data, outpath,filename) {
  pred1 <- predict(rf_model, test_data)
  print(pred1)
  write.csv(pred1, paste(outpath,filename))

}

# need to run the predictions on the full input data to test out visualizing on UMAP
outpath <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/"
filename <- "prediction_march25rf_dl2.csv"
  
rf_predictions(rf, df.2, outpath,filename)

filename <- "prediction_march15rf_dl2.csv"

rf_predictions(rf.old, df.2, outpath,filename)



##### first round of predicting the old dataset


# add some data
seu <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outs/clusters/Louvain/Allcellsretros_LouvainSeuratObject200.Rds")
df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))

dim(df)


# I need to know the names of the columns now
# from the input df
col.names <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")

colnames(df) <- col.names
print(col.names)
# add in the cell lables


outpath <- "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/"
filename <- "prediction_march15rf_allcells.csv"

rf_predictions(rf.old, df, outpath,filename)

filename <- "prediction_march25rf_allcells.csv"

rf_predictions(rf, df, outpath,filename)


