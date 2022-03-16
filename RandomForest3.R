library(randomForest)
library(caret)
library(Seurat)
library(data.table)



# add some data
seu <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/Louvkn60DifferentCellLabels09032022.Rds")

df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))
dim(df)


# I need to know the names of the columns now
# from the input df
col.names <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")

colnames(df) <- col.names
print(col.names)
# add in the cell lables

annotations <- seu$labels4

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

# like in the tutorial the highest max node has the highest value
# higher values might give better accuracy
# the accuracy is already very good.

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
# best number of trees is 600

# now fit the model with the best conditions
#  mtry = 5
# maxnodes = 29
# node size 15
# ntree = 600

fit_rf <- train(lables~.,
                train,
                method = "rf",
                metric = "Accuracy",
                tuneGrid = tuneGrid,
                trControl = trControl,
                importance = TRUE,
                nodesize = 15,
                ntree = 600,
                maxnodes = 29)
# not sure if mtry is being input as a variable or not

prediction.train <-predict(fit_rf, train)
prediction.test <-predict(fit_rf, test)

print("predict training data")
confusionMatrix(prediction.train, train$lables)

print("predict test data")
confusionMatrix(prediction.test, test$lables)

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
                   nodesize = 15,
                   ntree = 600,
                   maxnodes = 29)

# now the predictions
p1 <- predict(rf, train)
confusionMatrix(p1, train$lables)

p2 <- predict(rf, test)
confusionMatrix(p2, test$lables)

