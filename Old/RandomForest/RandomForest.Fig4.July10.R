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

# search the best maxnodes
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

fit_rf <- train(lables~.,
                train,
                method = "rf",
                metric = "Accuracy",
                tuneGrid = tuneGrid,
                trControl = trControl,
                importance = TRUE,
                nodesize = 25,
                ntree = 1000,
                maxnodes = 30)
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
c2 <- confusionMatrix(p2, test$lables)
c2.table <- as.data.frame(c2$table)

# try to plot the results

library(ggplot2)
ggplot(data =  c2.table, mapping = aes(x = Prediction, y = Reference, fill= Freq)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f",Freq)), vjust = 1) +
  scale_fill_gradient(low = "white", high = "steelblue") + theme(axis.text.x = element_text(angle = 90))


# save the model for later

saveRDS(rf, "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/rf_trained_march15.Rds")
write.csv(c2.table, "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/confusionm_values_march15.csv")