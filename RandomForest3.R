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


