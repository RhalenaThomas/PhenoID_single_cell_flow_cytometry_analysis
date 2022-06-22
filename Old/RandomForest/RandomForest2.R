
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
train.old <- df.l[ind==1,]
test.old <- df.l[ind==2,]

# now train the classifier and run the test
# 
rf <- randomForest(lables~., data=test, proximity=TRUE, ntrees = 50)  
# memory error with the full input - not sure about more trees n= 5 is very low

# try the training with tuning parmaters
t <- tuneRF(train[,-5], train[,5],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)


# save the model
saveRDS(rf, "/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/RandomForest/Louvkn60DifferentCellLabels09032022.Rds")


print(rf)
# predictions 
p1 <- predict(rf, test)
confusionMatrix(p1, test$lables)
table(p1)

p2 <- predict(rf, train)
confusionMatrix(p2, train$lables)


list.p2 <- confusionMatrix(p2, train$lables)

r.table <- as.data.frame(list.p2[["table"]])

r.overall <- as.data.frame(list.p2[["overall"]])

# try to improve parameters

t <- tuneRF(train[,-5], train[,5],
            stepFactor = 0.5,
            plot = TRUE,
            ntreeTry = 150,
            trace = TRUE,
            improve = 0.05)

# further look at the results
#MDSplot(rf, test$lables) # memore limit reach - 

hist(treesize(rf),
     main = "No. of Nodes for the Trees",
     col = "green")
#Variable Importance
varImpPlot(rf,
           sort = T,
           n.var = 13,
           main = "Variable Importance")
importance(rf)
#MeanDecreaseGini plot 

# In order of importance CD44, CD56, CD184, CD33 

# Partial plot
# put in the importance feature, the class to distinguish
partialPlot(rf, train, CD44, "Glia")

# which variables were actually used?

varUsed(rf, by.tree=FALSE, count=TRUE)

