
library(randomForest)
library(caret)
library(Seurat)
library(data.table)

# add some data
seu <- readRDS("/Users/rhalenathomas/Documents/Data/FlowCytometry/PhenoID/Analysis/9MBO/prepro_outsjan20-9000cells/Figure3/Louvkn60DifferentCellLabels09032022.Rds")



# follow a tutorial
library(datasets)

data <- iris
str(data)

class(data)

# the response variable is the label or "true"
data$Species <- as.factor(data$Species)

table(data$Species)



# set the seed
set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3)) # prop is the proportions
# split the data into train and test
train <- data[ind==1,]
test <- data[ind==2,]


# now tran the classifier and run the test
rf <- randomForest(Species~., data=train, proximity=TRUE) 
print(rf) # gives the random forest results
#Call:
 # randomForest(formula = Species ~ ., data = train)
#Type of random forest: classification
#Number of trees: 500
#No. of variables tried at each split: 2
#OOB estimate of  error rate: 2.83%

# predict and confusion marix on training data
p1 <- predict(rf, train)
confusionMatrix(p1, train$ Species)

# predict and confusion matrix on test data
p2 <- predict(rf, test)
confusionMatrix(p2, test$ Species)

plot(rf)


# https://www.r-bloggers.com/2021/04/random-forest-in-r/#:~:text=Random%20Forest%20in%20R%2C%20Random,to%20identify%20the%20important%20attributes.
# there are tuning parameters

# partial dependence plot
partialPlot(rf, train, Petal.Width, "setosa")
# shows petal width less than 1.5 increase changes of a calls of setosa

MDSplot(rf, train$Species)



### try with the real data

# make a df

df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))
dim(df)


# I need to know the names of the columns now
# from the input df
col.names <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")

colnames(df) <- col.names

# add in the cell lables

annotations <- seu$labels4

df.l <- cbind(df, lables = annotations)

table(df.l$lables)

# split the data
set.seed(222)
ind <- sample(2, nrow(data), replace = TRUE, prob = c(0.7, 0.3)) # prop is the proportions
# split the data into train and test
train <- df.l[ind==1,]
test <- df.l[ind==2,]

# now tran the classifier and run the test
rf <- randomForest(labels~., data=train) 
print(rf) # gives the random forest results

class(train)

