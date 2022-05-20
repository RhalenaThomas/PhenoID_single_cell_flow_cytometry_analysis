library(randomForest)
library(caret)
library(Seurat)
library(data.table)


seu <- readRDS("/Users/shumingli/Documents/Louvkn60DifferentCellLabels220220318.Rds")

df <- transpose(as.data.frame(GetAssayData(seu,slot = 'scale.data')))

col.names <- c("AQP4", "CD24", "CD44","CD184","CD15","HepaCAM","CD29","CD56", "O4","CD140a","CD133","GLAST","CD71")
colnames(df) <- col.names

annotations <- seu$labels6

df.2 <- cbind(df, label = annotations)
print(table(df.2$label))

df.2

ind <- sample(2, nrow(df.2), replace = TRUE, prob = c(0.5, 0.5)) # prop is the proportions

train <- df.2[ind==1,]
test <- df.2[ind==1,]


model <- randomForest(x=df[ind==1,],y=train$label)

model

typeof(data)

