# Annotate clusters functions:

# see_features                (Shuming/Rhalena *** new)
# CAM                         (Shuming)
# RFM_train                   (Rhalena)
# RFM_predict                 (Rhalena)
# seurat_transfer             (Rhalena)
# cluster_annotate            (Rhalena)



##############################################################################################

# see_features
# takes in a seurat object with clustering and UMAP run in explore_param or otherwise
# creates featuremaps for each antibody
# creates a UMAP with cluster numbers labelled
# creates a heatmap split by cluster - takes an argument for which cluster resolution is desired
# takes in a features list













##############################################################################################

# CAM
# input correlation matrix
# input seurat object
#Correlation assignment method, 
# predicts cell types based on each cells correlation to the matrix types
# creates plots and tables of the prediction outputs. 
# takes arguement for "unknown" threshold (default 0.4)  and "double-lable" thresholod(0.05)






##############################################################################################
# RFM_train
# Input annotated FC dataset 
# Random Forest Model internally optimizing parameters and saving the best model. 
# requires randomForst, caret, data.table




# prepare data
seu.r<- readRDS(pathway_to_data)
AB <- c("CD24","CD56","CD29","CD15","CD184","CD133","CD71","CD44","GLAST","AQP4","HepaCAM", "CD140a","O4")


df <- transpose(as.data.frame(GetAssayData(seu.r,slot = 'scale.data')))
dim(df)



mtry.range <- c(1:10)

set.seed(48)

tuneGrid <- expand.grid(.mtry = mtry.range)
rf_mtry <- train(lables~., 
                 data=train,
                 method = "rf",
                 metric = "Accuracy",
                 tuneGrid = tuneGrid,
                 trControl = trControl,
                 importance = TRUE,
                 nodesize = 15,
                 ntree = 300)
rf_mtry$bestTune$mtry
max(rf_mtry$results$Accuracy)

best_mtry <- rf_mtry$bestTune$mtry 
best_mtry

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






##############################################################################################
# RFM_predict
# take the saved trained RFM
# takes seurat data object
# predicts cell types
# creates table and umap
# outputs prediction vector





##############################################################################################
# seurat_transfer
# takes in a reference seurat object
# follows seurat workflow, finds anchors, predicts labels
# default of no threshold, a threshold can be applied
# outputs prediction vector









##############################################################################################

# cluster_annotate
# input vectors of cluster annotation, at least 1 is required, prediction vector from visualization
# type manually, be careful to match syntax to the other predictions
# input seurat data object
# finds most common label per cluster and applies that label


