# Create cell type clusters functions

# explore_param        (Shuming)
#Intrinsic stats       (Shuming)
# # clust_stability    (Shuming)


####################################################################################################

# explore_param
# reads in csv files with flow cytometry experiments or seurat data object depending on arguments
# runs FlowSom, Phenograph, Louvain (Seurat) 
# input arguments - cluster method, for flowsom k, for phenograph and louvain k paramater, for louvain resolution
# select outputs, generates - UMAPS, heatmaps, clustree  : save to folder or put in global environ
# creates list object to run stats

####################################################################################################


#Intrinsic stats
# takes in a dataframe or list of statistics
# plots intrinsic statistics from pararmater explorations




####################################################################################################

# clust_stability
# select cluster method and one pararmeter to vary (resolutions 0.1,0.3,0.5,0.8) 
# runs n iterations randomizing starting point - default 100
# calculates the RandIndex (RI) between each iteration of clustering for each resolution
# calculates the mean and standard deviation of the number of clusters and the RI for each resoloution
# outputs a table and plot of the results





