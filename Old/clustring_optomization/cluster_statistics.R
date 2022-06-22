# PhenoID
# identify best clustering method and conditions for multiple antibody panel in flow cytometry to identify cell types

#################################################################################################################################
########################                             clustering methods                                ##########################
# hierarchical clustering   
# k and x-means
# density based spatial clustering   - Dbscan in python (sklearn) R package: install.packages("dbscan")
# flowsom 
# network clustering:
#   phenograph (phenograph uses louvain but add in the jaccard index)
#   louvain
# get statistics on each method
# Start with flowsom


################################################################################################################################
#######################                   Statistics to compare clustering methods                    ##########################

# intrinsic methods we can use before we know the ground truth, which is our situation
# these methods will evaluate how well the clusters separated

#silhouette coefficient 
# calinski-harabasz index
# Davies-Bouldin Index


# methods that can be applied after expert labeling
# all of the following require ground truth/ labeled clusters (hung)
# hungarian index 
# RAND index   (RI)
# Adjusted RAND index (ARI)
# Adjusted mutal information (AMI)
# normaized mutual information (NMI)
# BCubed

###############################################################################################################################
# Install all required packages
install.packages("mclust") # clustering and hungarian index
install.packages("aricode")




###############################################################################################################################
############################################ Read in the starting data object #################################################
# starting with aligned normalized data from the pre-processing output

# for testing I'm using the csv file that is a subset of 5000 cells from the 9 samples of different MBOs

df <- read.csv("/home/rhalena/Documents/Documents/MyPapers/FACS/Data/Large9MBO/phenoID/June14_2021_testclustering/subset_clean_data.csv")
dim(df)

# I'm removing the cell indexes and the size measurements
df2 <- df[,4:16]



#################################################################################################################################
########################                             clustering methods                                ##########################
library(mclust)
# Hierarchical clustering - average method
# The following code creates HC 
ob_dist <- dist(df)  # create the distance matrix between rows in a matrix
library(stats)
hc <- hclust(ob_dist, method="average")
plot(hc)

#To cut for a give k for instance k=4 
cut1 <- cutree(hc, k=4)
#To cut for a give value of distance, 0.4 
cutree(hc, h=.4)


# find the number of clusters within the data
mclustBIC




# density based spatial clustering 






###############################################################################################################################
#########################                               Clustering statistics                    ###############################
###############################################################################################################################

######################################## Instrinsic statistics to use before we establish ground true #########################







########################################## Statistics after ground truth is established ########################################


# hungarian index
library(mclust) 
hung <-function(trueclu,obsclus)  1-classError(trueclu, obsclus)$errorRate

cu_or<-cutree(hclust(REDIST,method = "average"), k =6)
avcorrec2(cu_or,class0)

# # RAND index   (RI)
# Adjusted RAND index (ARI)
# Adjusted mutal information (AMI)
# normaized mutual information (NMI)

library(aricode)

ami <-AMI(label,true_class)
nmi <- NMI(label,true_class)
ri <- RI(label,true_class)
ari <- ARI(label,true_class)















