# preprocessing function
# open this file and run the function 
# will go into a library later

# input arguments
# takes in input pathway to fsc files folder can only contain fsc
# output folder for results
# if wanted : rename files for 'sample' names
# rename channels (column names)
# subset cells/down sample - input a value if left empty all cells will be taken from each file
# select desired processing - default all
# might need to take arguments for adjusting alignment 
# option to create fsc files - I'm not sure why can these be opened in flowjo?


# outputs
# csv files for: concationed and renamed only, aligned (transformed), aligned and retrotransformed
# png of density plots - different transformations




## set up

# setup the enviroment
# multiple packages may need to be installed

require("flowCore") #Used for reading the data
require("ggplot2")
require("ggridges") #visualization
require("stringr") #set of functions to manipulate strings type in order to have nice titles
require("rlist") #set of functions that allows to easily manipulate lists
require("reshape2") #visualization
require("flowStats") #Alignment functions
require("scales") #scale colour intensity for visualization
require("dplyr")

#libraries
library("flowCore")

# installations
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("flowCore")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("flowStats")


### functions
# functions Alex wrote
# will need to make these into a separate function with documentation

# function plotdensity_flow set

plotdensity_flowset <- function(flowset){ ggplot(melt(lapply(as.list(flowset@frames),function(x){x=as.data.frame(x@exprs)})), aes(x=value,y=L1,fill=L1)) + geom_density_ridges(alpha=.4,verbose=FALSE) +facet_wrap(~variable)+theme_light()} 

#defines a function for visualizing flowset with densityplots

# function renmame markers
# fcs files will have the marker names input during acquistion using flowjo

rename_markers<-function(flowset){#Defines a function to use marker names 
  copy_flowset=flowset[seq(along=flowset)]
  for (i in 1:length(copy_flowset)){
    marker.names=copy_flowset[[i]]@parameters@data$desc
    marker.names=lapply(marker.names,function(x){str_replace_all(x,"-","_")})
    colnames(copy_flowset[[i]]@exprs)<-unlist(lapply(marker.names, function(x){sapply(str_split(x,"_"),head,1)})) 
  }
  return(copy_flowset)
}







