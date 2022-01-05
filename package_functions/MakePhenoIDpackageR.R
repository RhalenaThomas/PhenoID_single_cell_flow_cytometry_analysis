# set up to make a package
# we need two packages devtools and roxygen2

install.packages("devtools")
library("devtools")
install.packages("roxygen2")
library("roxygen2")

# now create the frame work of the package
# this will create a folder name phenoID that will be the name of the package
devtools::create("phenoID")

# there will be 4 files inside the folder
# 1. Description: the meta-data about the package goes here
# 2. phenoID.Rproj: this is for R studio
# 3. NAMESPACE: indicates what needs to be exposed to users for the R package. Do not edit
# 4. R: this is where all the R code goes for the package


# start by filing out the description for the package

## Add functions

# let make one .R file for each function with the name of the file being the name of the function
