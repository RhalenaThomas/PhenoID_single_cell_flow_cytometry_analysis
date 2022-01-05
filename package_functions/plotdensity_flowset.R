#' plotdensity_flowset 
#'
#' This function creates a density plot for each channel in flow cytometery data 
#' The input in a data object containing the expression values from fsc files 
#' This matrix is obtain from flow cytometry and read in using the function 'flowset' from the package "FlowCore"
#' The input data is flowset object made with the make_flowset function (PhenoID package)
#' 
#'
#' @param infile flowset data object
#' @return Density plots
#' @export
#' 
#' 
plotdensity_flowset <- function(flowset){ ggplot(melt(lapply(as.list(flowset@frames),
                                                             function(x){x=as.data.frame(x@exprs)})), 
                                                 aes(x=value,y=L1,fill=L1)) + geom_density_ridges(alpha=.4,verbose=FALSE) 
  +facet_wrap(~variable)+theme_light()} #defines a function for visualizing flowset with densityplots