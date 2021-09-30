#' plotdensity_flowset 
#'
#' This function creates a density plot for each channel in flow cytometery data 
#' The input in a data object containg the expression values from fsc files 
#' This matrix is obtain from flow cytometry and read in using 'flowset'
#' This funciton will read the intensity of expression for each channel for each cell
#' The FSC-A values are taken (average forward side scatter)
#' The function than calles another function 'rename_markers' to change the antibody names
#' 
#' 
#'
#' @param infile Path to the input file
#' @return A matrix of the infile
#' @export