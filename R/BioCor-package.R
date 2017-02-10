#' BioCor: A package to calculate functional similarities
#'
#' Calculates a functional similarity measure between gene identifiers
#' based on the pathways described on KEGG and REACTOME.
#'
#' @section Important functions:
#' \describe{
#' \item{pathSim}{Calculates the similarity between two pathways}
#' \item{genesSim}{Calculates the similarity (based on pathSim)
#' between two genes}
#' \item{bioCor}{Allows to calculate the similarity for a large number of genes
#'  using parallel back-end}
#' \item{similarities}{Allows to combine the value of matrices of similarities}
#' \item{conversions}{Two functions to convert similarity measures}
#' \item{weigthed}{Functions provided to combine similarities}
#' }
#' @name BioCor-package
#' @aliases BioCor
#' @docType package
NULL
