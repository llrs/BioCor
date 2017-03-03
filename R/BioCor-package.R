#' BioCor: A package to calculate functional similarities
#'
#' Calculates a functional similarity measure between gene identifiers
#' based on the pathways described on KEGG and REACTOME.
#'
#' @section Important functions:
#' \describe{
#' \item{\code{\link{pathSim}}}{Calculates the similarity between two pathways}
#' \item{\code{\link{geneSim}}}{Calculates the similarity (based on pathSim)
#' between two genes}
#' \item{\code{\link{clusterSim}}}{Calculates the similarity between two
#' clusters of genes by joining pathways of each gene.}
#' \item{\code{\link{clusterGeneSim}}}{Calculates the similarity between two
#' clusters of genes by comparing the similarity between the genes of a cluster}
#' \item{\code{\link{similarities}}}{Allows to combine the value of matrices of
#' similarities}
#' \item{\code{\link{conversions}}}{Two functions to convert similarity
#' measures}
#' \item{\code{\link{weighted}}}{Functions provided to combine similarities}
#' }
#' @name BioCor-package
#' @aliases BioCor
#' @docType package
NULL
