#' BioCor: A package to calculate functional similarities
#'
#' Calculates a functional similarity measure between gene identifiers
#' based on the pathways described on KEGG and REACTOME.
#'
#' @section Important functions:
#' \describe{
#' \item{\code{\link{pathSim}}}{Calculates the similarity between two pathways}
#' \item{\code{\link{genesSim}}}{Calculates the similarity (based on pathSim)
#' between two genes}
#' \item{\code{\link{clusterSim}}}{Calculates the similarity between two
#' clusters of genes}
#' \item{\code{\link{clustersSim}}}{Calculates the similarity between two
#' clusters of genes}
#' \item{\code{\link{bioCor}}}{Allows to calculate the similarity for a large
#' number of genes
#'  using parallel back-end}
#' \item{\code{\link{similarities}}}{Allows to combine the value of matrices of
#' similarities}
#' \item{\code{\link{conversions}}}{Two functions to convert similarity
#' measures}
#' \item{\code{\link{weigthed}}}{Functions provided to combine similarities}
#' }
#' @name BioCor-package
#' @aliases BioCor
#' @docType package
NULL
