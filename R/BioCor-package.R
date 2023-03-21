#' BioCor: A package to calculate functional similarities
#'
#' Calculates a functional similarity measure between gene identifiers
#' based on the pathways described on KEGG and REACTOME.
#'
#' @section Important functions:
#' \describe{
#' \item{[pathSim()]}{Calculates the similarity between two pathways}
#' \item{[geneSim()]}{Calculates the similarity (based on pathSim)
#' between two genes}
#' \item{[clusterSim()]}{Calculates the similarity between two
#' clusters of genes by joining pathways of each gene.}
#' \item{[clusterGeneSim()]}{Calculates the similarity between two
#' clusters of genes by comparing the similarity between the genes of a cluster
#' }
#' \item{[similarities()]}{Allows to combine the value of matrices of
#' similarities}
#' \item{[conversions()]}{Two functions to convert similarity
#' measures}
#' \item{[weighted()]}{Functions provided to combine similarities}
#' }
#' @name BioCor-package
#' @aliases BioCor
#' @docType package
NULL
