#' BioCor: A package to calculate functional similarities
#'
#' Calculates a functional similarity measure between gene identifiers
#' based on the pathways described on KEGG and REACTOME.
#'
#' @section Important functions:
#'  - **[pathSim()]**: Calculates the similarity between two pathways.
#'  - **[geneSim()]**: Calculates the similarity (based on pathSim)
#'  between two genes.
#'  - **[clusterSim()]**: Calculates the similarity between two
#' clusters of genes by joining pathways of each gene.
#'  - **[clusterGeneSim()]**: Calculates the similarity between two
#' clusters of genes by comparing the similarity between the genes of a cluster.
#'  - **[similarities()]**: Allows to combine the value of matrices of
#' similarities.
#'  - **[conversions()]**: Two functions to convert similarity
#' measures.
#'  - **[weighted()]**: Functions provided to combine similarities.
'_PACKAGE'
