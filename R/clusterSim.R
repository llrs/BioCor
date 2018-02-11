# clusterSim ####
#' Similarity score between clusters of genes based on pathways similarity
#'
#' Looks for the similarity between genes in groups
#'
#' Once the pathways for each cluster are found they are combined using
#' \code{\link{combineScores}}.
#' @param cluster1,cluster2 A vector with genes.
#' @inheritParams geneSim
#' @inheritParams combineScores
#' @inheritParams pathSim
#' @export
#' @author Lluís Revilla
#' @seealso For a different approach see \code{\link{clusterGeneSim}},
#' \code{\link{combineScores}} and \code{\link{conversions}}
#' @return \code{clusterSim} returns a similarity score of the two clusters
#' @examples
#' if (require("org.Hs.eg.db")) {
#'     #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#'     # data of June 31st 2011)
#'     genes.kegg <- as.list(org.Hs.egPATH)
#'     clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg)
#'     clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg, NULL)
#'     clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg, "avg")
#' } else {
#'     warning('You need org.Hs.eg.db package for this example')
#' }
clusterSim <- function(cluster1, cluster2, info, method = "max", ...){

    if (length(unique(cluster1)) == 1L & length(unique(cluster2)) == 1L) {
        stop("Introduce several genes in each cluster!\n",
             "If you want to calculate similarities ",
             "between two genes use geneSim")
    }
    if (!all(is.character(cluster1)) | !all(is.character(cluster2))) {
        stop("The input genes should be characters")
    }
    cluster1 <- unique(cluster1)
    cluster2 <- unique(cluster2)

    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    clust1_logic <- !cluster1 %in% names(info)
    clust2_logic <- !cluster2 %in% names(info)

    if (all(clust1_logic) & all(clust2_logic)) {
        warning("At least one gene should be in the list provided")
        return(NA)
    } else if (any(clust1_logic) | any(clust2_logic)) {
        warning("Some genes are not in the list provided.")
    }

    # Extract all pathways for each gene
    pathways1 <- lapply(cluster1, getElement, object = info)
    pathways2 <- lapply(cluster2, getElement, object = info)

    # Remove duplicated and NA
    pathways1 <- unique(unlist(pathways1, use.names = FALSE))
    pathways2 <- unique(unlist(pathways2, use.names = FALSE))
    pathways1 <- pathways1[!is.na(pathways1)]
    pathways2 <- pathways2[!is.na(pathways2)]

    if (is.null(pathways1) & is.null(pathways2)) {
        return(NA)
    }

    pathways <- unique(c(pathways1, pathways2))

    sim_all <- mpathSim(pathways, info, NULL)
    sim <- sim_all[pathways1, pathways2]
    if (!is.null(method)) {
        combineScoresPar(sim, method, ...)
    } else {
        sim
    }
}


#' @describeIn clusterSim Calculates all the similarities of the
#' GeneSetCollection and combine them using \code{\link{combineScoresPar}}
#' @export
setMethod("clusterSim",
          c(info = "GeneSetCollection", cluster1 = "character",
            cluster2 = "character"),
          function(cluster1, cluster2, info, method, ...) {
              if (length(unique(cluster1)) == 1L & length(unique(cluster2)) == 1L) {
                  stop("Introduce several genes in each cluster!\n",
                       "If you want to calculate similarities ",
                       "between two genes use geneSim")
              }

              # Check they are unique
              cluster1 <- unique(cluster1)
              cluster2 <- unique(cluster2)

              # Extract the ids
              origGenes <- geneIds(info)
              # Check that the genes are in the GeneSetCollection
              genes <- unique(unlist(origGenes, use.names = FALSE))
              if (all(!cluster1 %in% genes)) {
                  warning("At least one gene should be in the list provided")
                  return(NA)
              }

              if (all(!cluster2 %in% genes)) {
                  warning("At least one gene should be in the list provided")
                  return(NA)
              }

              # Simplify the GeneSetCollection
              keep <- sapply(origGenes, function(x) {
                  any(c(cluster1, cluster2) %in% x)
              })
              gscGenes <- info[names(keep[keep])]

              # Search for the paths of each gene
              clusters <- list(cluster1 = cluster1, cluster2 = cluster2)
              ids <- geneIds(gscGenes)
              paths <- sapply(clusters, function(x){
                  keepPaths <- sapply(ids, function(y) {
                      any(x %in% y)
                  })
                  names(keepPaths[keepPaths])
              })

              # Calculate the pathSim of all the implied pathways
              pathsSim <- mpathSim(info = gscGenes, method = NULL)
              # Summarize the information
              out <- combineScoresPar(pathsSim, method, subSets = paths)
              out["cluster1", "cluster2"]
          }
)
