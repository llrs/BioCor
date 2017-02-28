# clusterSim ####
#' Compare two clusters of genes
#'
#' Looks for the similarity between genes in groups
#'
#' Once the pathways for each cluster are found they are combined using
#' combineScores.
#' @param cluster1,cluster2 A vector with genes in \code{id}.
#' @inheritParams geneSim
#' @export
#' @author Lluís Revilla
#' @seealso For a different approach see \code{\link{clustersSim}},
#' \code{\link{combineScores}} and \code{\link{conversions}}
#' @return \code{clusterSim} returns a similarity score of the two clusters
#' @examples
#' library("org.Hs.eg.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- select(org.Hs.eg.db, keys = entrezids, keytype = "ENTREZID",
#'                      columns = "PATH")
#' clusterSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH")
#' clusterSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", NULL)
#' clusterSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", "avg")
clusterSim <- function(cluster1, cluster2, genes, id, pathwayDB, method = "max"){

    # Convert data.frame into environment to speed the look up
    genes2pathways <- split(genes[ , pathwayDB], genes[, id])
    genes2pathways <- list2env(genes2pathways)

    # Extract all pathways for each gene
    pathways1 <- sapply(cluster1, function(x) {
            genes2pathways[[x]]
        }, simplify = FALSE)
    pathways2 <- sapply(cluster2, function(x) {
        genes2pathways[[x]]
    }, simplify = FALSE)

    pathways1 <- unlist(pathways1, use.names = FALSE)
    pathways2 <- unlist(pathways2, use.names = FALSE)

    mpathSim(pathways1, pathways2, genes, id, pathwayDB, method)
}

#' @param clusters A list of clusters of genes to be found in \code{id}.
#' @rdname clusterSim
#' @return \code{mclusterSim} returns a matrix with the similarity scores for
#' each cluster comparison.
#' @export
#' @examples
#'
#' clusters <- list(cluster1 = c("18", "81", "10"),
#'                  cluster2 = c("100", "10", "1"),
#'                  cluster3 = c("18", "10", "83"))
#' mclusterSim(clusters, genes.kegg, "ENTREZID", "PATH")
#' mclusterSim(clusters, genes.kegg, "ENTREZID", "PATH", "avg")
mclusterSim <- function(clusters, genes, id, pathwayDB, method = "max") {
    vc <- Vectorize(clusterSim, vectorize.args = c("cluster1", "cluster2"))
    outer(clusters, clusters, vc, genes = genes, id = id,
          pathwayDB = pathwayDB, method = method)

}
