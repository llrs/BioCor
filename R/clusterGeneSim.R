# clusterGeneSim ####
#' Compare two clusters of genes
#'
#' Looks for the similarity between genes of a group and then between each
#' group.
#'
#' Differs with clusterGeneSim that first each combination between genes is
#' calculated, and with this values then the comparison between the two
#' clusters is done. Thus applying combineScores twice, one at gene level and
#' another one at cluster level.
#' @param cluster1,cluster2 A vector with genes in \code{id}.
#' @inheritParams geneSim
#' @param method A vector with two  or one argument to be passed to combineScores the
#' first one is used to summarize the similarities of genes, the second one
#' for clusters.
#' @export
#' @author Lluís Revilla
#' @seealso \code{\link{clusterGeneSim}}, \code{\link{combineScores}} and \code{\link{conversions}}
#' @return \code{clusterGeneSim} returns a similarity score of the two clusters or
#' the similarity between the genes of the two clusters.
#' @examples
#' library("org.Hs.eg.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- select(org.Hs.eg.db, keys = entrezids, keytype = "ENTREZID",
#'                      columns = "PATH")
#' clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH")
#' clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", c("avg", "avg"))
#' clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", c("avg", "rcmax.avg"))
#' clus <- clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", "avg")
#' clus
#' combineScores(clus, "rcmax.avg")
clusterGeneSim <- function(cluster1, cluster2, genes, id, pathwayDB,
                        method = c("max", "rcmax.avg")) {
    if (length(method) > 2L | is.null(method) | any(is.na(method))) {
        stop("Please provide two  or one methods to combine scores.",
             "See Details")
    }

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

    vmpathSim <- Vectorize(mpathSim,
                           vectorize.args = c("pathways1", "pathways2"))

    out <- outer(pathways1, pathways2, vmpathSim, genes, id, pathwayDB, NULL)
    out <- apply(out, c(1,2), combineScores, method = method[1L])
    if (length(method) == 2) {
        combineScores(out, method[2L])
    } else {
        out
    }
}

#' @param clusters A list of clusters of genes to be found in \code{id}.
#' @rdname clusterGeneSim
#' @return \code{mclusterGeneSim} returns a matrix with the similarity scores for
#' each cluster comparison.
#' @export
#' @examples
#'
#' clusters <- list(cluster1 = c("18", "81", "10"),
#'                  cluster2 = c("100", "11", "1"),
#'                  cluster3 = c("18", "10", "83"))
#' mclusterGeneSim(clusters, genes.kegg, "ENTREZID", "PATH")
#' mclusterGeneSim(clusters, genes.kegg, "ENTREZID", "PATH", c("max", "avg"))
#' mclusterGeneSim(clusters, genes.kegg, "ENTREZID", "PATH", c("max", "BMA"))
mclusterGeneSim <- function(clusters, genes, id, pathwayDB,
                         method = c("max", "rcmax.avg")) {
    if (length(method) != 2) {
        stop("Please provide two methods to combine scores")
    }
    vc <- Vectorize(clusterGeneSim, vectorize.args = c("cluster1", "cluster2"))
    outer(clusters, clusters, vc, genes = genes, id = id,
          pathwayDB = pathwayDB, method = method)

}
