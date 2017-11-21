# clusterGeneSim ####
#' Similarity score between clusters of genes based on genes similarity
#'
#' Looks for the similarity between genes of a group and then between each
#' group.
#'
#' Differs with clusterGeneSim that first each combination between genes is
#' calculated, and with this values then the comparison between the two
#' clusters is done. Thus applying combineScores twice, one at gene level and
#' another one at cluster level.
#' @inheritParams clusterSim
#' @inheritParams geneSim
#' @inheritParams pathSim
# #' @import AnnotationDbi
# #' @import org.Hs.eg.db
# #' @import reactome.db
#' @param method A vector with two  or one argument to be passed to
#' combineScores the first one is used to summarize the similarities of genes,
#' the second one for clusters.
#' @export
#' @author Lluís Revilla
#' @seealso \code{\link{clusterGeneSim}}, \code{\link{combineScores}} and
#' \code{\link{conversions}}
#' @return \code{clusterGeneSim} returns a similarity score of the two clusters
#' or the similarity between the genes of the two clusters.
#' @examples
#' if (require("org.Hs.eg.db")) {
#'     #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#'     # data of June 31st 2011)
#'     genes.kegg <- as.list(org.Hs.egPATH)
#'     clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg)
#'     clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'                    c("avg", "avg"))
#'     clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'                    c("avg", "rcmax.avg"))
#'     (clus <- clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'                             "avg"))
#'     combineScores(clus, "rcmax.avg")
#' } else {
#'     warning('You need org.Hs.eg.db package for this example')
#' }
clusterGeneSim <- function(cluster1, cluster2, info,
                           method = c("max", "rcmax.avg"), ...) {
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

    if (any(!cluster1 %in% names(info)) | any(!cluster2 %in% names(info))) {
        warning("Some genes are not in the list provided.")
    }

    if (length(method) > 2L | is.null(method)) {
        stop("Please provide two  or one methods to combine scores.",
             "See Details")
    }

    # Extract all pathways for each gene
    pathways1.a <- lapply(cluster1, getElement, object = info)
    names(pathways1.a) <- cluster1
    pathways2.a <- lapply(cluster2, getElement, object = info)
    names(pathways2.a) <- cluster2

    # Remove duplicated and NA
    pathways1 <- unique(unlist(pathways1.a, use.names = FALSE))
    pathways2 <- unique(unlist(pathways2.a, use.names = FALSE))
    pathways1 <- pathways1[!is.na(pathways1)]
    pathways2 <- pathways2[!is.na(pathways2)]

    pathways <- unique(c(pathways1, pathways2))
    if (is.null(pathways1) || is.null(pathways2)) {
        return(NA)
    }
    simPaths <- mpathSim(pathways, info, method = NULL, ...)
    genes <- combineScoresPar(simPaths, method[1L],
                              c(pathways1.a, pathways2.a),
                              ... = ...)
    genes <- genes[names(pathways1.a), names(pathways2.a), drop = FALSE]

    if (length(method) == 2L) {
        combineScoresPar(as.matrix(genes), method = method[2L], ... = ...)
    } else {
        as.matrix(genes)
    }

}

vclusterGeneSim <- Vectorize(clusterGeneSim,
                             vectorize.args = c("cluster1", "cluster2"))

#' @param clusters A list of clusters of genes to be found in \code{id}.
#' @rdname clusterGeneSim
#' @return \code{mclusterGeneSim} returns a matrix with the similarity scores
#' for each cluster comparison.
#' @export
#' @examples
#'
#' clusters <- list(cluster1 = c("18", "81", "10"),
#'                  cluster2 = c("100", "594", "836"),
#'                  cluster3 = c("18", "10", "83"))
#' mclusterGeneSim(clusters, genes.kegg)
#' mclusterGeneSim(clusters, genes.kegg, c("max", "avg"))
#' mclusterGeneSim(clusters, genes.kegg, c("max", "BMA"))
mclusterGeneSim <- function(clusters, info, method = c("max", "rcmax.avg"),
                            ...) {

    if (!is.list(clusters)) {
        stop("Please use a list to introduce the clusters.")
    }
    if (length(clusters) == 1) {
        stop("Introduce several clusters!\n",
             "If you want to calculate the similarity ",
             "between genes use mgeneSim")
    }
    if (!all(sapply(clusters, is.character))) {
        stop("The input genes should be characters")
    }
    # Remove duplicate genes in each cluster
    clusters <- lapply(clusters, unique)

    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    if (any(!unlist(clusters, use.names = FALSE) %in% names(info))) {
        warning("Some genes are not in the list provided.")
    }
    if (length(method) != 2) {
        stop("Please provide two methods to combine scores")
    }

    # Extract all pathways for each gene
    pathways <-lapply(unlist(clusters, use.names = FALSE), function(x) {
        info[[x]]
    })
    pathwaysl <- unique(unlist(pathways, use.names = FALSE))
    pathwaysl <- pathwaysl[!is.na(pathwaysl)]

    # Calculate similarities between pathways
    pathSims <- mpathSim(pathwaysl, info, NULL)

    # Calculate similarities between genes
    names(pathways) <- unlist(clusters, use.names = FALSE) # give the name of the genes
    genesSims <- combineScoresPar(pathSims, method[1L], pathways, ... = ...)

    # Calculate similarities between clusters of genes
    as.matrix(combineScoresPar(genesSims, method[2L], clusters, ... = ...))
}
