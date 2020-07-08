# mclusterGeneSim ####
#' Similarity score between clusters of genes based on genes similarity
#'
#' Looks for the similarity between genes of a group and then between each
#' group's genes.
#'
#' Differs with clusterSim that first each combination between genes is
#' calculated, and with this values then the comparison between the two
#' clusters is done. Thus applying combineScores twice, one at gene level and
#' another one at cluster level.
#' @inheritParams clusterSim
#' @inheritParams geneSim
#' @inheritParams pathSim
#' @param method A vector with two  or one argument to be passed to
#' combineScores the first one is used to summarize the similarities of genes,
#' the second one for clusters.
#' @export
#' @author Lluís Revilla
#' @seealso \code{\link{mclusterGeneSim}}, \code{\link{combineScores}} and
#' \code{\link{clusterSim}}
#' @return Returns a similarity score between the genes of the two clusters.
#' @examples
#' if (require("org.Hs.eg.db")) {
#'     # Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#'     # data of June 31st 2011)
#'     genes.kegg <- as.list(org.Hs.egPATH)
#'     clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg)
#'     clusterGeneSim(
#'         c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'         c("avg", "avg")
#'     )
#'     clusterGeneSim(
#'         c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'         c("avg", "rcmax.avg")
#'     )
#'     (clus <- clusterGeneSim(
#'         c("18", "81", "10"), c("100", "10", "1"),
#'         genes.kegg, "avg"
#'     ))
#'     combineScores(clus, "rcmax.avg")
#' } else {
#'     warning("You need org.Hs.eg.db package for this example")
#' }
clusterGeneSim <- function(cluster1, cluster2, info,
    method = c("max", "rcmax.avg"), ...) {
    if (length(unique(cluster1)) == 1L & length(unique(cluster2)) == 1L) {
        stop(
            "Introduce several genes in each cluster!\n",
            "If you want to calculate similarities ",
            "between two genes use geneSim"
        )
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
        stop(
            "Please provide two  or one methods to combine scores.",
            "See Details"
        )
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
        ... = ...
    )
    genes <- genes[names(pathways1.a), names(pathways2.a), drop = FALSE]

    if (length(method) == 2L) {
        combineScoresPar(as.matrix(genes), method = method[2L], ... = ...)
    } else {
        as.matrix(genes)
    }
}

#' @describeIn clusterGeneSim Calculates the gene similarities in a
#' GeneSetCollection and combine them using \code{\link{combineScoresPar}}
#' @export
setMethod(
    "clusterGeneSim",
    c(
        info = "GeneSetCollection", cluster1 = "character",
        cluster2 = "character"
    ),
    function(cluster1, cluster2, info, method, ...) {
        if (length(unique(cluster1)) == 1L & length(unique(cluster2)) == 1L) {
            stop(
                "Introduce several genes in each cluster!\n",
                "If you want to calculate similarities ",
                "between two genes use geneSim"
            )
        }
        if (!all(is.character(cluster1)) | !all(is.character(cluster2))) {
            stop("The input genes should be characters")
        }
        cluster1 <- unique(cluster1)
        cluster2 <- unique(cluster2)

        # Revert back to list
        list_info <- inverseList(GSEABase::geneIds(info))

        if (any(!cluster1 %in% names(list_info)) | any(!cluster2 %in% names(list_info))) {
            cluster1 <- cluster1[cluster1 %in% names(list_info)]
            cluster2 <- cluster2[cluster2 %in% names(list_info)]
            warning("Some genes are not in the GeneSetCollection provided.")
        }

        if (length(method) > 2L | is.null(method)) {
            stop(
                "Please provide two  or one method to combine scores.",
                "See Details"
            )
        }

        # FIXME: Should take advantage of GSC object (if any)
        # Extract all pathways for each gene
        pathways1.a <- lapply(cluster1, getElement, object = list_info)
        names(pathways1.a) <- cluster1
        pathways2.a <- lapply(cluster2, getElement, object = list_info)
        names(pathways2.a) <- cluster2

        # Remove duplicated and NA
        pathways1 <- unique(unlist(pathways1.a, use.names = FALSE))
        pathways2 <- unique(unlist(pathways2.a, use.names = FALSE))
        pathways1 <- pathways1[!is.na(pathways1)]
        pathways2 <- pathways2[!is.na(pathways2)]

        if (is.null(pathways1) || is.null(pathways2)) {
            return(NA)
        }
        pathways <- unique(c(pathways1, pathways2))

        simPaths <- mpathSim(pathways, list_info, method = NULL, ...)
        genes <- combineScoresPar(simPaths, method[1L],
            c(pathways1.a, pathways2.a),
            ... = ...
        )
        genes <- genes[names(pathways1.a), names(pathways2.a), drop = FALSE]

        if (length(method) == 2L) {
            combineScoresPar(as.matrix(genes), method = method[2L], ... = ...)
        } else {
            as.matrix(genes)
        }
    }
)
