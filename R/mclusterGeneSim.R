#' Similarity score between clusters of genes based on genes similarity
#'
#' Looks for the similarity between genes of a group and then between each
#' group's genes.
#'
#' @inheritParams clusterSim
#' @inheritParams geneSim
#' @inheritParams pathSim
#' @inheritParams clusterGeneSim
#' @param method A vector with two  or one argument to be passed to
#' combineScores the first one is used to summarize the similarities of genes,
#' the second one for clusters.
#' @export
#' @author Lluís Revilla
#' @seealso \code{\link{clusterGeneSim}}, \code{\link{clusterSim}} and
#' \code{\link{combineScores}}
#' @param clusters A list of clusters of genes to be found in \code{id}.
#' @return Returns a matrix with the similarity scores for each cluster
#' comparison.
#' @export
#' @examples
#' if (require("org.Hs.eg.db")) {
#'     genes.kegg <- as.list(org.Hs.egPATH)
#'     clusters <- list(
#'         cluster1 = c("18", "81", "10"),
#'         cluster2 = c("100", "594", "836"),
#'         cluster3 = c("18", "10", "83")
#'     )
#'     mclusterGeneSim(clusters, genes.kegg)
#'     mclusterGeneSim(clusters, genes.kegg, c("max", "avg"))
#'     mclusterGeneSim(clusters, genes.kegg, c("max", "BMA"))
#' } else {
#'     warning("You need org.Hs.eg.db package for this example")
#' }
mclusterGeneSim <- function(clusters, info, method = c("max", "rcmax.avg"),
    ...) {
    if (!is.list(clusters)) {
        stop("Please use a list to introduce the clusters.")
    }
    if (length(clusters) == 1) {
        stop(
            "Introduce several clusters!\n",
            "If you want to calculate the similarity ",
            "between genes use mgeneSim"
        )
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
    pathways <- lapply(unlist(clusters, use.names = FALSE), function(x) {
        info[[x]]
    })
    pathwaysl <- unique(unlist(pathways, use.names = FALSE))
    pathwaysl <- pathwaysl[!is.na(pathwaysl)]

    # Calculate similarities between pathways
    pathSims <- mpathSim(pathwaysl, info, NULL)

    # Calculate similarities between genes
    # give the name of the genes
    names(pathways) <- unlist(clusters, use.names = FALSE)
    genesSims <- combineScoresPar(pathSims, method[1L], pathways, ... = ...)

    # Calculate similarities between clusters of genes
    as.matrix(combineScoresPar(genesSims, method[2L], clusters, ... = ...))
}


#' @describeIn mclusterGeneSim Calculates all the similarities of the
#' GeneSetCollection and combine them using \code{\link{combineScoresPar}}
#' @export
setMethod(
    "mclusterGeneSim",
    c(info = "GeneSetCollection", clusters = "list"),
    function(clusters, info, method, ...) {
        if (length(clusters) == 1) {
            stop(
                "Introduce several clusters!\n",
                "If you want to calculate the similarity ",
                "between genes use mgeneSim"
            )
        }

        if (!all(sapply(clusters, is.character))) {
            stop("The input genes should be characters")
        }
        # Remove duplicate genes in each cluster
        clusters <- lapply(clusters, unique)

        # Convert to a list
        info_list <- inverseList(GSEABase::geneIds(info))

        if (any(!unlist(clusters, use.names = FALSE) %in% names(info_list))) {
            clusters <- sapply(clusters, function(x) {
                x[x %in% names(info_list)]
            })
            warning("Some genes are not in the GeneSetCollection provided.")
        }
        if (any(lengths(clusters) == 0)) {
            clusters <- clusters[lengths(clusters) != 0]
            warning("Some clusters are not present in the GeneSetCollection provided")
        }
        if (length(method) != 2) {
            stop("Please provide two methods to combine scores")
        }

        # Extract all pathways for each gene
        pathways <- lapply(unlist(clusters, use.names = FALSE), function(x) {
            info_list[[x]]
        })
        pathwaysl <- unique(unlist(pathways, use.names = FALSE))
        pathwaysl <- pathwaysl[!is.na(pathwaysl)]

        # Calculate similarities between pathways
        pathSims <- mpathSim(pathwaysl, info_list, NULL)

        # Calculate similarities between genes
        # give the name of the genes
        names(pathways) <- unlist(clusters, use.names = FALSE)
        genesSims <- combineScoresPar(pathSims, method[1L], pathways, ... = ...)

        # Calculate similarities between clusters of genes
        as.matrix(combineScoresPar(genesSims, method[2L], clusters, ... = ...))
    }
)
