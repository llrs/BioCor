# clusterSim ####
#' Compare two clusters of genes
#'
#' Looks for the similarity between genes in groups
#'
#' Once the pathways for each cluster are found they are combined using
#' combineScores.
#' @param cluster1,cluster2 A vector with genes in \code{id}.
#' @inheritParams geneSim
#' @inheritParams combineScores
#' @export
#' @author Lluís Revilla
#' @seealso For a different approach see \code{\link{clustersSim}},
#' \code{\link{combineScores}} and \code{\link{conversions}}
#' @return \code{clusterSim} returns a similarity score of the two clusters
#' @examples
#' library("org.Hs.eg.db")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- as.list(org.Hs.egPATH)
#' clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg)
#' clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg, NULL)
#' clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg, "avg")
clusterSim <- function(cluster1, cluster2, info, method = "max"){

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
        warning("Some genes are not in the list you provided.")
    }

    # Extract all pathways for each gene
    pathways1 <- sapply(cluster1, function(x) {
            info[[x]]
        }, simplify = FALSE)
    pathways2 <- sapply(cluster2, function(x) {
        info[[x]]
    }, simplify = FALSE)

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
        combineScores(sim, method)
    } else {
        sim
    }
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
#' mclusterSim(clusters, genes.kegg)
#' mclusterSim(clusters, genes.kegg, "avg")
mclusterSim <- function(clusters, info, method = "max") {

    if (!is.list(cluster)) {
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
    clusters <- sapply(clusters, unique, simplify = FALSE)

    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    if (any(!unlist(clusters) %in% names(info))) {
        warning("Some genes are not in the list you provided.")
    }

    if (is.null(method)) {
        method <- "max"
        warning("Method to combine pathways can't be null, set to 'max'")
    }

    pathways <- unique(unlist(sapply(unlist(clusters), getElement,
                                     object = info)))

    # Depending how big the pathways are we do one or other strategy
    if (sum(!is.na(pathways)) >= 30) {
        # Using precalculated pathway similarities
        pathSim <- pathSims_matrix(info)

        nas <- sapply(info, function(y){all(is.na(y))})
        lge2 <- info[!nas]
        sim <- outer(lge2, lge2, vcombineScoresPrep, prep = pathSim,
                     method = method)
    } else {

        names(pathways) <- pathways
        sim <- outer(pathways, pathways, vgeneSim, info, method = method)
    }

    sim_all <- matrix(NA, ncol = length(clusters), nrow = length(clusters),
                      dimnames = list(names(clusters), names(clusters)))
    AintoB(sim, sim_all)

    # vc <- Vectorize(clusterSim, vectorize.args = c("cluster1", "cluster2"))
    # suppressWarnings(outer(clusters, clusters, vc, info, method = method))

}
