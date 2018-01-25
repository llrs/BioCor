# clusterSim ####
#' Similarity score between clusters of genes based on pathways similarity
#'
#' Looks for the similarity between genes in groups
#'
#' Once the pathways for each cluster are found they are combined using
#' \code{combineScores}. In \code{mclusterSim} the function
#' \code{combineScoresPar} is used instead.
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
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- as.list(org.Hs.egPATH)
#' clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg)
#' clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg, NULL)
#' clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg, "avg")
#' } else {
#' warning('You need org.Hs.eg.db package for this example')
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

vclusterSim <- Vectorize(clusterSim,
                         vectorize.args = c("cluster1", "cluster2"))

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
mclusterSim <- function(clusters, info, method = "max", ...) {

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

    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    if (any(!unlist(clusters, use.names = FALSE) %in% names(info))) {
        warning("Some genes are not in the list provided.")
    }

    if (is.null(method)) {
        method <- "max"
        warning("Method to combine pathways can't be null, set to 'max'")
    }

    # Find the pathways for each cluster
    pathsGenes <- info[unlist(clusters, use.names = FALSE)]
    pathsGenes <- pathsGenes[!is.na(names(pathsGenes))]
    cluster2pathways <- lapply(clusters, function(genes){
        x <- unlist(pathsGenes[genes], use.names = FALSE)
        x[!is.na(x)]})

    pathways <- unique(unlist(cluster2pathways, use.names = FALSE)) # Total pathways
    pathways <- pathways[!is.na(pathways)]

    # Calculates similarities between pathways
    names(pathways) <- pathways

    if (!is.null(pathways)) { # check that there is at least one pathway
        pathSims <- mpathSim(pathways, info, method = NULL)

        # Calculates similarities between clusters
        sim <- combineScoresPar(pathSims, method, cluster2pathways, ... = ...)
    } else {
        sim <- as.matrix(NA)
    }

    # In case any cluster don't have any relevant data
    sim_all <- matrix(NA, ncol = length(clusters), nrow = length(clusters),
                      dimnames = list(names(clusters), names(clusters)))
    AintoB(as.matrix(sim), sim_all)

}
