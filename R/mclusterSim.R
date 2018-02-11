#' Similarity score between clusters of genes based on pathways similarity
#'
#' Looks for the similarity between genes in groups. Once the pathways for each
#' cluster are found they are combined using code{\link{combineScores}}.
#' @param clusters A list of clusters of genes to be found in \code{id}.
#' @inheritParams geneSim
#' @inheritParams combineScores
#' @inheritParams pathSim
#' @export
#' @author Lluís Revilla
#' @seealso For a different approach see \code{\link{clusterGeneSim}},
#' \code{\link{combineScores}} and \code{\link{conversions}}
#' @return \code{mclusterSim} returns a matrix with the similarity scores for
#' each cluster comparison.
#' @examples
#' if (require("org.Hs.eg.db")) {
#'     #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#'     # data of June 31st 2011)
#'     genes.kegg <- as.list(org.Hs.egPATH)
#'
#'     clusters <- list(cluster1 = c("18", "81", "10"),
#'                      cluster2 = c("100", "10", "1"),
#'                      cluster3 = c("18", "10", "83"))
#'     mclusterSim(clusters, genes.kegg)
#'     mclusterSim(clusters, genes.kegg, "avg")
#' } else {
#'     warning('You need org.Hs.eg.db package for this example')
#' }
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


#' @describeIn mclusterSim Calculates all the similarities of the GeneSetCollection
#' and combine them using \code{\link{combineScoresPar}}
#' @export
setMethod("mclusterSim",
          c(info = "GeneSetCollection", clusters = "list"),
          function(clusters, info, method, ...) {
              if (length(clusters) <= 2 & all(lengths(clusters) == 1)) {
                  warnings("Using mgeneSim")
                  return(mgeneSim(clusters, info, method, ...))
              }

              # Check they are unique
              clusters <- lapply(clusters, unique)

              # Extract the ids
              origGenes <- geneIds(info)

              # Check that the genes are in the GeneSetCollection
              genes <- unique(unlist(origGenes, use.names = FALSE))
              if (all(!unlist(clusters, use.names = FALSE) %in% genes)) {
                  warning("At least one gene should be in the list provided")
                  return(NA)
              }

              # Simplify the GeneSetCollection
              keep <- sapply(origGenes, function(x) {
                  keepPaths <- sapply(clusters, function(y){
                      any(y %in% x)
                  })
                  unique(keepPaths[keepPaths])
              })

              gscGenes <- info[names(keep[keep])]

              # Search for the paths of each gene
              ids <- geneIds(gscGenes)
              paths <- lapply(clusters, function(x){
                  keepPaths <- sapply(ids, function(y) {
                      any(x %in% y)
                  })
                  names(keepPaths[keepPaths])
              })

              # Calculate the pathSim of all the implied pathways
              pathsSim <- mpathSim(info = gscGenes, method = NULL)
              # Summarize the information
              combineScoresPar(pathsSim, method, subSets = paths, ...)
          }
)
