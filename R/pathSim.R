# pathSim ####
#' Calculates the Dice similarity between pathways
#'
#' Calculates the similarity between pathways using dice similarity score.
#'
#' \code{diceSim} is used to calculate similarities between the two pathways.
#'
#' \code{mpathSim} compares the similarity between several pathways and can use
#' \code{\link{combineScores}} to extract the similarity between those
#' pathways.
#' If one needs the matrix of similarities between pathways set the argument
#' methods to \code{NULL}.
#' @param pathway1,pathway2 A single pathway to calculate the similarity
#' @param info A list of genes and the pathways they are involved.
#' @param method To combine the scores of each pathway, one of \code{c("avg",
#' "max", "rcmax", "rcmax.avg", "BMA")}, if NULL returns the matrix of
#' similarities.
#' @param ... Other arguments passed to \code{\link{combineScores}}
#' @return The similarity between those pathways or all the similarities
#' between each comparison.
#' @seealso \code{\link{diceSim}} and \code{\link{combineScores}} and
#' \code{\link{conversions}} help page to transform Dice score to Jaccard
#' score.
#' @author Lluís Revilla
#' @export
#' @examples
#' library("reactome.db")
#' # Extracts the paths of all genes of org.Hs.eg.db from reactome
#' genes.react <- as.list(reactomeEXTID2PATHID)
#' (paths <- sample(unique(unlist(genes.react)), 2))
#' pathSim(paths[1], paths[2], genes.react)
#'
#' (paths <- sample(unique(unlist(genes.react)), 10))
#' mpathSim(pathways, genes.react, NULL)
#' named_paths <- structure(c("R-HSA-112310", "R-HSA-112316", "R-HSA-112315"),
#' .Names = c("Neurotransmitter Release Cycle",
#' "Neuronal System", "Transmission across Chemical Synapses"))
#' mpathSim(named_paths, genes.react, NULL)
pathSim <- function(pathway1, pathway2, info) {
    if (length(pathway1) != 1 | length(pathway2) != 1) {
        stop("Introduce just one pathway!\n",
             "If you want to calculate several similarities ",
             "between pathways use mpathSim")
    }
    if (!is.character(pathway1)  | !is.character(pathway2)) {
        stop("The input pathways should be characters")
    }
    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    if (any(!c(pathway1, pathway2) %in% unlist(info))) {
        return(NA)
    }

    # Invert the list
    rId <- unlist(info, use.names = FALSE)
    lId <- rep(names(info), sapply(info, length))
    pathways2genes <- split(lId, rId)

    # Convert the list
    pathways2genes <- list2env(pathways2genes)

    # Extract the gene ids for each pathway
    g1 <- pathways2genes[[pathway1]]
    g2 <- pathways2genes[[pathway2]]

    diceSim(g1, g2)
}

vpathSim <- Vectorize(pathSim, vectorize.args = c("pathway1", "pathway2"))

#' @rdname pathSim
#' @param pathways Pathways to calculate the similarity for
#' @note pathways accept named characters, and then the output will have
#' the names
#' @export
mpathSim <- function(pathways, info, method = NULL, ...) {

    if (length(unique(pathways)) == 1 ) {
        stop("Introduce several unique pathways!\n",
             "If you want to calculate one similarity ",
             "between pathways use pathSim")
    }

    if (!all(is.character(pathways))) {
        stop("The input pathways should be characters")
    }
    nam <- names(pathways)
    pathways <- unique(pathways)

    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    if (any(!pathways %in% unlist(info))) {
        warning("Some pathways are not in the list you provided.")
    }

    # If the number of pathways is quite big use another strategy
    if (length(pathways) >= 150) {

        sim <- pathSims_matrix(info)
        sim <- sim[pathways, pathways]

    } else {
        # Invert the list
        rId <- unlist(info, use.names = FALSE)
        lId <- rep(names(info), sapply(info, length))
        pathways2genes <- split(lId, rId)

        # Extract the gene ids for each pathway
        g1 <- sapply(pathways, function(x) {
            pathways2genes[[x]]
        }, simplify = FALSE)
        g2 <- g1

        # Calculate similarities
        sim <- outer(g1, g2, vdiceSim)
    }
    if (!is.null(nam)) {
        dimnames(sim) <- list(nam, nam)
    }

    # Calculate the similarity between the two genes
    if (is.null(method)) {
        return(sim)
    } else {
        combineScores(sim, method, ... = ...)
    }
}

# pathSims_matrix ####
# Uses linear algebra to speed the caluclations
# x is a list of genes to pathways
# Omits pathways with no gene
pathSims_matrix <- function(x) {
    nas <- sapply(x, function(y){all(is.na(y))})
    lge2 <- x[!nas]
    pathways <- unique(unlist(lge2))
    mat <- as.matrix(sapply(names(lge2), function(y){
        ifelse(pathways %in% lge2[[y]], TRUE, FALSE)
    }
    ))
    rownames(mat) <- pathways
    overPath <- crossprod(t(mat))
    genesPerPathway <- diag(overPath)
    genesPerPathway <- matrix(genesPerPathway, ncol(overPath), ncol(overPath))
    2*overPath/(t(genesPerPathway) + genesPerPathway)
}
