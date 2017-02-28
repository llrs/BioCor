# pathSim ####
#' Compare pathways
#'
#' Given two vectors of pathways calculates the similarity between them.
#'
#' \code{diceSim} is used to calculate similarities between each pathway and
#' combineScores to extract the similarity between those pathways. If one need
#' the matrix of similarities set methods to \code{NULL}.
#' @param pathways1,pathways2 Pathways to be found in
#' \code{pathwayDB}.
#' @inheritParams geneSim
#' @return The similarity between those pathways or all the similarities
#' between each comparison.
#' @seealso \code{\link{pathSim}} and \code{\link{combineScores}}
#' @author Lluís Revilla
#' @export
#' @examples
#' library("org.Hs.eg.db")
#' library("reactome.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' # Extracts the paths of all genes of org.Hs.eg.db from reactome
#' genes.react <- select(reactome.db, keys = entrezids, keytype = "ENTREZID",
#'                       columns = "REACTOMEID")
#' pathways1 <- c("112310", "112315", "112316", "373753", "916853")
#' pathways2 <- c("109582", "114608", "1500931", "888590", "76002", "76005")
#' pathSim(pathways1[1], pathways2[1], genes.react, "ENTREZID", "REACTOMEID")
#' mpathSim(pathways1, pathways2, genes.react, "ENTREZID", "REACTOMEID", NULL)
pathSim <- function(pathways1, pathways2, genes, id, pathwayDB) {

    # Convert data.frame into environment to speed the look up
    pathways2genes <- split(genes[ , id], genes[, pathwayDB])
    pathways2genes <- list2env(pathways2genes)

    # Extract the gene ids for each pathway
    g1 <- pathways2genes[[pathways1]]
    g2 <- pathways2genes[[pathways2]]

    diceSim(g1, g2)
}

#' @rdname pathSim
#' @export
mpathSim <- function(pathways1, pathways2, genes, id, pathwayDB,
                     method = "max") {

    # Convert data.frame into environment to speed the look up
    pathways2genes <- split(genes[ , id], genes[, pathwayDB])
    pathways2genes <- list2env(pathways2genes)

    # Extract the gene ids for each pathway
    g1 <- sapply(pathways1, function(x) {
        pathways2genes[[x]]
    }, simplify = FALSE)
    g2 <- sapply(pathways2, function(x) {
        pathways2genes[[x]]
    }, simplify = FALSE)


    vp <- Vectorize(diceSim)
    react <- outer(g1, g2, vp)

    # Calculate the similarity between the two genes
    if (is.null(method)) {
        react
    } else {
        combineScores(react, method)
    }
}
# pathSims_precalculate ####
pathSims_precalculate <- function(x) {
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

pairN <- function(x, y, prep, method){
    combineScores(prep[x, y, drop = FALSE], method)
}
vpairN <- Vectorize(pairN, vectorize.args = c("x", "y"))
