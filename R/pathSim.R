# pathSim ####
#' Calculates the Dice similarity between pathways
#'
#' Calculates the similarity between pathways using dice similarity score.
#' \code{diceSim} is used to calculate similarities between the two pathways.
#' @param pathway1,pathway2 A single pathway to calculate the similarity
#' @param info A GeneSetCollection or a list of genes and the pathways they are
#' involved.
#' @return The similarity between those pathways or all the similarities
#' between each comparison.
#' @seealso
#' \code{\link{conversions}} help page to transform Dice score to Jaccard
#' score.
#' \code{\link{mpathSim}} for multiple pairwise comparison of pathways.
#' @author Lluís Revilla
#' @export
#' @examples
#' if (require("reactome.db")){
#'     # Extracts the paths of all genes of org.Hs.eg.db from reactome
#'     genes.react <- as.list(reactomeEXTID2PATHID)
#'     (paths <- sample(unique(unlist(genes.react)), 2))
#'     pathSim(paths[1], paths[2], genes.react)
#' } else {
#'     warning('You need reactome.db package for this example')
#' }
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
        stop("info should be a list or a GeneSetCollection.\n
             See documentation.")
    }

    if (any(!c(pathway1, pathway2) %in% unlist(info, use.names = FALSE))) {
        return(NA)
    }

    # Invert the list
    rId <- unlist(info, use.names = FALSE)
    lId <- rep(names(info), lengths(info))
    pathways2genes <- split(lId, rId)

    # Convert the list
    pathways2genes <- list2env(pathways2genes)

    # Extract the gene ids for each pathway
    g1 <- pathways2genes[[pathway1]]
    g2 <- pathways2genes[[pathway2]]

    diceSim(g1, g2)
}

#' @describeIn pathSim Calculates all the similarities of a GeneSetCollection
#' and combine them using \code{combineScoresPar}
#' @export
setMethod("pathSim",
          c(info = "GeneSetCollection", pathway1 = "character",
            pathway2 = "character"),
          function(pathway1, pathway2, info) {
              if (length(pathway1) != 1 | length(pathway2) != 1) {
                  stop("Introduce just one pathway!\n",
                       "If you want to calculate several similarities ",
                       "between pathways use mpathSim")
              }
              if (any(!c(pathway1, pathway2) %in% names(info))) {
                  return(NA)
              }
              genes <- geneIds(info[c(pathway1, pathway2)])
              diceSim(genes[[1]], genes[[2]])
          }
)
