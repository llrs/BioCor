# diceSim ####
#' Compare pathways
#'
#' Function to estimate how much two list of genes overlap by looking
#' how much of the nodes are shared. Calculates the Dice similarity
#'
#' It requires a vector of characters otherwise will return an `NA`.
#' @param g1,g2 A character list with the names of the proteins in each pathway.
#' @return A score between 0 and 1 calculated as the double of the proteins
#' shared by g1 and g2 divided by the number of genes in both groups.
#' @export
#' @author Lluís Revilla
#' @seealso Used for [geneSim()], see [conversions()] help
#' page to transform Dice score to Jaccard score.
#' @examples
#' genes.id2 <- c("52", "11342", "80895", "57654", "548953", "11586", "45985")
#' genes.id1 <- c(
#'     "52", "11342", "80895", "57654", "58493", "1164", "1163",
#'     "4150", "2130", "159"
#' )
#' diceSim(genes.id1, genes.id2)
#' diceSim(genes.id2, genes.id2)
diceSim <- function(g1, g2) {
  if (is.character(g1) && is.character(g2)) {
    prot1 <- g1
    prot2 <- g2
  } else if (is.null(g1) || is.null(g2)) { # Check if length is 0
    return(NA)
  }
  else {
    warning("g1 or g2 is not character")
    return(NA)
  }
  # If there isn't any information of a pathway for a gene then then
  # functional similarity is 0
  if (length(prot1) == 0L || length(prot2) == 0L) {
    return(0L)
  }
  score <- (length(intersect(prot1, prot2))) * 2L / (
    length(prot2) + length(prot1))
  score
}

vdiceSim <- Vectorize(diceSim)
