# geneSim ####
#' Calculates the Dice similarity score of two genes
#'
#' For the genes given calculates the Dice similarity between each pathway
#' which is combined to obtain a similarity between the genes.
#'
#' Given the information about the genes and their pathways, uses the ids
#' of the genes to find the Dice similarity score for each pathway comparison
#' between the genes. Later this similarities are combined using
#' \code{\link{combineScores}}.
#' @param gene1,gene2 Ids of the genes to calculate the similarity, to be found
#' in genes.
#' @inheritParams pathSim
#' @inheritParams combineScores
#' @return
#' The highest Dice score of all the combinations of pathways between
#' the two ids compared if a method to combine scores is provided or NA if
#' there isn't information for one gene.
#' If an \code{NA} is returned this means that there isn't information
#' available for any pathways for one of the genes. Otherwise a number
#' between 0 and 1 (both included) is returned. Note that there isn't a
#' negative value of similarity.
#' @export
#' @author Lluis Revilla
#' @seealso \code{\link{conversions}} help page to transform Dice
#' score to Jaccard score. For the method to combine the scores see
#' \code{\link{combineScores}}.
#' @examples
#' library("org.Hs.eg.db")
#' library("reactome.db")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- as.list(org.Hs.egPATH)
#' # Extracts the paths of all genes of org.Hs.eg.db from reactome
#' genes.react <- as.list(reactomeEXTID2PATHID)
#' geneSim("81", "18", genes.react)
#' geneSim("81", "18", genes.kegg)
#' geneSim("81", "18", genes.react, NULL)
#' geneSim("81", "18", genes.kegg, NULL)
geneSim <- function(gene1, gene2, info, method = "max") {

    if (length(unique(gene1)) != 1L | length(unique(gene2)) != 1L) {
        stop("Introduce just one gene!\n",
             "If you want to calculate several similarities ",
             "between genes use mgeneSim")
    }

    if (!is.character(gene1) | !is.character(gene2)) {
        stop("Please use character")
    }

    if (!is.list(info)) {
        stop("Please introduce info as a list")
    }

    if (any(!c(gene1, gene2) %in% names(info))) {
        stop("A gene is not in the list you provided.")
    }

    comb <- c(gene1, gene2)

    if (any(is.na(comb))) {
        return(NA)
    }
    if (gene1 == gene2) {
        return(1)
    }

    # Convert list into environment to speed the look up
    genes2pathways <- list2env(info)

    # Extract all pathways for each gene
    pathways <- sapply(comb, function(x) {
        genes2pathways[[x]]
    }, simplify = TRUE)

    # Check that we have pathways info for this combination
    if (any(lengths(pathways) == 0L)) {
        return(NA)
    }

    # Subseting just the important pathways
    pathways_all <- unique(unlist(pathways))
    sim <- mpathSim(pathways_all, info = info, method = NULL)
    sim <- sim[pathways[[1]], pathways[[2]], drop = FALSE]

    # Combine or not
    if (is.null(method)) {
        sim
    }
    else {
        combineScores(sim, method = method)
    }
}

vgeneSim <- Vectorize(geneSim, vectorize.args = c("gene1", "gene2"))

#' @rdname geneSim
#' @export
#' @param genes A vector of genes.
#' @return \code{mgeneSim} returns the matrix of similarities between the genes
#' in the vector
#' @examples
#'
#' mgeneSim(c("81", "18", "10"), genes.react)
#' mgeneSim(c("81", "18", "10"), genes.react, "avg")
mgeneSim <- function(genes, info, method = "max") {

    if (length(unique(genes)) == 1) {
        stop("Introduce several unique genes!\n",
             "If you want to calculate one similarity ",
             "between pathways use geneSim")
    }
    if (!all(is.character(genes))) {
        stop("The input genes should be characters")
    }

    genes <- unique(genes)

    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    if (any(!genes %in% names(info))) {
        warning("Some genes are not in the list you provided.")
    }

    if (is.null(method)) {
        method <- "max"
        warning("Method to combine pathways can't be null, set to 'max'")
    }
    pathways <- unique(unlist(sapply(genes, getElement, object = info)))

    # Depending how big the pathways are we do one or other strategy
    if (sum(!is.na(pathways)) >= 30) { #TODO Improve this section not correctly done?
        # Using precalculated pathway similarities
        pathSim <- pathSims_matrix(info)

        nas <- sapply(info, function(y) {all(is.na(y))})
        lge2 <- info[!nas]
        pathsGenes <- sapply(genes, getElement, object = lge2)
        sim <- outer(pathsGenes, pathsGenes, vcombineScoresPrep, prep = pathSim,
                     method = method)
    } else {
        names(genes) <- genes # To have the output with names
        sim <- outer(genes, genes, vgeneSim, info, method = method)
    }

    sim_all <- matrix(NA, ncol = length(genes), nrow = length(genes),
           dimnames = list(genes, genes))
    AintoB(sim, sim_all)

}

