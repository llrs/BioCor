# geneSim ####
#' Similarity score genes based on pathways similarity
#'
#' Given two genes, calculates the Dice similarity between each pathway
#' which is combined to obtain a similarity between the genes.
#'
#' Given the information about the genes and their pathways, uses the ids
#' of the genes to find the Dice similarity score for each pathway comparison
#' between the genes. Later this similarities are combined using
#' \code{\link{combineScoresPar}}.
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
#' \code{\link{combineScoresPar}}.
#' @examples
#' if (require("org.Hs.eg.db") & require("reactome.db")) {
#'     # Extract the paths of all genes of org.Hs.eg.db from KEGG
#'     # (last update in data of June 31st 2011)
#'     genes.kegg <- as.list(org.Hs.egPATH)
#'     # Extracts the paths of all genes of org.Hs.eg.db from reactome
#'     genes.react <- as.list(reactomeEXTID2PATHID)
#'     geneSim("81", "18", genes.react)
#'     geneSim("81", "18", genes.kegg)
#'     geneSim("81", "18", genes.react, NULL)
#'     geneSim("81", "18", genes.kegg, NULL)
#' } else {
#'     warning('You need reactome.db and org.Hs.eg.db package for this example')
#' }
geneSim <- function(gene1, gene2, info, method = "max", ...) {

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
        return(NA)
    }

    comb <- c(gene1, gene2)

    if (any(is.na(comb))) {
        return(NA)
    }

    # Extract all pathways for each gene
    pathways <- sapply(comb, function(x) {
        y <- info[[x]]
        y[!is.na(y)]
    }, simplify = FALSE)

    # Check that we have pathways info for this combination
    if (any(lengths(pathways) == 0L)) {
        return(NA)
    }
    # Subseting just the important pathways
    pathways_all <- unique(unlist(pathways, use.names = FALSE))


    sim <- mpathSim(pathways_all, info = info, method = NULL, ...)
    sim <- sim[pathways[[1]], pathways[[2]], drop = FALSE]

    # Combine or not
    if (is.null(method)) {
        sim
    }
    else {
        combineScoresPar(sim, method = method, ...)
    }
}

vgeneSim <- Vectorize(geneSim, vectorize.args = c("gene1", "gene2"))

#' @rdname geneSim
#' @export
#' @param genes A vector of genes.
#' @return \code{mgeneSim} returns the matrix of similarities between the genes
#' in the vector
#' @note genes accept named characters and the output will use the names of the
#' genes.
#' @examples
#'
#' mgeneSim(c("81", "18", "10"), genes.react)
#' mgeneSim(c("81", "18", "10"), genes.react, "avg")
#' named_genes <- structure(c("81", "18", "10"),
#' .Names = c("ACTN4", "ABAT", "NAT2"))
#' mgeneSim(named_genes, genes.react, "max")
mgeneSim <- function(genes, info, method = "max", ...) {

    if (length(unique(genes)) == 1) {
        stop("Introduce several unique genes!\n",
             "If you want to calculate one similarity ",
             "between pathways use geneSim")
    }
    if (!all(is.character(genes))) {
        stop("The input genes should be characters")
    }
    namgenes <- names(genes)
    genes <- unique(genes)

    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    if (all(!genes %in% names(info))) {
        stop("Check genes are in the list provided.")
    } else if (any(!genes %in% names(info))) {
        warning("Some genes are not in the list provided.")
    }

    if (is.null(method)) {
        method <- "max"
        warning("Method to combine pathways can't be null, set to 'max'")
    }

    pathways <- info[names(info) %in% genes]
    pathwaysl <- unique(unlist(pathways, use.names = FALSE))
    pathwaysl <- pathwaysl[!is.na(pathwaysl)]

    pathsSims <- mpathSim(pathwaysl, info, NULL)
    sim <- combineScoresPar(pathsSims, method, pathways, ... = ...)

    sim_all <- matrix(NA, ncol = length(genes), nrow = length(genes),
           dimnames = list(genes, genes))
    sim <- AintoB(as.matrix(sim), sim_all)
    if (!is.null(namgenes)) {
        if (length(namgenes) != nrow(sim)) {
            warning("Omitting gene names: duplicated names")
        } else {
            dimnames(sim) <- list(namgenes, namgenes)
        }
    }
    sim


}

