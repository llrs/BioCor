# geneSim ####
#' Similarity score genes based on pathways similarity
#'
#' Given two genes, calculates the Dice similarity between each pathway
#' which is combined to obtain a similarity between the genes.
#'
#' Given the information about the genes and their pathways, uses the ids
#' of the genes to find the Dice similarity score for each pathway comparison
#' between the genes. Later this similarities are combined using
#' [combineScoresPar()].
#' @param gene1,gene2 Ids of the genes to calculate the similarity, to be found
#' in genes.
#' @inheritParams pathSim
#' @inheritParams combineScores
#' @return
#' The highest Dice score of all the combinations of pathways between
#' the two ids compared if a method to combine scores is provided or NA if
#' there isn't information for one gene.
#' If an `NA` is returned this means that there isn't information
#' available for any pathways for one of the genes. Otherwise a number
#' between 0 and 1 (both included) is returned. Note that there isn't a
#' negative value of similarity.
#' @export
#' @author Lluís Revilla
#' @seealso [mgeneSim()], [conversions()] help page to transform Dice
#' score to Jaccard score. For the method to combine the scores see
#' [combineScoresPar()].
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
#'     warning("You need reactome.db and org.Hs.eg.db package for this example")
#' }
geneSim <- function(gene1, gene2, info, method = "max", ...) {
    if (length(unique(gene1)) != 1L | length(unique(gene2)) != 1L) {
        stop(
            "Introduce just one gene!\n",
            "If you want to calculate several similarities ",
            "between genes use mgeneSim"
        )
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
    pathways <- lapply(comb, function(x) {
        y <- info[[x]]
        y[!is.na(y)]
    })
    names(pathways) <- comb

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


#' @describeIn geneSim Calculates all the similarities of the GeneSetCollection
#' and combine them using [combineScoresPar()]
#' @export
setMethod(
    "geneSim",
    c(
        info = "GeneSetCollection", gene1 = "character",
        gene2 = "character"
    ),
    function(gene1, gene2, info, method, ...) {
        if (length(gene1) != 1 | length(gene2) != 1) {
            stop(
                "Introduce just one gene!\n",
                "If you want to calculate several similarities ",
                "between genes use mgeneSim"
            )
        }
        # Extract the ids
        origGenes <- geneIds(info)
        # Check that the genes are in the GeneSetCollection
        genes <- unique(unlist(origGenes, use.names = FALSE))
        if (any(!c(gene1, gene2) %in% genes)) {
            return(NA)
        }
        # Simplify the GeneSetCollection
        keep <- sapply(origGenes, function(x) {
            any(c(gene1, gene2) %in% x)
        })
        gscGenes <- info[names(keep[keep])]

        # Search for the paths of each gene
        paths <- sapply(c(gene1, gene2), function(x) {
            keepPaths <- sapply(geneIds(gscGenes), function(y) {
                any(x %in% y)
            })
            names(keepPaths[keepPaths])
        })

        # Calculate the pathSim of all the implied pathways
        pathsSim <- mpathSim(info = gscGenes, method = NULL)

        # Summarize the information
        if (is.null(method)) {
            pathsSim[paths[[1]], paths[[2]], drop = FALSE]
        } else {
            out <- combineScoresPar(pathsSim, method, subSets = paths, ...)
            out[gene1, gene2]
        }
    }
)
