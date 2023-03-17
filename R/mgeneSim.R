#' Similarity score genes based on pathways similarity
#'
#' Given two genes, calculates the Dice similarity between each pathway
#' which is combined to obtain a similarity between the genes.
#'
#' Given the information about the genes and their pathways, uses the ids
#' of the genes to find the Dice similarity score for each pathway comparison
#' between the genes. Later this similarities are combined using
#' \code{\link{combineScoresPar}}.
#' @inheritParams pathSim
#' @inheritParams combineScores
#' @export
#' @param genes A vector of genes.
#' @return \code{mgeneSim} returns the matrix of similarities between the genes
#' in the vector
#' @note genes accept named characters and the output will use the names of the
#' genes.
#' @seealso \code{\link{geneSim}}, \code{\link{conversions}} help page to transform Dice
#' score to Jaccard score. For the method to combine the scores see
#' \code{\link{combineScoresPar}}.
#' @examples
#' if (require("org.Hs.eg.db") & require("reactome.db")) {
#'     # Extract the paths of all genes of org.Hs.eg.db from KEGG
#'     # (last update in data of June 31st 2011)
#'     genes.kegg <- as.list(org.Hs.egPATH)
#'     # Extracts the paths of all genes of org.Hs.eg.db from reactome
#'     genes.react <- as.list(reactomeEXTID2PATHID)
#'     mgeneSim(c("81", "18", "10"), genes.react)
#'     mgeneSim(c("81", "18", "10"), genes.react, "avg")
#'     named_genes <- structure(c("81", "18", "10"),
#'         .Names = c("ACTN4", "ABAT", "NAT2")
#'     )
#'     mgeneSim(named_genes, genes.react, "max")
#' } else {
#'     warning("You need reactome.db and org.Hs.eg.db package for this example")
#' }
mgeneSim <- function(genes, info, method = "max", ...) {
    if (length(unique(genes)) == 1) {
        stop(
            "Introduce several unique genes!\n",
            "If you want to calculate one similarity ",
            "between pathways use geneSim"
        )
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

    sim_all <- matrix(NA,
        ncol = length(genes), nrow = length(genes),
        dimnames = list(genes, genes)
    )
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

#' @describeIn mgeneSim Calculates all the similarities of the list and
#' combine them using \code{\link{combineScoresPar}}
#' @export
setMethod(
    "mgeneSim",
    c(info = "GeneSetCollection", genes = "character"),
    function(genes, info, method, ...) {
        if (length(genes) < 2) {
            stop(
                "Introduce several genes!\n",
                "If you want to calculate one similarities ",
                "between genes use geneSim"
            )
        }

        if (length(unique(genes)) == 1) {
            stop(
                "Introduce several unique genes!\n",
                "If you want to calculate one similarity ",
                "between pathways use geneSim"
            )
        }

        # Extract the ids
        origGenes <- geneIds(info)
        # Check that the genes are in the GeneSetCollection
        genesU <- unique(unlist(origGenes, use.names = FALSE))

        if (any(!genes %in% genesU)) {
            return(NA)
        }
        # Simplify the GeneSetCollection
        keep <- sapply(origGenes, function(x) {
            any(genes %in% x)
        })
        gscGenes <- info[names(keep[keep])]

        ids <- origGenes[keep]
        # Search for the paths of each gene
        paths <- lapply(genes, function(x) {
            keepPaths <- sapply(ids, function(y) {
                any(x %in% y)
            })
            names(keepPaths[keepPaths])
        })
        if (is.null(names(genes))) {
            names(paths) <- genes
        } else {
            names(paths) <- names(genes)
        }
        # Calculate the pathSim of all the implied pathways
        pathsSim <- mpathSim(info = gscGenes, method = NULL)
        # Summarize the information
        combineScoresPar(pathsSim, method, subSets = paths, ...)
    }
)


#' @describeIn mgeneSim Calculates all the similarities of the list and
#' combine them using \code{\link{combineScoresPar}}
#' @export
setMethod(
    "mgeneSim",
    c(info = "GeneSetCollection", genes = "missing"),
    function(info, method, ...) {

        # Extract the name of the pathways per gene
        paths <- inverseList(GSEABase::geneIds(info))

        # Calculate the pathSim of all the implied pathways
        pathsSim <- mpathSim(info = info, method = NULL)

        # Summarize the information
        combineScoresPar(pathsSim, method, subSets = paths, ...)
    }
)
