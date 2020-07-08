# mpathSim ####
#' Calculates the Dice similarity between pathways
#'
#' Calculates the similarity between several pathways using dice similarity score.
#' If one needs the matrix of similarities between pathways set the argument
#' methods to \code{NULL}.
#' @param pathways Pathways to calculate the similarity for
#' @param info A list of genes and the pathways they are involved or a
#' GeneSetCollection object
#' @param method To combine the scores of each pathway, one of \code{c("avg",
#' "max", "rcmax", "rcmax.avg", "BMA")}, if NULL returns the matrix of
#' similarities.
#' @param ... Other arguments passed to \code{\link{combineScoresPar}}
#' @return The similarity between those pathways or all the similarities
#' between each comparison.
#' @note \code{pathways} accept named characters, and then the output will have
#' the names
#' @seealso \code{\link{pathSim}} For single pairwise comparison.
#' \code{\link{conversions}} To convert the Dice similarity to Jaccard similarity
#' @export
#' @examples
#' if (require("reactome.db")) {
#'     genes.react <- as.list(reactomeEXTID2PATHID)
#'     (pathways <- sample(unique(unlist(genes.react)), 10))
#'     mpathSim(pathways, genes.react, NULL)
#'     named_paths <- structure(
#'         c("R-HSA-112310", "R-HSA-112316", "R-HSA-112315"),
#'         .Names = c(
#'             "Neurotransmitter Release Cycle",
#'             "Neuronal System",
#'             "Transmission across Chemical Synapses"
#'         )
#'     )
#'     mpathSim(named_paths, genes.react, NULL)
#'     many_pathways <- sample(unique(unlist(genes.react)), 152)
#'     mpathSim(many_pathways, genes.react, "avg")
#' } else {
#'     warning("You need reactome.db package for this example")
#' }
mpathSim <- function(pathways, info, method = NULL, ...) {
    if (length(unique(pathways)) == 1) {
        stop(
            "Introduce several unique pathways!\n",
            "If you want to calculate one similarity ",
            "between pathways use pathSim"
        )
    }

    if (!all(is.character(pathways))) {
        stop("The input pathways should be characters")
    }
    nam <- names(pathways)
    pathways <- unique(pathways)

    if (!is.list(info)) {
        stop("info should be a list. See documentation.")
    }

    if (any(!pathways %in% unlist(info, use.names = FALSE))) {
        warning("Some pathways are not in the list provided.")
    }

    # If the number of pathways is quite big uses matrix properties
    # Calculate just the pathways needed
    if (length(pathways) >= 150) {

        # Keep only the pathways of interest
        pathways2genes <- inverseList(info)
        keep <- pathways %in% unique(unlist(info, use.names = FALSE))
        info <- inverseList(pathways2genes[pathways[keep]])

        sim <- pathSims_matrix(info)
    } else {
        # Invert the list
        pathways2genes <- inverseList(info)

        # Extract the gene ids for each pathway
        g1 <- lapply(pathways, function(x) {
            pathways2genes[[x]]
        })
        names(g1) <- pathways
        g2 <- g1

        # Calculate similarities
        sim <- outer(g1, g2, vdiceSim)
    }

    if (!is.null(nam)) {
        if (length(nam) != nrow(sim)) {
            warning("Omitting pathway names: duplicated names")
        } else {
            dimnames(sim) <- list(nam, nam)
        }
    }

    # Calculate the similarity between the two genes
    if (is.null(method)) {
        return(sim)
    } else {
        combineScoresPar(sim, method, ... = ...)
    }
}


#' Creates the incidence matrix
#'
#' Given a list of pathways and its genes creates an incidence matrix.
#' @note Designed to be easier to work with list and GeneSetCollection
#' @param x A list
#' @return A matrix with pathways as rows and genes in columns.
#' @author Lluís Revilla
#' @keywords internal
setMethod(
    "incidence",
    signature(x = "list"),
    function(x) {
        # Remove empty genes
        nas <- sapply(x, function(y) {
            all(is.na(y))
        })
        lge2 <- x[!nas]
        # Extract all pathways
        pathways <- unique(unlist(lge2, use.names = FALSE))
        # Create the incidence matrix
        mat <- as.matrix(sapply(names(lge2), function(y) {
            ifelse(pathways %in% lge2[[y]], TRUE, FALSE)
        }))
        rownames(mat) <- pathways
        mat
    }
)

# pathSims_matrix ####
# Uses linear algebra to speed the caluclations
# x is a list of genes to pathways or a GeneSetCollection
# Omits pathways with no gene
#' @importMethodsFrom GSEABase incidence
#' @import GSEABase
#' @keywords internal
pathSims_matrix <- function(x) {
    mat <- incidence(x)

    # Calculate genes in common between pathways
    overPath <- tcrossprod(mat)
    # Extract the genes per pathway
    genesPerPathway <- rowSums(mat)
    genesPerPathway <- matrix(genesPerPathway, ncol(overPath), ncol(overPath))
    # Calculate the dice similarity
    2 * overPath / (t(genesPerPathway) + genesPerPathway)
}

#' @describeIn mpathSim Calculates the similarity between the provided pathways
#' of the GeneSetCollection using \code{combineScoresPar}
#' @export mpathSim
setMethod(
    "mpathSim",
    c(info = "GeneSetCollection", pathways = "character", method = "ANY"),
    function(pathways, info, method = NULL, ...) {
        if (length(unique(pathways)) == 1) {
            stop(
                "Introduce several unique pathways!\n",
                "If you want to calculate one similarity ",
                "between pathways use pathSim"
            )
        }

        if (!all(is.character(pathways))) {
            stop("The input pathways should be characters")
        }
        nam <- names(pathways)
        pathways <- unique(pathways)


        if (any(!pathways %in% names(info))) {
            warning("Some pathways are not in the GeneSetCollection provided.")
            m <- matrix(
                nrow = length(pathways), ncol = length(pathways),
                dimnames = list(pathways, pathways)
            )
            pathways <- pathways[pathways %in% names(info)]

            sim <- pathSims_matrix(info[pathways])

            sim <- AintoB(sim, m)
        } else {
            sim <- pathSims_matrix(info[pathways])
        }



        if (is.null(method)) {
            return(sim)
        } else {
            combineScoresPar(sim, method, ... = ...)
        }
    }
)

#' @describeIn mpathSim Calculates all the similarities of the
#' GeneSetCollection and combine them using \code{combineScoresPar}
#' @export mpathSim
setMethod(
    "mpathSim",
    c(info = "GeneSetCollection", pathways = "missing"),
    function(pathways, info, method = NULL, ...) {
        sim <- pathSims_matrix(info)

        if (is.null(method)) {
            return(sim)
        } else {
            combineScoresPar(sim, method, ... = ...)
        }
    }
)

#' @describeIn mpathSim Calculates all the similarities of the list and
#' combine them using \code{combineScoresPar}
#' @export mpathSim
setMethod(
    "mpathSim",
    c(info = "list", pathways = "missing"),
    function(pathways, info, method = NULL, ...) {
        sim <- pathSims_matrix(info)

        if (is.null(method)) {
            return(sim)
        } else {
            combineScoresPar(sim, method, ... = ...)
        }
    }
)

#' @describeIn mpathSim Calculates all the similarities of the list
#' @export mpathSim
setMethod(
    "mpathSim",
    c(info = "list", pathways = "missing", method = "missing"),
    function(pathways, info, method = NULL, ...) {
        sim <- pathSims_matrix(info)

        return(sim)
    }
)
