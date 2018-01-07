IC <- function(x){
    freq <- table(x)/length(x)
    -sum(freq*log2(freq))
}

maxIC <- function(x) {
    tab <- table(x)
    freq <- rep(1/length(tab), length(tab))
    -sum(freq*log2(freq))
}

#' @importMethodsFrom GSEABase collectionType
check_collectionType <- function(object) {
    typeGS <- lapply(object, collectionType)
    cl <- sapply(typeGS, class)
    equal <- length(table(cl))

    if (equal > 1) {
        warning("Several collections origins detected")
    }

    if (any(cl %in% "GOCollection")) {
        warning("Gene Ontologies are not pathways, proceed with caution.")
    }
}

#' @importMethodsFrom GSEABase geneIdType
check_geneIdType <- function(object) {
    typeGS <- lapply(object, geneIdType)
    cl <- sapply(typeGS, class)
    equal <- length(table(cl))

    if (equal > 1) {
        stop("Several gene ID types detected. Use the same identifier for all
             the gene sets.")
    }
}

# Checks that the GeneSetCollection is not from class GeneOntology
# Checks that the identifiers are not mixed as self reported
check_gsc <- function(object) {
    check_collectionType(object)
    check_geneIdType(object)

}


#' Invert a list
#'
#' Calculate the pathways per gene of an incidence matrix
#'
#' @param incidence Incidence matrix
#' @return The number of pathways each gene has
#' @keywords internal
inverseList <- function(x){
    genes <- unlist(x, use.names = FALSE)
    pathways <- rep(names(x), lengths(x))
    split(pathways, genes)
}

#' Calculates H-index
#'
#' Given a table find which is the h-index of the data
#' @param x The table or values to calculate the h-index from it.
#' @return The numeric value of the h.index
#' @keywords internal
h_index <- function(x) {
    if (class(x) != "table") {
        x <- table(x)
    }
    position <- x - as.numeric(names(x))

    return(max(0, as.numeric(names(x)[position >= 0])))

}
