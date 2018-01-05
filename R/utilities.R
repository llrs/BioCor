IC <- function(x){
    freq <- table(x)/length(x)
    -sum(freq*log2(freq))
}

#' @importMethodsFrom GSEABase collectionType
#' @importMethodsFrom GSEABase GeneIdentifierType
check_collectionType <- function(object) {
    typeGS <- lapply(object, collectionType)
    cl <- sapply(typeGS, class)
    equal <- lenth(table(cl))

    if (equal > 1) {
        warning("Several collections origins detected")
    }

    if (any(cl %in% "GOCollection")) {
        warning("Gene Ontologies are not pathways, proceed with caution.")
    }
}

check_geneIdType <- function(object) {
    typeGS <- lapply(object, geneIdType)
    cl <- sapply(typeGS, class)
    equal <- lenth(table(cl))

    if (equal > 1) {
        stop("Several gene ID types detected. Use the same identifier for all
             the gene sets.")
    }
}

# Checks that the GeneSetCollection is not from class GeneOntology
# Checks that the identifiers are not mixed as self reported
#' @importMethodsFrom GSEABase collectionType
#' @importMethodsFrom GSEABase GeneIdentifierType
check_gsc <- function(object) {
    check_collectionType(object)
    check_geneIdType(object)

}
