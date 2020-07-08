#' Combine different sources of pathways
#'
#' Given several sources of pathways with the same for the same id of the genes
#' it merge them.
#'
#' It assumes that the identifier of the genes are the same for both sources
#' but if many aren't equal it issues a warning. Only unique pathways
#' identifiers are returned.
#' @param ... Lists of genes and their pathways.
#' @return A single list with the pathways of each source on the same gene.
#' @importFrom stats setNames
#' @examples
#' DB1 <- list(g1 = letters[6:8], g2 = letters[1:5], g3 = letters[4:7])
#' DB2 <- list(
#'     g1 = c("one", "two"), g2 = c("three", "four"),
#'     g3 = c("another", "two")
#' )
#' combineSources(DB1, DB2)
#' combineSources(DB1, DB1)
#' DB3 <- list(
#'     g1 = c("one", "two"), g2 = c("three", "four"),
#'     g4 = c("five", "six", "seven"), g5 = c("another", "two")
#' )
#' combineSources(DB1, DB3) # A warning is expected
#' @export
combineSources <- function(...) {
    sources <- list(...)
    geneIds <- lapply(sources, names)
    uGeneIds <- unique(unlist(geneIds, use.names = FALSE))

    # Compare the percentatge of uniques of each one
    uniquess <- sapply(sources, function(x) {
        nam <- names(x)
        length(setdiff(uGeneIds, nam)) / length(uGeneIds)
    })

    if (any(uniquess >= 0.5)) {
        warning(
            "More than 50% of genes identifiers of a source are unique\n",
            "Check the identifiers of the genes"
        )
    } else if (any(uniquess >= 0.25)) {
        warning(
            "More than 25% of genes identifiers of a source are unique\n",
            "Check the identifiers of the genes"
        )
    }
    # see http://stackoverflow.com/a/18539199/2886003
    out <- setNames(do.call(mapply, c(
        FUN = c,
        lapply(sources, `[`, uGeneIds)
    )), uGeneIds)
    sapply(out, unique)
}
