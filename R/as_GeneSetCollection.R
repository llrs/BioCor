#' @importClassesFrom GSEABase GeneSetCollection
#' @importClassesFrom GSEABase GeneSet
#' @importMethodsFrom GSEABase GeneSet GeneSetCollection
NULL

#' Convert a list or a Bimap interface into GeneSetCollections
#'
#' @param object A list of genes and their pathways or an AnnDbBimap.
#' @return A GeneSetCollection
#' @author Lluís Revilla
#' @exportMethod as.GeneSetCollection
#' @seealso \code{\link{GeneSetCollection}}
setGeneric("as.GeneSetCollection", function(object)
    standardGeneric("as.GeneSetCollection")
)

#' @describeIn as.GeneSetCollection Convert list to GeneSetCollections
#' @exportMethod as.GeneSetCollection
setMethod("as.GeneSetCollection",
          signature(object = "list"),
          function(object) {
              paths2gene <- inverseList(object)
              gsl <- sapply(names(paths2gene), function(x){
                  GeneSet(unique(paths2gene[[x]]), setName = x)})

              GeneSetCollection(gsl)
          }
)

#' @describeIn as.GeneSetCollection Convert AnnDbBimap to GeneSetCollections
#' @exportMethod as.GeneSetCollection
setMethod("as.GeneSetCollection",
          signature(object = "AnnDbBimap"),
          function(object) {
              paths2gene <- inverseList(as.list(object))
              gsl <- sapply(names(paths2gene), function(x){
                  GeneSet(unique(paths2gene[[x]]), setName = x)})

              GeneSetCollection(gsl)
          }
)

setAs("list", "GeneSetCollection",
      function(from)
          as.GeneSetCollection(from)
)
