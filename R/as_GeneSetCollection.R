#' @importClassesFrom GSEABase GeneSetCollection
#' @importClassesFrom GSEABase GeneSet
#' @importMethodsFrom GSEABase GeneSet GeneSetCollection
NULL

#' @exportMethod as.GeneSetCollection
setGeneric("as.GeneSetCollection", function(object) {
    standardGeneric("as.GeneSetCollection")
})


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
