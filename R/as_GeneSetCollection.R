


#' @exportMethod as.GeneSetCollection
setGeneric("as.GeneSetCollection", function(object) {
    standardGeneric("as.GeneSetCollection")
})


# Table of Length pathways vs Number of pathways a gene is in

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
