#' @importClassesFrom GSEABase GeneSetCollection
#' @importMethodsFrom GSEABase geneIds
#' @importFrom GSEABase getBroadSets
NULL

#' Pathways per gene
#'
#' Calculates the pathways per gene of a GeneSetCollection
#' @param object A GeneSetCollection object
#' @return a list of pathways per genes
#' @author Lluís Revilla
#' @rdname PathwaysPerGene-GeneSetCollection-method
#' @aliases PathwaysPerGene
#' @exportMethod PathwaysPerGene
setGeneric("PathwaysPerGene", function(object) {
    standardGeneric("PathwaysPerGene")
})

#' @exportMethod PathwaysPerGene
setMethod("PathwaysPerGene",
          "GeneSetCollection",
          function(object) {
              check_gsc(object)
              pathways2genes <- GSEABase::geneIds(object)
              lengths(inverseList(pathways2genes))
          }
)

#' Genes per pathway
#'
#' Calculates the genes per pathway of a GeneSetCollection
#' @param object A GeneSetCollection object
#' @return a list of genes per pathway
#' @author Lluís Revilla
#' @rdname GenesPerPathway-GeneSetCollection-method
#' @aliases GenesPerPathway
#' @exportMethod GenesPerPathway
setGeneric("GenesPerPathway", function(object) {
    standardGeneric("GenesPerPathway")
})

#' @exportMethod GenesPerPathway
setMethod("GenesPerPathway",
          "GeneSetCollection",
          function(object) {
              check_gsc(object)
              lengths(GSEABase::geneIds(object))
          }
)
