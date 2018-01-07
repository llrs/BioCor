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
#' @aliases PathwaysPerGene
#' @exportMethod PathwaysPerGene
#' @seealso \code{\link{GenesPerPathway}}
setGeneric("PathwaysPerGene", function(object) {
    standardGeneric("PathwaysPerGene")
})

#' @describeIn PathwaysPerGene Pathways per gene in the GeneSetCollection
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
#' @return A list of genes per pathway.
#' @author Lluís Revilla
#' @aliases GenesPerPathway
#' @exportMethod GenesPerPathway
#' @seealso \code{\link{PathwaysPerGene}}
setGeneric("GenesPerPathway", function(object) {
    standardGeneric("GenesPerPathway")
})

#' @describeIn GenesPerPathway Calculates genes per pathway of the
#' GeneSetCollection
#' @exportMethod GenesPerPathway
setMethod("GenesPerPathway",
          "GeneSetCollection",
          function(object) {
              check_gsc(object)
              lengths(GSEABase::geneIds(object))
          }
)
