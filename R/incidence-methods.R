#' @importClassesFrom GSEABase GeneSetCollection
#' @importMethodsFrom GSEABase geneIds
#' @importFrom GSEABase getBroadSets
NULL

#' Pathways per gene
#'
#' Calculate the pathways per gene of an incidence matrix
#'
#' @param incidence Incidence matrix
#' @return The number of pathways each gene has
#' @keywords internal
PathwaysPerGene <- function(incidence){
    colSums(incidence)
}

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
#' Calculate the genes per pathway of an incidence matrix
#'
#' @param incidence Incidence matrix
#' @return The number of genes each pathway has
#' @keywords internal
GenesPerPathway <- function(incidence){
    rowSums(incidence)
}


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
