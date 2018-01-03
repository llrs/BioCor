#' @importClassesFrom GSEABase GeneSetCollection
NULL

#' Calculate the pathways per gene of an incidence matrix
#'
#'
#' @param incidence Incidence matrix
#' @return The number of pathways each gene has
PathwayPerGene <- function(incidence){
    colSums(incidence)
}

#' @exportMethod PathwayPerGene
setGeneric("PathwayPerGene", function(object) {
    standardGeneric("PathwayPerGene")
})

#' @exportMethod PathwayPerGene
setMethod("PathwayPerGene",
          "GeneSetCollection",
          function(object) {
            colSums(GSEABase::incidence(object))
          }
)

#' Calculate the genes per pathway of an incidence matrix
#'
#'
#' @param incidence Incidence matrix
#' @return The number of genes each pathway has
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
              rowSums(GSEABase::incidence(object))
          }
)
