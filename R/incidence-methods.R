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
PathwayPerGene <- function(incidence){
    colSums(incidence)
}

#' @exportMethod PathwayPerGene
setGeneric("PathwayPerGene", function(object) {
    standardGeneric("PathwayPerGene")
})

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


#' @exportMethod PathwayPerGene
setMethod("PathwayPerGene",
          "GeneSetCollection",
          function(object) {
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
              lengths(GSEABase::geneIds(object))
          }
)
