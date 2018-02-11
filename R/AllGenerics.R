#' @export pathSim
setGeneric("pathSim", function(pathway1, pathway2, info, ...)
    standardGeneric("pathSim"))

#' @export mpathSim
setGeneric("mpathSim", function(pathways, info, method, ...)
    standardGeneric("mpathSim"))

#' @export geneSim
setGeneric("geneSim", function(gene1, gene2, info, method, ...)
    standardGeneric("geneSim"))

#' @export mgeneSim
setGeneric("mgeneSim", function(genes, info, method, ....)
    standardGeneric("mgeneSim"))

#' @export clusterSim
setGeneric("clusterSim", function(cluster1, cluster2, info, method, ...)
    standardGeneric("clusterSim"))

#' @export mclusterSim
setGeneric("mclusterSim", function(clusters, info, method, ....)
    standardGeneric("mclusterSim"))

#' @export clusterGeneSim
setGeneric("clusterGeneSim", function(cluster1, cluster2, info, method, ...)
    standardGeneric("clusterGeneSim"))

#' @export mclusterGeneSim
setGeneric("mclusterGeneSim", function(clusters, info, method, ....)
    standardGeneric("mclusterGeneSim"))
