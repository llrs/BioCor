#' @importClassesFrom GSEABase GeneSetCollection
NULL


#' @exportMethod summary
setGeneric("summary", function(object) {
    standardGeneric("summary")
})

#' @exportMethod summary
setMethod("summary",
          "GeneSetCollection",
          function(object) {
              gpp <- GenesPerPathway(object)
              ppg <- PathwayPerGene(object)
              out <- c(ICgpp = IC(gpp),
                       ICppg = IC(ppg),
                       genes = length(ppg),
                       maxGPP = max(gpp),
                       maxPPG = max(ppg),
                       pathways = length(gpp)
              )
              cat("Genes ", out["genes"], "\n\tGene in more pathways:", out["maxPPG"], "pathways\n")
              cat("Pathways ", out["pathways"], "\n\tBiggest pathway:", out["maxGPP"], "genes\n")
              cat("IC(GenesPerPathway) ", round(out["ICgpp"], 2), "\n")
              cat("IC(PathwaysPerGene) ", round(out["ICppg"], 2))

            invisible(out)
          }
)
