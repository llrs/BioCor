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
              ppg <- PathwaysPerGene(object)

              tab_gpp <- table(gpp)
              tab_ppg <- table(ppg)

              h_gpp <- tab_gpp[as.numeric(names(tab_gpp)) == tab_gpp]
              h_ppg <- tab_ppg[as.numeric(names(tab_ppg)) == tab_ppg]

              out <- c(ICgpp = IC(gpp),
                       ICppg = IC(ppg),
                       maxICgpp = maxIC(gpp),
                       maxICppg = maxIC(ppg),
                       genes = length(ppg),
                       maxGPP = max(gpp),
                       maxPPG = max(ppg),
                       pathways = length(gpp)
              )
              cat("Genes:", out["genes"], "\n")
              cat("\tGene in more pathways:", out["maxPPG"], "pathways\n")
              cat("\th-index:", h_ppg, "genes with", h_ppg, "pathways\n")
              cat("Pathways:", out["pathways"], "\n")
              cat("\tBiggest pathway:", out["maxGPP"], "genes\n")
              cat("\th-index:", h_gpp, "pathways with", h_gpp, "genes.\n")
              cat("IC(GenesPerPathway):", round(out["ICgpp"], 2), "close to the maximum by", round(out["ICgpp"]/out["maxICgpp"], 2), "\n")
              cat("IC(PathwaysPerGene)", round(out["ICppg"], 2), "close to the maximum by", round(out["ICppg"]/out["maxICppg"], 2), "\n")

            invisible(out)
          }
)
