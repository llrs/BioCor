#' @importClassesFrom GSEABase GeneSetCollection
NULL

#' @exportMethod summary
setMethod("summary",
          "GeneSetCollection",
          function(object) {
              check_gsc(object)
              gpp <- GenesPerPathway(object)
              ppg <- PathwaysPerGene(object)

              h_gpp <- h_index(gpp)
              h_ppg <- h_index(ppg)

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
              cat("\th-index:", h_gpp, "genes with at least", h_gpp, "pathways.\n")
              cat("Pathways:", out["pathways"], "\n")
              cat("\tBiggest pathway:", out["maxGPP"], "genes\n")
              cat("\th-index:", h_ppg, "pathways with at least", h_ppg, "genes.\n")
              cat("IC(GenesPerPathway):", round(out["ICgpp"], 2), "(", round(out["ICgpp"]/out["maxICgpp"], 2), "% of the maximum)\n")
              cat("IC(PathwaysPerGene)", round(out["ICppg"], 2), "(", round(out["ICppg"]/out["maxICppg"], 2), "% of the maximum)\n")

            invisible(out)
          }
)
