#' @importClassesFrom GSEABase GeneSetCollection
NULL

setGeneric("ratio", function(object) standardGeneric("ratio"))


setMethod("ratio",
          signature(object = "GeneSetCollection"),
          function(object) {
              check_gsc(object)

              gpp <- GenesPerPathway(object)
              ppg <- PathwaysPerGene(object)

              uPPG <- unique(ppg)
              uGPP <- unique(gpp)
              # Create the matrix to fill
              m <- matrix(0, ncol = length(uPPG), nrow = length(uGPP),
                          dimnames =
                              list("GenesPerPathway" = uGPP[order(uGPP)],
                                          "PathwaysPerGene" = uPPG[order(uPPG)]))

              paths2genes <- GSEABase::geneIds(object)
              genes2paths <- inverseList(paths2genes)

              # Find the genes with those pathways
              for (g in rownames(m)) {
                  p <- names(gpp)[gpp == g]
                  genes <- unique(unlist(paths2genes[p], use.names = FALSE))
                  subPPG <- ppg[genes]
                  # Find the pathways with those number of genes
                  for (pSize in colnames(m)) {
                      dupli <- names(subPPG)[subPPG == pSize]
                      m[g, pSize] <- length(dupli)
                  }
              }
              m
          }
)




setGeneric("condGene", function(object) standardGeneric("condGene"))

#' @export
setMethod("condGene",
          signature(object = "GeneSetCollection"),
          function(object) {
              check_gsc(object)

              gpp <- GenesPerPathway(object)
              ppg <- PathwaysPerGene(object)

              uPPG <- unique(ppg)
              uGPP <- unique(gpp)
              # Create the matrix to fill
              m <- matrix(0, ncol = length(uPPG), nrow = length(uGPP),
                          dimnames =
                              list("GenesPerPathway" = uGPP[order(uGPP)],
                                   "PathwaysPerGene" = uPPG[order(uPPG)]))

              paths2genes <- GSEABase::geneIds(object)
              genes2paths <- inverseList(paths2genes)

              # Find the genes with those pathways
              for (g in colnames(m)) {
                  genes <- names(ppg)[ppg == g]
                  paths <- unique(unlist(genes2paths[genes], use.names = FALSE))
                  subGPP <- gpp[paths]
                  # Find the pathways with those number of genes
                  for (pSize in rownames(m)) {
                      if (length(pSize) == 0) {
                          m[pSize, g] <- 0
                      }
                      m[pSize, g] <- sum(subGPP == pSize, na.rm = TRUE)/length(subGPP)
                  }
              }
              m
          }
)

#' @export
setGeneric("condPath", function(object) standardGeneric("condPath"))

#' @export
setMethod("condPath",
          signature(object = "GeneSetCollection"),
          function(object) {
              check_gsc(object)

              gpp <- GenesPerPathway(object)
              ppg <- PathwaysPerGene(object)

              uPPG <- unique(ppg)
              uGPP <- unique(gpp)
              # Create the matrix to fill
              m <- matrix(0, ncol = length(uPPG), nrow = length(uGPP),
                          dimnames =
                              list("GenesPerPathway" = uGPP[order(uGPP)],
                                   "PathwaysPerGene" = uPPG[order(uPPG)]))

              paths2genes <- GSEABase::geneIds(object)
              genes2paths <- inverseList(paths2genes)

              # Find the genes with those pathways
              for (p in rownames(m)) {
                  paths <- names(gpp)[gpp == p]
                  genes <- unlist(paths2genes[paths], use.names = FALSE)
                  subPPG <- ppg[genes]
                  # Find the pathways with those number of genes
                  for (gSize in colnames(m)) {
                      if (length(gSize) == 0) {
                          m[p, gSize] <- 0
                      }
                      m[p, gSize] <- sum(subPPG == gSize, na.rm = TRUE)/length(subPPG)
                  }
              }
              m
          }
)

#' @export
setGeneric("GenesSizePathways", function(object) standardGeneric("GenesSizePathways"))



#' @export
setMethod("GenesSizePathways",
          signature(object = "GeneSetCollection"),
          function(object) {
              check_gsc(object)

              genes2paths <- inverseList(GSEABase::geneIds(object))
              gpp <- GenesPerPathway(object)
              u_gpp <- unique(gpp)
              u_gpp <- u_gpp[order(u_gpp)]
              dimn <- list("SizePaths" = u_gpp, "Genes" = names(genes2paths))
              out <- matrix(0, ncol = length(genes2paths),
                            nrow = length(unique(gpp)),
                            dimnames = dimn)
              for (g in names(genes2paths)) {
                  calc <- table(gpp[genes2paths[[g]]])
                  out[names(calc), g] <- calc
              }
              out
          }
)
