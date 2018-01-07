#' @importClassesFrom GSEABase GeneSetCollection
NULL

#' Statistics about pathways
#'
#' Calculates some statistics of the genes, like the biggest
#' pathway, the gene in more pathways, h-index and entropy or information
#' content for the genes per pathway and pathways for genes and they percentage
#' of the maximum.
#' @param object A GeneSetCollection object
#' @param pathway A character of the pathway to analyze. If missing it will analyze
#' all pathways
#' @return If a single pathway is provided invisible a list of statistics (the
#' statistics will be printed on the screen), otherwise the data.frame of
#' statistics.
#' @author Lluís Revilla
#' @exportMethod pathway
setGeneric("pathway", function(object, pathway) {
    standardGeneric("pathway")
})

#' @export
setGeneric("pathway",
           function(object, pathway) {
               stop("Incorrect input object should be a GeneSetCollection and
                    pathway a character")
           }
)

#' @describeIn pathway Calculates statistics for a single pathway
#' @exportMethod pathway
setMethod("pathway",
          signature(object = "GeneSetCollection", pathway = "character"),
          function(object, pathway) {
              check_gsc(object)
              ppg <- PathwaysPerGene(object)
              gpp <- GenesPerPathway(object)

              paths2genes <- GSEABase::geneIds(object)
              genes2paths <- inverseList(paths2genes)

              totalGenes <- length(genes2paths)
              totalPathways <- length(paths2genes)

              Tgpp <- prop.table(table(gpp))
              Tppg <- prop.table(table(ppg))


              if (!pathway %in% names(paths2genes)) {
                  stop("Pathway provided is not in the GeneSetCollection.")
              }

              genes <- unique(unlist(paths2genes[pathway], use.names = FALSE))
              Genes_pathways <- lengths(genes2paths[genes])

              ic <- IC(Genes_pathways)
              maxIC <- maxIC(ppg)

              out <- c(ic, # Information content
                       maxIC(Genes_pathways), # Relative to its possibilities
                       length(genes), # Number of pathways per gene
                       Tgpp[as.character(length(genes))], # Probability to find a gene with those pathways
                       prod(Tgpp[Genes_pathways], na.rm = TRUE)) # Probability to find a gene with pathways of these length
              names(out) <- c("IC", "maxIC", "genes", "Prob_#_genes", "Prob_pathwas_#_genes")

              # cat("Genes ", out["genes"], "\n")
              # cat("Pathways ", out["pathways"], "\n")
              # cat("IC(GenesPerPathway) ", round(out["ICgpp"], 2), "\n")
              cat("IC(PathwaysPerGene) ", round(out["ICppg"], 2), "\n")
              cat("Pathways", out["pathways"], "\n")
              cat("Probability", out["Prob"], "\n")

              invisible(out)
          }
)

#' @describeIn pathway Calculates statistics for all pathways
#' @exportMethod pathway
setMethod("pathway",
          signature(object = "GeneSetCollection", pathway = "missing"),
          function(object, pathway) {
              ppg <- PathwaysPerGene(object)
              gpp <- GenesPerPathway(object)

              paths2genes <- GSEABase::geneIds(object)
              genes2paths <- inverseList(paths2genes)

              totalGenes <- length(genes2paths)
              totalPathways <- length(paths2genes)

              Tgpp <- prop.table(table(gpp))
              Tppg <- prop.table(table(ppg))

              maxIC <- maxIC(ppg)
              out <- sapply(names(paths2genes), function(pathway) {
                  genes <- unique(unlist(paths2genes[pathway], use.names = FALSE))
                  Genes_length <- lengths(genes2paths[genes])

                  if (length(genes) == 1) {
                      ic <- 0
                  } else {
                      ic <- IC(Genes_length)
                  }
                  out <- c(ic, # Information content
                           maxIC(Genes_length), # Relative to its possibilities
                           length(genes), # Number of pathways per gene
                           Tgpp[as.character(length(genes))], # Probability to find a gene with those pathways
                           prod(Tppg[as.character(Genes_length)], na.rm = TRUE)) # Probability to find a pathway with genes at n-pathways these
                  names(out) <- c("IC", "maxIC", "genes", "Prob_#_genes", "Prob_pathways_#_genes")
                  out
              })
              out <- t(out)
              as.data.frame(out)
          }
)


# IC_{ppg}/IC_{ppg_{total}} vs number of Pathways shows that there is more
#  information for pathways bigger than 25 and the assimptodata is close to 1.5
# ggplot(o, aes(pathways, ICppg/4.49)) + geom_count()

