#' @importClassesFrom GSEABase GeneSetCollection
NULL

#' Statistics about genes
#'
#' Calculates some statistics of the genes, like the biggest
#' pathway, the gene in more pathways, h-index and entropy or information
#' content for the genes per pathway and pathways for genes and they percentage
#' of the maximum.
#' @param object A GeneSetCollection object
#' @param gene A character of the gene to analyze. If missing it will analyze
#' all genes
#' @return If a single gene is provided invisible a list of statistics (the
#' statistics will be printed on the screen), otherwise the data.frame of
#' statistics.
#' @author Lluís Revilla
#' @exportMethod gene
setGeneric("gene", function(object, gene) standardGeneric("gene") )

#' @describeIn gene Calculates statistics for a single gene
#' @exportMethod gene
setMethod("gene",
          signature(object = "GeneSetCollection", gene = "character"),
          function(object, gene) {
              check_gsc(object)
              ppg <- PathwaysPerGene(object)
              gpp <- GenesPerPathway(object)

              paths2genes <- GSEABase::geneIds(object)
              genes2paths <- inverseList(paths2genes)

              totalGenes <- length(genes2paths)
              totalPathways <- length(paths2genes)

              Tgpp <- prop.table(table(gpp))
              Tppg <- prop.table(table(ppg))


              if (! gene %in% names(genes2paths)) {
                  stop("Gene provided is not in the GeneSetCollection.")
              }

              pathways <- unique(unlist(genes2paths[gene], use.names = FALSE))
              Pathways_length <- lengths(paths2genes[pathways])

              ic <- IC(Pathways_length)
              maxIC <- maxIC(ppg)

              out <- c(ic, # Information content
                       maxIC(Pathways_length), # Relative to its possibilities
                       length(pathways), # Number of pathways per gene
                       Tppg[as.character(length(pathways))], # Probability to find a gene with those pathways
                       prod(Tgpp[as.character(Pathways_length)], na.rm = TRUE)) # Probability to find a gene with pathways of these length
              names(out) <- c("IC", "maxIC", "pathways", "Prob_#_pathways", "Prob_genes_#_pathways")

                            # cat("Genes ", out["genes"], "\n")
              # cat("Pathways ", out["pathways"], "\n")
              # cat("IC(GenesPerPathway) ", round(out["ICgpp"], 2), "\n")
              cat("IC(PathwaysPerGene) ", round(out["ICppg"], 2), "\n")
              cat("Pathways", out["pathways"], "\n")
              cat("Probability", out["Prob"], "\n")

              invisible(out)
          }
)

#' @describeIn gene Calculates statistics for all genes
#' @exportMethod gene
setMethod("gene",
          signature(object = "GeneSetCollection", gene = "missing"),
          function(object, gene) {
              ppg <- PathwaysPerGene(object)
              gpp <- GenesPerPathway(object)

              paths2genes <- GSEABase::geneIds(object)
              genes2paths <- inverseList(paths2genes)

              totalGenes <- length(genes2paths)
              totalPathways <- length(paths2genes)

              Tgpp <- prop.table(table(gpp))
              Tppg <- prop.table(table(ppg))

              maxIC <- maxIC(ppg)
              out <- sapply(names(genes2paths), function(gene) {
                  pathways <- unique(unlist(genes2paths[gene], use.names = FALSE))
                  Pathways_length <- lengths(paths2genes[pathways])

                  if (length(pathways) == 1) {
                      ic <- 0
                  } else {
                      ic <- IC(Pathways_length)
                  }
                  out <- c(ic, # Information content
                           maxIC(Pathways_length), # Relative to its possibilities
                           length(pathways), # Number of pathways per gene
                           Tppg[as.character(length(pathways))], # Probability to find a gene with those pathways
                           prod(Tgpp[as.character(Pathways_length)], na.rm = TRUE)) # Probability to find a gene with pathways of these length
                  names(out) <- c("IC", "maxIC", "pathways", "Prob_#_pathways", "Prob_genes_#_pathways")
                  out
              })
              out <- t(out)
              as.data.frame(out)
          }
)


# IC_{ppg}/IC_{ppg_{total}} vs number of Pathways shows that there is more
#  information for pathways bigger than 25 and the assimptodata is close to 1.5
# ggplot(o, aes(pathways, ICppg/4.49)) + geom_count()

