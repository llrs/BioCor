#' @importClassesFrom GSEABase GeneSetCollection
NULL

# Table of Length pathways vs Number of pathways a gene is in
#' Distribution of genes per pathway in pathways
#'
#' Calculates how many pathways per gene are distributed in genes per pathways.
#' @param object A GeneSetCollection.
#' @return A matrix with the number of pathways per genes as columns and the
#' genes per pathways as rows.
#' @author Lluís Revilla
#' @exportMethod distribution
setGeneric("distribution", function(object) {
    standardGeneric("distribution")
})

#' @exportMethod distribution
setMethod("distribution",
          signature(object = "GeneSetCollection"),
          function(object) {
              check_gsc(object)

              gpp <- GenesPerPathway(object)
              ppg <- PathwaysPerGene(object)

              m <- matrix(0, ncol = length(table(ppg)),
                          nrow = length(table(gpp)),
                          dimnames = list("GenesPerPathway" = names(table(gpp)),
                                          "PathwaysPerGene" = names(table(ppg))))

              paths2genes <- GSEABase::geneIds(object)
              genes2paths <- inverseList(paths2genes)

              for (pp in colnames(m)) {
                  genes <- names(genes2paths)[lengths(genes2paths) == as.numeric(pp)]
                  paths <- unique(unlist(genes2paths[genes], use.names = FALSE))
                  tab <- as.matrix(table(lengths(paths2genes[paths])))
                  colnames(tab) <- pp
                  m <- AintoB(tab, m)

              }
              m
          }
)

# # Proportion relative to pathways per gene
# rel <- sweep(m, 2, colSums(m), "/")
# rel <- reshape2::melt(rel)
# library("ggplot2")
# library("reshape2")
#
# rel$GenesPerPathway <- as.factor(rel$GenesPerPathway)
# rel$PathwaysPerGene <- as.factor(rel$PathwaysPerGene)
# ggplot(rel, aes(x=value, fill=PathwaysPerGene)) +
#     geom_density(alpha=0.25) + theme_bw() + guides(fill = FALSE)
#
# ggplot(rel, aes(x=value, fill=GenesPerPathway)) +
#     geom_density(alpha=0.25) + theme_bw() + guides(fill = FALSE)
#
# # Proportion relative to pathways per gene
# rel <- sweep(m, 2, sum(colSums(m)), "/")
# rel <- reshape2::melt(rel)
# library("ggplot2")
# library("reshape2")
#
# rel$GenesPerPathway <- as.factor(rel$GenesPerPathway)
# rel$PathwaysPerGene <- as.factor(rel$PathwaysPerGene)
# ggplot(rel, aes(x=value, fill=PathwaysPerGene)) +
#     geom_density(alpha=0.25) + theme_bw() + guides(fill = FALSE)
#
# ggplot(rel, aes(x=value, fill=GenesPerPathway)) +
#     geom_density(alpha=0.25) + theme_bw() + guides(fill = FALSE)
