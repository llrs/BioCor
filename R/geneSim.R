# geneSim ####
#' Calculates the Dice similarity score of two genes
#'
#' For the genes given calculates the Dice similarity between each pathway
#' which is combined to obtain a similarity between the genes.
#'
#' Given the information about the genes and their pathways, uses the ids
#' of the genes to find the Dice similarity score for each pathway comparison
#' between the genes. Later this similarities are combined using
#' \code{\link{combineScores}}.
#' @param gene1,gene2 Ids of the genes to calculate the similarity, to be found
#' in genes.
#' @param genes is the matrix with the information to calculate the similarity.
#' It is created using select, should contain the column "id" and "pathwayDB".
#' @param id is the column of \code{genes} where \code{gene1} and \code{gene2}
#' are to be found
#' @param pathwayDB is the column where pathways should be found. It is usually
#' the name of the database where they come from.
#' @param method To combine the scores of each pathway, one of \code{c("avg",
#' "max", "rcmax", "rcmax.avg", "BMA")}, if NULL returns the matrix.
#' @return
#' The highest Dice score of all the combinations of pathways between
#' the two ids compared if a method to combine scores is provided or NA if
#' there isn't information for one gene.
#' If an \code{NA} is returned this means that there isn't information
#' available for any pathways for one of the genes. Otherwise a number
#' between 0 and 1 (both included) is returned. Note that there isn't a
#' negative value of similarity.
#' @export
#' @author Lluis Revilla
#' @seealso \code{\link{conversions}} help page to transform Dice
#' score to Jaccard score. For the method to combine the scores see
#' \code{\link{combineScores}}.
#' @examples
#' library("org.Hs.eg.db")
#' library("reactome.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- select(org.Hs.eg.db, keys = entrezids, keytype = "ENTREZID",
#'                      columns = "PATH")
#' # Extracts the paths of all genes of org.Hs.eg.db from reactome
#' genes.react <- select(reactome.db, keys = entrezids, keytype = "ENTREZID",
#'                       columns = "REACTOMEID")
#' geneSim("81", "18", genes.react, "ENTREZID", "REACTOMEID")
#' geneSim("81", "18", genes.kegg, "ENTREZID", "PATH")
#' geneSim("81", "18", genes.react, "ENTREZID", "REACTOMEID", NULL)
#' geneSim("81", "18", genes.kegg, "ENTREZID", "PATH", NULL)
geneSim <- function(gene1, gene2, genes, id, pathwayDB, method = "max") {
    comb <- c(gene1, gene2)
    if (!pathwayDB %in% colnames(genes)) {
        stop("Please check which type of pathway do you want")
    }
    if (any(is.na(comb))) {
        return(NA)
    }
    if (comb[1L] == comb[2L]) {
        return(1)
    }
    if (length(comb) > 2L) {
        stop("comb can only be of length 2, to compare pairs of genes")
    }

    # Convert data.frame into environment to speed the look up
    genes2pathways <- split(genes[ , pathwayDB], genes[, id])
    genes2pathways <- list2env(genes2pathways)

    pathways2genes <- split(genes[ , id], genes[, pathwayDB])
    pathways2genes <- list2env(pathways2genes)

    # Extract all pathways for each gene
    pathways <- sapply(comb, function(x) {
        genes2pathways[[x]]
    }, simplify = FALSE)

    # Check that we have pathways info for this combination
    if (any(lengths(pathways) == 0L)) {
        return(NA)
    }
    # Extract the gene ids for each pathway
    g1 <- sapply(pathways[[1]], function(x) {
        pathways2genes[[x]]
    }, simplify = FALSE)
    g2 <- sapply(pathways[[2]], function(x) {
        pathways2genes[[x]]
    }, simplify = FALSE)

    # Calculate the dice index
    vdiceSim <- Vectorize(diceSim)
    react <- outer(g1, g2, vdiceSim)

    # Calculate the similarity between the two genes
    if (is.null(method)) {
        react
    } else {
        combineScores(react, method)
    }
}

#' @rdname geneSim
#' @export
#' @param gene.list Given a list of vectors return the similarities
#' @return \code{mgeneSim} returns the matrix of similarities between the genes
#' in the vector
#' @examples
#'
#' mgeneSim(c("81", "18", "10"), genes.react, "ENTREZID", "REACTOMEID")
#' mgeneSim(c("81", "18", "10"), genes.react, "ENTREZID", "REACTOMEID", "avg")
mgeneSim <- function(gene.list, genes, id, pathwayDB, method = "max") {
    vg <- Vectorize(geneSim, vectorize.args = c("gene1", "gene2"))
    names(gene.list) <- gene.list
    outer(gene.list, gene.list, vg, genes = genes, id = id,
          pathwayDB = pathwayDB, method = method)
}
