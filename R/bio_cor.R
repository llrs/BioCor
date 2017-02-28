# Functions required for the bio.cor or bio.cor2 computation.


# diceSim ####
#' Compare pathways
#'
#' Function to estimate how much two graphs or list of genes overlap by looking
#' how much of the nodes are shared.
#' @param g1,g2 Graph in GraphNEL format, or a character list with the names of
#' the proteins in each pathway.
#' @return A score between 0 and 1 calculated as the doble of the proteins
#' shared by g1 and g2 divided by the number of genes in both groups.
#' @importFrom methods is
#' @importFrom graph nodes
#' @export
#' @author Lluís Revilla
#' @seealso Used for \code{\link{genesSim}}, see \code{\link{conversions}} help
#' page to transform Dice score to Jaccard score.
#' @examples
#' genes.id2 <- c("52", "11342", "80895", "57654", "548953", "11586", "45985")
#' genes.id1 <- c("52", "11342", "80895", "57654", "58493", "1164", "1163",
#' "4150", "2130", "159")
#' diceSim(genes.id1, genes.id2)
#' diceSim(genes.id2, genes.id2)
diceSim <- function(g1, g2) {
    # Check which case are we using
    if (is(g1, "graph") & is(g2, "graph")) {
        prot1 <- nodes(g1)
        prot2 <- nodes(g2)
    } else if (is(g1, "graph") & is.character(g2)) {
        prot1 <- nodes(g1)
        prot2 <- g2
    } else if (is(g2, "graph") & is.character(g1)) {
        prot2 <- nodes(g2)
        prot1 <- g1
    } else {
        prot1 <- g1
        prot2 <- g2
    }
    # If there isn't any information of a pathway for a gene then then
    # functional similarity is 0
    if (length(prot1) == 0L | length(prot2) == 0L) {
        return(0L)
    }
    score <- (length(intersect(prot1, prot2)))*2L/(
        length(prot2) + length(prot1))
    score
}

# pathSim ####
#' Compare pathways
#'
#' Given two vectors of pathways calculates the similarity between them.
#'
#' \code{diceSim} is used to calculate similarities between each pathway and
#' combineScores to extract the similarity between those pathways. If one need
#' the matrix of similarities set methods to \code{NULL}.
#' @param pathways1,pathways2 Pathways to be found in
#' \code{pathwayDB}.
#' @inheritParams genesSim
#' @return The similarity between those pathways or all the similarities
#' between each comparison.
#' @seealso \code{\link{pathSim}} and \code{\link{combineScores}}
#' @author Lluís Revilla
#' @export
#' @examples
#' library("org.Hs.eg.db")
#' library("reactome.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' # Extracts the paths of all genes of org.Hs.eg.db from reactome
#' genes.react <- select(reactome.db, keys = entrezids, keytype = "ENTREZID",
#'                       columns = "REACTOMEID")
#' pathways1 <- c("112310", "112315", "112316", "373753", "916853")
#' pathways2 <- c("109582", "114608", "1500931", "888590", "76002", "76005")
#' pathSim(pathways1[1], pathways2[1], genes.react, "ENTREZID", "REACTOMEID")
#' mpathSim(pathways1, pathways2, genes.react, "ENTREZID", "REACTOMEID", NULL)
pathSim <- function(pathways1, pathways2, genes, id, pathwayDB) {

    # Extract the gene ids for each pathway
    g1 <- genesInfo(genes, pathwayDB, pathways1, id)[[1]]
    g2 <- genesInfo(genes, pathwayDB, pathways2, id)[[1]]

    diceSim(g1, g2)
}

#' @rdname pathSim
#' @export
mpathSim <- function(pathways1, pathways2, genes, id, pathwayDB,
                     method = "max") {

    # Extract the gene ids for each pathway
    g1 <- genesInfo(genes, pathwayDB, pathways1, id)
    g2 <- genesInfo(genes, pathwayDB, pathways2, id)

    vp <- Vectorize(diceSim)
    react <- outer(g1, g2, vp)

    # Calculate the similarity between the two genes
    if (is.null(method)) {
        react
    } else {
        combineScores(react, method)
    }
}

# removeDup ####
#' Remove duplicated rows and columns
#'
#' Given the indices of the duplicated entries remove the columns and rows
#' until just one is left, it keeps the duplicated with the highest absolute
#' mean value.
#'
#' @param cor_mat List of matrices
#' @param dupli List of indicies with duplicated entries
#' @return A matrix with only one of the columns and rows duplicated
#' @export
#' @author Lluís Revilla
#' @seealso \code{\link{duplicateIndices}} to obtain the list of indicies with
#' duplicated entries.
#' @examples
#' a <- seq2mat(c("52", "52", "53", "55"), runif(choose(4, 2)))
#' b <- seq2mat(c("52", "52", "53", "55"), runif(choose(4, 2)))
#' mat <- list("kegg" = a, "react" = b)
#' mat
#' dupli <- duplicateIndices(rownames(a))
#' remat <- removeDup(mat, dupli)
#' remat
removeDup <- function(cor_mat, dupli) {
    if (!all(sapply(cor_mat, isSymmetric))) {
        stop("All the matrices of mat should be symmetric and with the same ",
             "column names and rownames")
    }
    cor_mat <- Map(function(mat, x = dupli) {
        rem.colum <- sapply(x, function(y, m) {
            mean.column <- apply(m[, y], 2L, mean, na.rm = TRUE)
            i <- which.max(abs(mean.column))
            # Select those who don't bring more information
            rem.colum <- setdiff(y, y[i])
        }, m = mat)

        mat[-rem.colum, -rem.colum]
    }, cor_mat)
    return(cor_mat)
}


# bioCor ####
#' Comparing information in databases
#'
#' Calculates a functional similarity of genes  using information available in
#' several databases.
#'
#' For Gene Ontologies, the DAG path structure is used to compute how similar
#' two genes are. For metabolic pathways the max number of proteins involved in
#' a pathway for each gene is calculated.
#'
#' @param genes_id is vector of ids to compare
#' @param ids indicate if the id is eihter ENTREZID or SYMBOL
#' @param react logical; indicates if the similarities in Reactome pathway is
#' calculated
#' @param kegg logical; indicates if the similarities in Kegg database
#' is calculated
#' @param all logical; indicates if all the previous (go, react, and kegg)
#' similarity measures, should be set to TRUE.
#' @inheritParams BiocParallel::bpvec
#' @return A list where each element is a matrix with the similarity for such
#' database \code{NA} indicates that there isn't any information of one of
#' those genes for that pathway database. 0 If there is information of both
#' genes but no overlap between those genes.
#' @importFrom reactome.db reactome.db
#' @importFrom AnnotationDbi select
#' @importFrom AnnotationDbi keys
# #' @importFrom biganalytics apply
#' @import BiocParallel
#' @import org.Hs.eg.db
# #' @import foreach
# #' @import bigmemory
#' @export
#' @author Lluís Revilla
#' @seealso \code{\link{genesSim}}
#' @examples
#' bioCor(c("18", "81"))
bioCor <- function(genes_id, ids = "ENTREZID", react = TRUE, kegg = FALSE,
                   all = FALSE, BPPARAM = bpparam()) {
    if (!ids %in% c("ENTREZID", "SYMBOL")) {
        stop("Please check the input of genes in ENTREZID or SYMBOL format")
    }
    if (all) {
        kegg <- react <- all
    }

    # Obtain data from the annotation packages
    if (ids == "SYMBOL") {
        gene.symbol <- suppressMessages(select(org.Hs.eg.db, keys = genes_id,
                              keytype = "SYMBOL", columns = "ENTREZID"))
    } else {
        gene.symbol <- suppressMessages(select(org.Hs.eg.db, keys = genes_id,
                              keytype = "ENTREZID", columns = "SYMBOL"))
    }

    n.combin <- choose(length(gene.symbol$ENTREZID), 2L)
    orig.ids <- gene.symbol$ENTREZID

    if (sum(is.na(gene.symbol$ENTREZID)) >= 1L) {
        message("Some symbols are not mapped to Entrez Genes IDs")
    }
    dup_symb <- duplicated(gene.symbol$SYMBOL[
        !is.na(gene.symbol$ENTREZID)])

    if (sum(dup_symb) >= 1L & ids == "SYMBOL") {
        message("Some symbols are mapped to several Entrez Genes IDs.")
    } else if (sum(dup_symb) >= 1L & ids == "Entrez Gene") {
        message("Some Entrez Genes IDs are mapped to several symbols.")
    }
    entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")
    # Obtain the data of kegg and Reactome pathways
    if (kegg) {
        # Obtain data
        gene.kegg <- suppressMessages(select(org.Hs.eg.db,
                                             keys = entrez,
                            keytype = "ENTREZID", columns = "PATH"))
        # Merge data
        genes <- unique(merge(gene.symbol, gene.kegg, all = TRUE,
                              sort = FALSE))
    }

    if (react) {
        if (!kegg) {
            genes <- gene.symbol
        }
        if (all(!gene.symbol$ENTREZID %in% keys(reactome.db))) {
            gene.reactome <- cbind(ENTREZID = gene.symbol$ENTREZID,
                                   REACTOMEID = NA)
        } else {
            gene.reactome <- suppressMessages(select(reactome.db,
                                    keys = entrez,
                                    keytype = "ENTREZID",
                                    columns = "REACTOMEID"))
        }

        genes <- unique(merge(genes, gene.reactome, all = TRUE, sort = FALSE))
    }

    if (kegg) {  # parallel # to run non parallel transform the %dopar% into
        # %do%
        # kegg_mat <- big.matrix(length(orig.ids), length(orig.ids), init = NA,
        #                        dimnames = list(orig.ids, orig.ids),
        #                        backingfile = "kegg.bin",
        #                        descriptorfile = "kegg.desc")
        # datadesc <- describe(kegg_mat)

        message("Calculating KEGG similarities")
        # kegg.bio <- foreach(i = orig.ids, .combine = c,
        #                     .verbose = FALSE) %dopar% {
        #                         for (j in orig.ids) {
        #                             kegg_mat <- attach.big.matrix(datadesc)
        #                             kegg_mat[i, j] <- genesSim(
        #                                 c(i, j), genes,"ENTREZID", "PATH")
        #                         }
        #                     }
        kegg.bio <- bpmapply(function(x){
            comb <- combinadic(n = orig.ids, r = 2, i = x)
            genesSim(comb[1], comb[2], genes, "ENTREZID", "PATH")},
            seq_len(n.combin), BPPARAM = BPPARAM)
        message("KEGG similarities has been calculated")

        kegg_mat <- seq2mat(orig.ids, kegg.bio)
        if (all(apply(kegg_mat, 1L, function(x){
            sum(is.na(x)) == length(x) - 1L }))) {
            warning("KEGG didn't found relevant similarities!")
        }
    }

    if (react) {  # parallel # to run non parallel transform the %dopar% into
        message("Calculating REACTOME similarities")
        react.bio <- bpmapply( function(x){
            comb <- combinadic(n = orig.ids, r = 2, i = x)
            genesSim(comb[1], comb[2], genes, "ENTREZID", "REACTOMEID")},
            seq_len(n.combin), BPPARAM = BPPARAM)

        message("REACTOME similarities has been calculated")
        react_mat <- seq2mat(orig.ids, react.bio)

        if (all(apply(react_mat, 1L, function(x){
            sum(is.na(x)) == length(x) - 1L}))) {
            warning("REACTOME didn't found relevant similarities!")
        }
    }

    if (kegg & react) {
        cor_mat <- list(reactome = as.matrix(react_mat),
                        kegg = as.matrix(kegg_mat))
    } else if (react) {
        cor_mat <- list(react = as.matrix(react_mat))
    } else if (kegg) {
        cor_mat <- list(kegg = as.matrix(kegg_mat))
    }

    # To match the ncol and nrow of the input with the initial input
    dupli <- duplicateIndices(gene.symbol[, ids])

    # Keep the interesting columns
    if (length(dupli) >= 1L) {
        cor_mat <- removeDup(cor_mat, dupli)
    }
    # Rename to the original ids
    cor_mat <- lapply(cor_mat, function(x) {
        colnames(x) <- genes_id
        rownames(x) <- genes_id
        x})

    return(cor_mat)
}

# genesSim ####
#' Calculates the Dice similarity score of two genes
#'
#' For the genes given calculates the Dice similarity between each pathway
#' which is combined to obtain a similarity between the genes.
#'
#' Given the information about the genes and their pathways, uses the ids
#' of the genes to find the Dice similarity score for each pathway comparison
#' between the genes. Later this similarities are combined using
#' \code{\link{combineScores}}.
#' @param gene1,gene2 Entrez gene id.
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
#' genesSim("81", "18", genes.react, "ENTREZID", "REACTOMEID")
#' genesSim("81", "18", genes.kegg, "ENTREZID", "PATH")
#' genesSim("81", "18", genes.react, "ENTREZID", "REACTOMEID", NULL)
genesSim <- function(gene1, gene2, genes, id, pathwayDB, method = "max") {
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

    # Extract all pathways for each gene
    pathways <- genesInfo(genes, id, comb, pathwayDB)
    # Check that we have pathways info for this combination
    if (any(lengths(pathways) == 0L)) {
        return(NA)
    }
    # Extract the gene ids for each pathway
    g1 <- genesInfo(genes, pathwayDB, pathways[[1]], id)
    g2 <- genesInfo(genes, pathwayDB, pathways[[2]], id)

    vdiceSim <- Vectorize(diceSim)

    react <- outer(g1, g2, vdiceSim)

    # Calculate the similarity between the two genes
    if (is.null(method)) {
        react
    } else {
        combineScores(react, method)
    }
}
#' @rdname genesSim
#' @export
#' @param gene.list Given a list of vectors return the similarities
#' @return \code{mgeneSim} returns the matrix of similarities between the genes
#' in the vector
#' @examples
#'
#' mgeneSim(c("81", "18", "10"), genes.react, "ENTREZID", "REACTOMEID")
#' mgeneSim(c("81", "18", "10"), genes.react, "ENTREZID", "REACTOMEID", "avg")
mgeneSim <- function(gene.list, genes, id, pathwayDB, method = "max") {
    vg <- Vectorize(genesSim, vectorize.args = c("gene1", "gene2"))
    names(gene.list) <- gene.list
    outer(gene.list, gene.list, vg, genes = genes, id = id,
          pathwayDB = pathwayDB, method = method)
}
# clusterSim ####
#' Compare two clusters of genes
#'
#' Looks for the similarity between genes in groups
#'
#' Once the pathways for each cluster are found they are combined using
#' combineScores.
#' @param cluster1,cluster2 A vector with genes in \code{id}.
#' @inheritParams genesSim
#' @export
#' @author Lluís Revilla
#' @seealso For a different approach see \code{\link{clustersSim}},
#' \code{\link{combineScores}} and \code{\link{conversions}}
#' @return \code{clusterSim} returns a similarity score of the two clusters
#' @examples
#' library("org.Hs.eg.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- select(org.Hs.eg.db, keys = entrezids, keytype = "ENTREZID",
#'                      columns = "PATH")
#' clusterSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH")
#' clusterSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", NULL)
#' clusterSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", "avg")
clusterSim <- function(cluster1, cluster2, genes, id, pathwayDB, method = "max"){
    pathways1 <- genesInfo(genes, id, cluster1, pathwayDB)
    pathways2 <- genesInfo(genes, id, cluster2, pathwayDB)

    pathways1 <- unlist(pathways1, use.names = FALSE)
    pathways2 <- unlist(pathways2, use.names = FALSE)

    mpathSim(pathways1, pathways2, genes, id, pathwayDB, method)
}

#' @param clusters A list of clusters of genes to be found in \code{id}.
#' @rdname clusterSim
#' @return \code{mclusterSim} returns a matrix with the similarity scores for
#' each cluster comparison.
#' @export
#' @examples
#'
#' clusters <- list(cluster1 = c("18", "81", "10"),
#'                  cluster2 = c("100", "10", "1"),
#'                  cluster3 = c("18", "10", "83"))
#' mclusterSim(clusters, genes.kegg, "ENTREZID", "PATH")
#' mclusterSim(clusters, genes.kegg, "ENTREZID", "PATH", "avg")
mclusterSim <- function(clusters, genes, id, pathwayDB, method = "max") {
    vc <- Vectorize(clusterSim, vectorize.args = c("cluster1", "cluster2"))
    outer(clusters, clusters, vc, genes = genes, id = id,
          pathwayDB = pathwayDB, method = method)

}

# clustersSim ####
#' Compare two clusters of genes
#'
#' Looks for the similarity between genes of a group and then between each
#' group.
#'
#' Differs with clusterSim that first each combination between genes is
#' calculated, and with this values then the comparison between the two
#' clusters is done. Thus applying combineScores twice, one at gene level and
#' another one at cluster level.
#' @param cluster1,cluster2 A vector with genes in \code{id}.
#' @inheritParams genesSim
#' @param method A vector with two  or one argument to be passed to combineScores the
#' first one is used to summarize the similarities of genes, the second one
#' for clusters.
#' @export
#' @author Lluís Revilla
#' @seealso \code{\link{clusterSim}}, \code{\link{combineScores}} and \code{\link{conversions}}
#' @return \code{clustersSim} returns a similarity score of the two clusters or
#' the similarity between the genes of the two clusters.
#' @examples
#' library("org.Hs.eg.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- select(org.Hs.eg.db, keys = entrezids, keytype = "ENTREZID",
#'                      columns = "PATH")
#' clustersSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH")
#' clustersSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", c("avg", "avg"))
#' clustersSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", c("avg", "rcmax.avg"))
#' clus <- clustersSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
#'            "ENTREZID", "PATH", "avg")
#' clus
#' combineScores(clus, "rcmax.avg")
clustersSim <- function(cluster1, cluster2, genes, id, pathwayDB,
                        method = c("max", "rcmax.avg")) {
    if (length(method) > 2L | is.null(method) | any(is.na(method))) {
        stop("Please provide two  or one methods to combine scores.",
             "See Details")
    }

    pathways1 <- genesInfo(genes, id, cluster1, pathwayDB)
    pathways2 <- genesInfo(genes, id, cluster2, pathwayDB)

    vmpathSim <- Vectorize(mpathSim,
                           vectorize.args = c("pathways1", "pathways2"))

    out <- outer(pathways1, pathways2, vmpathSim, genes, id, pathwayDB, NULL)
    out <- apply(out, c(1,2), combineScores, method = method[1L])
    if (length(method) == 2) {
        combineScores(out, method[2L])
    } else {
        out
    }
}

#' @param clusters A list of clusters of genes to be found in \code{id}.
#' @rdname clustersSim
#' @return \code{mclustersSim} returns a matrix with the similarity scores for
#' each cluster comparison.
#' @export
#' @examples
#'
#' clusters <- list(cluster1 = c("18", "81", "10"),
#'                  cluster2 = c("100", "11", "1"),
#'                  cluster3 = c("18", "10", "83"))
#' mclustersSim(clusters, genes.kegg, "ENTREZID", "PATH")
#' mclustersSim(clusters, genes.kegg, "ENTREZID", "PATH", c("max", "avg"))
#' mclustersSim(clusters, genes.kegg, "ENTREZID", "PATH", c("max", "BMA"))
mclustersSim <- function(clusters, genes, id, pathwayDB,
                         method = c("max", "rcmax.avg")) {
    if (length(method) != 2) {
        stop("Please provide two methods to combine scores")
    }
    vc <- Vectorize(clustersSim, vectorize.args = c("cluster1", "cluster2"))
    outer(clusters, clusters, vc, genes = genes, id = id,
          pathwayDB = pathwayDB, method = method)

}

