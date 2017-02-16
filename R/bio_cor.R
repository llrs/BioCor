# Functions required for the bio.cor or bio.cor2 computation.


# pathSim ####
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
#' @seealso Used for \code{\link{geneSim}}, see \code{\link{conversions}} help
#' page to transform Dice score to Jaccard score.
#' @examples
#' genes.id2 <- c("52", "11342", "80895", "57654", "548953", "11586", "45985")
#' genes.id1 <- c("52", "11342", "80895", "57654", "58493", "1164", "1163",
#' "4150", "2130", "159")
#' pathSim(genes.id1, genes.id2)
#' pathSim(genes.id2, genes.id2)
pathSim <- function(g1, g2) {
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

    # intersect <- crossprod(table(stack(x)))
    # y <- matrix(lengths(x), ncol = ncol(intersect), nrow = nrow(intersect))
    # 2*intersect/(y+t(y))
}

# distCor ####
# Not useful because it doesn't have any real causation 02/08/2016
# Keeped in case it might be deemed useful
distCor <- function(a, b, info) {
    use_info <- c("chromosome_name", "strand", "start_position",
                  "end_position", "gene_biotype")
    info_a <- unique(info[info$affy_hg_u133_plus_2 == a, use_info])
    info_b <- unique(info[info$affy_hg_u133_plus_2 == b, use_info])
    if (nrow(info_a) != 1L | nrow(info_b) != 1L) {
        return(NA)
    }
    # Using the position and the type of genes output a cor
    if (info_a["chromosome_name"] != info_b["chromosome_name"]) {
        score <- 0
    }
    score <- 0.5 # If in the same chromosome at least 0.5 distance correlation

    # Maybe the default for the same chromosome could be substituted by the
    # CM distance
    # 5*10^5 is the region of upstream/downstream where regulation usually
    # occurs
    if (info_b["strand"] == info_a["strand"]) {
        start_dist <- info_a["start_position"] - info_b["start_position"]
        end_dist <- info_a["end_position"] - info_b["end_position"]
    } else {
        start_dist <- info_a["start_position"] - info_a["end_position"]
        end_dist <-  info_b["start_position"] - info_b["end_position"]
    }

    if (abs(start_dist) < 5*10 ^ 5) {
        score <- score + 0.25
    } else if (abs(end_dist) < 5*10 ^ 5) {
        score <- score + 0.25
    }

    # Using the cathegory of biotypes to score them
    biotypes <- c("miRNA", "snRNA", "snoRNA", "scaRNA", "lincRNA", "lncRNA")
    if (((info_a["gene_biotype"] %in% biotypes)  &
         (info_b["gene_biotype"] == "protein_coding")) |
        ((info_b["gene_biotype"] %in% biotypes) |
         (info_a["gene_biotype"] == "protein_coding"))) {
        score <- score + 0.25
    }
    if ( (info_a["gene_biotype"] == "processed_pseudogene") |
         (info_b["gene_biotype"] == "processed_pseudogene")) {
        score <- score - 0.25
    }
    return(score)
}


# combBiopath ####
#' Wrapper to expand.grid for pathways of genes
#'
#' Given a data.frame with the information and the genes we want to extract,
#' finds all the pathway for each gene and return a data.frame with the
#' combinations of all pathways of gene 1 with all pathways of gene 2
#' @param comb Pair of ids of genes to find the pathways of. They should be on
#' the column indicated by \code{by} to compare
#' @param info The data.fram with a column for the genes id, and another for
#' pathway
#' @param by Type of columns to subset
#' @param biopath Column of info we want to compare
#' @return A matrix of combinations of all the ids selected from biopath
#' and compare them all
#' @author Lluís Revilla
#' @note Documented but not exported. Used internally for \code{\link{genesSim}}
#' @examples
#' library("org.Hs.eg.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- select(org.Hs.eg.db, keys = entrezids, keytype = "ENTREZID",
#'                      columns = "PATH")
#' combBiopath(c("81", "18"), genes.kegg, "ENTREZID", "PATH")
combBiopath <- function(comb, info, by, biopath) {
    a <- unique(info[info[[by]] == comb[1L], biopath, drop = TRUE])
    a <- a[a != ""]
    a <- a[!is.na(a)]

    b <- unique(info[info[[by]] == comb[2L], biopath, drop = TRUE])
    b <- b[b != ""]
    b <- b[!is.na(b)]

    if (all(sapply(a, is.na)) | all(sapply(b, is.na))) {
        return(0L)
    } else if (length(a) == 0L | length(b) == 0L) {
        return(0L)
    }
    out <- expand.grid(a, b)
    colnames(out) <- comb
    out
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
#' @seealso \code{\link{geneSim}}
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

# genesInfo  ####
#' Extract which genes are from which reactome
#' @param genes is the data.frame with information
#' @param colm is the colum where \code{id} is found
#' @param id is the ids we are looking for in column \code{type} of
#' \code{genes}  data.frame
#' @param type is the column we are looking to keep.
#' @return A vector with the unique identifiers of genes of the \code{type}
#' column
genesInfo <- function(genes, colm, id, type) {
    sapply(id, function(x){
        out <- unique(as.vector(genes[genes[, colm, drop = TRUE] == x,
                                      type, drop = TRUE]))
        out[!is.na(out)]})
}

# genesSim ####
#' Calculates the Dice similarity score of two genes
#'
#' Given the information about the two genes via their pathways, uses the ids
#' of the genes in comb to find the Dice similarity score.
#'
#' If an \code{NA} is returned this means that there isn't information
#' available for any pathways for one of those two genes. Otherwise a number
#' between 0 and 1 (both included) is returned. Note that there isn't a
#' negative value of similarity for pathways correlation.
#' @param gene1 Entrez gene id.
#' @param gene2 Entrez gene id.
#' @param genes is the matrix with the information to calculate the similarity.
#' It is created using select, should contain the column "id" and "pathwayDB".
#' @param id is the column of "genes" where \code{gene1} and \code{gene2} are
#' to be found
#' @param pathwayDB is the column where pathways should be found. It is usually
#'  the name of the database where they come from.
#' @return The highest Dice score of all the combinations of pathways between
#' the two ids compared.
#' @export
#' @author Lluis Revilla
#' @seealso See also \code{\link{conversions}} help page to transform Dice
#' score to Jaccard score.
#' @examples
#' library("org.Hs.eg.db")
#' library("reactome.db")
#' entrezids <- keys(org.Hs.eg.db, keytype = "ENTREZID")
#' #Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
#' # data of June 31st 2011)
#' genes.kegg <- select(org.Hs.eg.db, keys = entrezids, keytype = "ENTREZID",
#'                      columns = "PATH")
#' genesSim("81", "18", genes.kegg, "ENTREZID", "PATH")
#' # Extracts the paths of all genes of org.Hs.eg.db from reactome
#' genes.react <- select(reactome.db, keys = entrezids, keytype = "ENTREZID",
#'                       columns = "REACTOMEID")
#' genesSim("81", "18", genes.react, "ENTREZID", "REACTOMEID")
genesSim <- function(gene1, gene2, genes, id, pathwayDB) {
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

    # Find all the combinations of pathways of the two genes
    react_path <- combBiopath(comb, genes, id, pathwayDB)

    # Check that we have pathways info for this combination
    if (is.null(react_path)) {
        return(NA)
    } else if (length(react_path) == 2L) {
        if (nrow(react_path) == 0L) {
            return(NA)
        }
    } else if (is.na(react_path)) {
        return(NA)
    } else if (react_path == 0L) {
        return(NA)
    }

    # calculate the similarity between each pathway combination
    react <- apply(react_path, 1L, function(x){
        genes_1 <- genesInfo(genes, pathwayDB, x[1L], id)
        genes_2 <- genesInfo(genes, pathwayDB, x[2L], id)
        out <- pathSim(genes_1, genes_2)
        out
    })

    # Calculate the max similarity between the two genes
    if (length(react) != sum(is.na(react))) {
        out <- max(react, na.rm = TRUE)
    } else {
        out <- 0L
    }
    out
}
