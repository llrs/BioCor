# Functions required for the bio.cor or bio.cor2 computation.


# comparePathways ####
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
comparePathways <- function(g1, g2) {

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
    if (length(prot1) == 0 | length(prot2) == 0) {
        return(0)
    }
    score <- (length(intersect(prot1, prot2)))*2/(
        length(prot2) + length(prot1))
    score
}

# distCor ####
# Not useful because it doesn't have any real causation 02/08/2016
# Keeped in case it might be deemed useful
distCor <- function(a, b, info) {
    use_info <- c("chromosome_name", "strand", "start_position",
                  "end_position", "gene_biotype")
    info_a <- unique(info[info$affy_hg_u133_plus_2 == a, use_info])
    info_b <- unique(info[info$affy_hg_u133_plus_2 == b, use_info])
    if (nrow(info_a) != 1 | nrow(info_b) != 1) {
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
#' Finds all the pathways of genes
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
#' @export
combBiopath <- function(comb, info, by, biopath) {
    a <- unique(info[info[[by]] == comb[1], biopath])
    a <- a[a != ""]
    a <- a[!is.na(a)]

    b <- unique(info[info[[by]] == comb[2], biopath])
    b <- b[b != ""]
    b <- b[!is.na(b)]

    if (all(sapply(a, is.na)) | all(sapply(b, is.na))) {
        return(0)
    } else if (length(a) == 0 | length(b) == 0) {
        return(0)
    }
    expand.grid(a, b)
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
removeDup <- function(cor_mat, dupli) {
    if (!all(sapply(cor_mat, isSymmetric))) {
        stop("All the matrices of mat should be symmetric and with the same ",
             "column names and rownames")
    }
    cor_mat <- Map(function(mat, x = dupli) {
        rem.colum <- sapply(x, function(y, m) {
            mean.column <- apply(m[, y], 2, mean, na.rm = TRUE)
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
#' @param ids indicate if the id is eihter Entre Gene or Symbol
#' @param react logical; indicates if the similarities in Reactome pathway is
#' calculated
#' @param kegg logical; indicates if the similarities in Kegg database
#' is calculated
#' @param all logical; indicates if all the previous (go, react, and kegg)
#' similarity measures, should be set to TRUE.
#' @return A list where each element is a matrix with the similarity for such
#' database \code{NA} indicates that there isn't any information of one of
#' those genes.
#' @importFrom reactome.db reactome.db
#' @importFrom AnnotationDbi select
#' @importFrom AnnotationDbi keys
#' @import org.Hs.eg.db
#' @import foreach
#' @export
bioCor <- function(genes_id, ids = "ENTREZID", react = TRUE, kegg = FALSE,
                   all = FALSE) {
    if (!ids %in% c("ENTREZID", "SYMBOL")) {
        stop("Please check the input of genes in ENTREZID or SYMBOL format")
    }
    if (all) {
        kegg <- react <- all
    }

    # Obtain data from the annotation packages
    if (ids == "Symbol") {
        gene.symbol <- suppressMessages(select(org.Hs.eg.db, keys = genes_id,
                              keytype = "SYMBOL", columns = "ENTREZID"))
    } else {
        gene.symbol <- suppressMessages(select(org.Hs.eg.db, keys = genes_id,
                              keytype = "ENTREZID", columns = "SYMBOL"))
    }

    n.combin <- choose(length(gene.symbol$ENTREZID), 2)
    orig.ids <- gene.symbol$ENTREZID

    if (sum(is.na(gene.symbol$ENTREZID)) >= 1) {
        message("Some symbols are not mapped to Entrez Genes IDs")
    }
    dup_symb <- duplicated(gene.symbol$SYMBOL[
        !is.na(gene.symbol$ENTREZID)])

    if (sum(dup_symb) >= 1 & ids == "Symbol") {
        message("Some symbols are mapped to several Entrez Genes IDs.")
    } else if (sum(dup_symb) >= 1 & ids == "Entrez Gene") {
        message("Some Entrez Genes IDs are mapped to several symbols.")
    }

    # Obtain the data of kegg and Reactome pathways
    if (kegg) {
        # Obtain data
        gene.kegg <- suppressMessages(select(org.Hs.eg.db,
                                             keys = gene.symbol$ENTREZID,
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
                                    keys = gene.symbol$ENTREZID,
                                    keytype = "ENTREZID",
                                    columns = "REACTOMEID"))
        }

        genes <- unique(merge(genes, gene.reactome, all = TRUE, sort = FALSE))
    }

    if (kegg) {  # parallel # to run non parallel transform the %dopar% into
        # %do%
        message("Calculating KEGG similarities")
        kegg.bio <- foreach(i = seq_len(n.combin), .combine = c,
                            .verbose = FALSE) %dopar% {
                                comb <- combinadic(orig.ids, 2, i)
                                corPathways(comb, genes, "ENTREZID", "PATH")
                            }

        if (all(is.na(kegg.bio))) {
            warning("KEGG didn't found relevant similarities!")
        }
        kegg_mat <- seq2mat(orig.ids, kegg.bio)
        message("KEGG similarities has been calculated")
    }

    if (react) {  # parallel # to run non parallel transform the %dopar% into
        # %do%
        message("Calculating REACTOME similarities")
        react.bio <- foreach(i = seq_len(n.combin),
                             .combine = c, .verbose = FALSE) %dopar% {

                                 comb <- combinadic(orig.ids, 2, i)
                                 corPathways(comb, genes, "ENTREZID",
                                             "REACTOMEID")
                             }

        if (all(is.na(react.bio))) {
            warning("REACTOME didn't found relevant similarities!")
        }
        react_mat <- seq2mat(orig.ids, react.bio)
        message("REACTOME similarities has been calculated")
    }

    if (kegg & react) {
        cor_mat <- list(reactome = react_mat, kegg = kegg_mat)
    } else if (react) {
        cor_mat <- list(react = react_mat)
    } else if (kegg) {
        cor_mat <- list(kegg = kegg_mat)
    }

    # To match the ncol and nrow of the input with the initial input
    dupli <- duplicateIndices(gene.symbol[, ids])

    # Keep the interesting columns
    if (length(dupli) >= 1) {
        cor_mat <- removeDup(cor_mat, dupli)
    }

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
#' @export
genesInfo <- function(genes, colm, id, type) {
    out <- unique(genes[genes[[colm]] == id, type])
    out[!is.na(out)]
}
# corPathways ####
#' Calculates the similarity score of patwhays
#'
#' Given the information about the relationship between the genes and the
#' pathways, uses the ids of the genes in comb to find the similarity score in
#' them.
#'
#' If an \code{NA} is returned this means that there isn't information
#' available for any pathways for one of those two genes. Otherwise a number
#' between 0 and 1 (both included) is returned. Note that there isn't a
#' negative value of similarity for pathways correlation.
#' @param comb are the ids to compare, it is expected two ids of genes.
#' @param genes is the matrix with the information about "id" and "react"
#' @param id is the column of "genes" where \code{comb} are to be found
#' @param pathwayDB is the column of \code{genes} where pathways should be
#' found. It is usually the name of the database where they come from.
#' @return The highest score of all the combinations of pathways between the
#' two ids compared.
#' @export
corPathways <- function(comb, genes, id, pathwayDB) {
    if (!pathwayDB %in% colnames(genes)) {
        stop("Please check which type of pathway do you want")
    }
    if (any(is.na(comb))) {
        return(NA)
    }
    if (comb[1] == comb[2]) {
        return(1)
    }
    if (length(comb) > 2) {
        stop("comb can only be of length 2, to compare pairs of genes")
    }

    # Find all the combinations of pathways of the two genes
    react_path <- combBiopath(comb, genes, id, pathwayDB)

    # Check that we have pathways info for this combination
    if (is.null(react_path)) {
        return(NA)
    } else if (length(react_path) == 2) {
        if (nrow(react_path) == 0) {
            return(NA)
        }
    } else if (is.na(react_path)) {
        return(NA)
    } else if (react_path == 0) {
        return(NA)
    }

    # calculate the similarity between each pathway combination
    react <- apply(react_path, 1, function(x){
        genes_1 <- genesInfo(genes, pathwayDB, x[1], id)
        genes_2 <- genesInfo(genes, pathwayDB, x[2], id)
        out <- comparePathways(genes_1, genes_2)
        out
    })

    # Calculate the max similarity between the two genes
    if (length(react) != sum(is.na(react))) {
        out <- max(react, na.rm = TRUE)
    } else {
        out <- 0
    }
    out
}
