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
#' @import BiocParallel
#' @import org.Hs.eg.db
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
        #                             kegg_mat[i, j] <- geneSim(
        #                                 c(i, j), genes,"ENTREZID", "PATH")
        #                         }
        #                     }
        kegg.bio <- bpmapply(function(x){
            comb <- combinadic(n = orig.ids, r = 2, i = x)
            geneSim(comb[1], comb[2], genes, "ENTREZID", "PATH")},
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
            geneSim(comb[1], comb[2], genes, "ENTREZID", "REACTOMEID")},
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
