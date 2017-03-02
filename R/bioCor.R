# bioCor ####
#' Comparing information in databases
#'
#' Calculates a functional similarity of genes  using information available in
#' several databases.
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
        message("Some symbols are not mapped to Entrez Genes IDs ",
                "and will be omitted")
    }
    dup_symb <- duplicated(gene.symbol$SYMBOL[
        !is.na(gene.symbol$ENTREZID)])

    if (sum(dup_symb) >= 1L & ids == "SYMBOL") {
        message("Some symbols are mapped to several Entrez Genes IDs.")
    } else if (sum(dup_symb) >= 1L & ids == "Entrez Gene") {
        message("Some Entrez Genes IDs are mapped to several symbols.")
    }

    if (kegg) {
        gene.kegg <- as.list(org.Hs.egPATH)
    }

    if (react) {
        gene.reactome <- as.list(reactomeEXTID2PATHID)
    }

    if (kegg) {  # parallel # to run non parallel transform the %dopar% into
        message("Calculating KEGG similarities")
        kegg.bio <- bpmapply(function(x){
            comb <- combinadic(n = orig.ids, r = 2, i = x)
            geneSim(comb[1], comb[2], gene.kegg)},
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
            geneSim(comb[1], comb[2], gene.reactome)},
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
