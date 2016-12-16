# Functions required for the bio.cor or bio.cor2 computation.

# combinadic ####
#' i-th combination of n elements taken from r
#'
#' Function similar to combn but for larger vectors
#' @param n Elements to extract the combination from
#' @param r Number of elements per combination
#' @param i ith combination
#' @rdname combinadic
#' @seealso \code{\link{combn}}
#' @export
#' @import doParallel
#' @import foreach
combinadic <- function(n, r, i) {

  # http://msdn.microsoft.com/en-us/library/aa289166(VS.71).aspx
  # http://en.wikipedia.org/wiki/Combinadic
  n0 <- length(n)
  if (i < 1 | i > choose(n0,r)) {
    stop("'i' must be 0 < i <= n0!/(n0-r)!")
  }
  largestV <- function(n, r, i) {
    v <- n # Adjusted for one-based indexing
    while (choose(v,r) >= i) { # Adjusted for one-based indexing
      v <- v - 1
    }
    return(v)
  }

  res <- rep(NA,r)
  for (j in 1:r) {
    res[j] <- largestV(n0,r,i)
    i <- i - choose(res[j],r)
    n0 <- res[j]
    r <- r - 1
  }
  res <- res + 1
  res <- n[res]
  return(res)
}

# comparePathways ####
#' Compare pathways
#'
#' Function to estimate how much two graphs overlap by looking how much of the
#' nodes are shared. It also works with list of genes.
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

# function to  calculate the arrows from top to bottom
# ig a graph
# return all the paths from top to bottom

# Why just BP and note CC and MF ?? Because we want to join by function
# Calculates the degree of overlap of the GO BP ontologies of entrez ids.
# test genes
# 52 11342
# 52 80895
# 57654 58493
# 1164 1163
# 4150 2130
# 159 52


# Not useful because it doesn't have any real causation 02/08/2016
# Keeped in case it might be deemed useful
dist_cor <- function(a, b, info) {
  use_info <- c("chromosome_name", "strand", "start_position", "end_position",
                "gene_biotype")
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
  # 5*10^5 is the region of upstream/downstream where regulation usually occurs
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

# seq2mat ####
#' Transform a vector to a symetric matrix
#'
#' The matrix should be of \code{ncol = length(x)} and \code{nrow = length(x)},
#' so \code{dat} is at least \code{choose(length(x), 2)} of length.
#'
#' It assumes that the data provided comes from using the row and column id to
#' obtain it.
#' @param x names of columns and rows, used to define the sieze of the matrix
#' @param dat Data to fill with the matrix
#' @return A matrix of length equal to \code{x} with the diagonal set to 1, and
#' \code{dat} on the upper and lower triangle.
#' @examples
#' seq2mat(LETTERS[1:5], 1:10)
#' @export
seq2mat <- function(x, dat) {
  if (length(dat) != choose(length(x), 2)) {
    stop("Data is not enough big to populate the matrix")
  }
  out <- matrix(ncol = length(x), nrow = length(x))
  out[upper.tri(out)] <- unlist(dat)
  out[lower.tri(out)] <- t(out)[lower.tri(t(out))]
  diag(out) <- 1
  rownames(out) <- colnames(out) <- x
  return(out)
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

# weighted ####
#' Calculates the weighted sum of the values
#'
#' Each values should have its weight, otherwise it will throw an error.
#' @param x Vector of numbers
#' @param weights Vector of weights
#' @return A number product of x*weights removing all \code{NA} values
#' @export
weighted <- function(x, weights) {
  if (length(x) != length(weights)) {
    stop("Weights and data don't match the length: ", length(x), " != ",
         length(weights))
  }

  if (!is.numeric(x) | !is.numeric(weights)) {
    stop("weights and x should be numeric")
  }
  if (sum(weights) > 1) {
    warning("The sum of the weights is above 1")
  }
  sum(x*weights, na.rm = TRUE)
}

# removeDup ####
#' Remove duplicated rows and columns
#'
#' Given the indices of the duplicated entries remove the columns and rows until
#'  just one is left, it keeps the duplicated with the highest absolute mean
#' value.
#'
#' @param cor_mat List of matrices
#' @param dupli List of indicies with duplicated entries
#' @return A matrix with only one of the columns and rows duplicated
#' @export
removeDup <- function(cor_mat, dupli) {
  if (!all(sapply(cor_mat, isSymmetric))) {
    stop("All the matrices of mat should be symmetric and with the same column",
         "names and rownames")
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

# addSimilarities ####
#' Additive integration of similarities
#'
#' Function that used the previously calculated similarities into a single
#' similarity matrix.
#'
#' The total weight sum can't be higher than 1 to prevent values above 1 but can be
#' below 1.
#' @param x A matrix with the similarity of expression
#' @param bio_mat A list of matrices of the same dimension
#' @param weights A numeric vector of weight to multiply each similarity
#' @return A similarity matrix
#' @export
addSimilarities <- function(x, bio_mat, weights = c(0.5, 0.18, 0.10, 0.22)){
  # exp, reactome, kegg, go
  # cor_mat <- cor(x, use = "p")
  if (sum(weights) > 1) {
    stop("Weights are too big. The sum must be equal to 1")
  } else if (sum(weights) < 1) {
    warning("Weights are smaller than 1.")
  }

  if (!is.matrix(x)) {
    stop("Expected a matrix, generally a similarity measure from expression")
  }
  if (!all(dim(x) == dim(bio_mat[[1]]))) {
    stop("Dimensions of x and bio_mat matrices is different")
  }
  cors <- c(list(exp = x), bio_mat)

  # Apply weighted to each cell position of each similarity measure
  apply(simplify2array(cors), c(1,2), weighted, w = weights)
}
# duplicateIndices ####
#' Finds the indices of the duplicated events of a vector
#'
#' Finds the indices of duplicated elements in the vector given.
#'
#' For each duplication it can return a list or if all the duplication events
#' are of the same length it returns a matrix, where each column is duplicated.
#' @param vec Vector of identififiers presumably duplicated
#' @return The format is determined by the simplify2array
#' @export
duplicateIndices <- function(vec) {
  if (!is.character(vec)) {
    stop("Expected a list of characters to find duplicates on it")
  }
  sapply(unique(vec[duplicated(vec)]), function(x){
     b <- 1:length(vec)
     b[vec == x]}, simplify = FALSE)
}

# BioCor ####
#' Comparing information in databases
#'
#' Calculates a functional similarity of genes  using information available in
#' several databases.
#'
#' For Gene Ontologies, the DAG path structure is used to compute how similar two
#' genes are. For metabolic pathways the max number of proteins involved in a
#' pathway for each gene is calculated.
#'
#' @param genes_id is vector of ids to compare
#' @param ids indicate if the id is eihter Entre Gene or Symbol
#' @param go logical; indicates if the overlap in terms of GO is calculated
#' @param react logical; indicates if the overlap in Reactome pathway is
#' calculated
#' @param kegg logical; indicates if the overlap in Kegg database is calculated
#' @param all logical; indicates if all the previous (go, react, and kegg),
#' should be set to TRUE
#' @return A list where each element is a matrix with the similarity for such
#' database
#' @export
#' @importFrom reactome.db reactome.db
#' @import org.Hs.eg.db
#' @importFrom AnnotationDbi select
bioCor <- function(genes_id, ids = "Entrez Gene",
                     go = FALSE, react = TRUE, kegg = FALSE, all = FALSE) {
  if (!ids %in% c("Entrez Gene", "Symbol")) {
    stop("Please check the input of genes in Symbol or Entrez Gene format")
  }
  if (all) {
    go <- kegg <- react <- all
  }

  # Obtain data from the annotation packages
  if (ids == "Symbol") {
    gene.symbol <- suppressMessages(select(org.Hs.eg.db, keys = genes_id,
                                           keytype = "SYMBOL",
                                           columns = "ENTREZID"))
    colnames(gene.symbol) <- c("Symbol", "Entrez Gene")
  } else {
    gene.symbol <- suppressMessages(select(org.Hs.eg.db, keys = genes_id,
                                           keytype = "ENTREZID",
                                           columns = "SYMBOL"))
    colnames(gene.symbol) <- c("Entrez Gene", "Symbol")
  }
  n.combin <- choose(length(gene.symbol$`Entrez Gene`), 2)
  # FIXME the length of input and for the calculus may differ!!!!
  if (sum(is.na(gene.symbol$`Entrez Gene`)) >= 1) {
    message("Some symbols are not mapped to Entrez Genes IDs")
  }
  dup_symb <- duplicated(gene.symbol$Symbol[!is.na(gene.symbol$`Entrez Gene`)])
  if (sum(dup_symb) >= 1 & ids == "Symbol") {
    message("Some symbols are mapped to several Entrez Genes IDs.")
  } else if (sum(dup_symb) >= 1 & ids == "Entrez Gene") {
    message("Some Entrez Genes IDs are mapped to several symbols.")
  }
  # Obtain the data of kegg and Reactome pathways
  if (kegg) {
    # Obtain data
    gene.kegg <- suppressMessages(
      select(org.Hs.eg.db, keys = gene.symbol$`Entrez Gene`,
             keytype = "ENTREZID", columns = "PATH"))
    colnames(gene.kegg) <- c("Entrez Gene", "KEGG") # Always check it!
    # Merge data
    genes <- unique(merge(gene.symbol, gene.kegg, all = TRUE, sort = FALSE))
  }

  if (react) {
    if (!kegg) {
      genes <- gene.symbol
    }
    gene.reactome <- suppressMessages(select(reactome.db,
                                             keys = gene.symbol$`Entrez Gene`,
                                             keytype = "ENTREZID",
                                             columns = "REACTOMEID"))
    colnames(gene.reactome) <- c("Entrez Gene", "Reactome")
    genes <- unique(merge(genes, gene.reactome, all = TRUE, sort = FALSE))
  }

  # Calculate the GO, or pathways overlap
  if (go) {  # parallel # to run non parallel transform the %dopar% into %do%
    message("Calculating GO information")
    # go.mat <- foreach(i = seq_len(n.combin), .verbose = TRUE) %dopar% {
    #                     comb <- .combinadic(gene.symbol$`Entrez Gene`, 2, i)
    #                     # message("new comb ", paste(comb))
    #                     score <- go_cor(comb[1], comb[2], mapfun = TRUE,
    #                                     Ontology = "BP")
    #                     score
      # go.mat <- c(go.mat, score)
    # }
    # registerDoParallel(4)
    # print(go.mat)
    # go_mat <- seq2mat(gene.symbol$`Entrez Gene`, go.mat)
    # go_mat <- comb2mat(gene.symbol$`Entrez Gene`, func = go_cor, mapfun = TRUE,
                       # Ontology = "BP")
    # go_mat.mf <- comb2mat(genes_id, go_cor, mapfun = TRUE, Ontology = "MF")
    # go_mat.cc <- comb2mat(genes_id, go_cor, mapfun = TentrRUE, Ontology = "CC")
    # if (sum(!is.na(go_mat))  == length(gene.symbol$`Entrez Gene`)) {
    #   warning("GO didn't found relevant information!")
    # }
    message("GO information has been calculated")
  }

  if (kegg) {  # parallel # to run non parallel transform the %dopar% into %do%
    message("Calculating KEGG information")
    kegg.bio <- foreach(i = seq_len(n.combin), .combine = c,
                        .verbose = F) %dopar% {
      comb <- .combinadic(gene.symbol$`Entrez Gene`, 2, i)
      corPathways(comb, genes, "Entrez Gene", "KEGG")
    }

    if (sum(!is.na(kegg.bio)) == length(genes_id)) {
      warning("React didn't found relevant information!\n")
    }
    kegg_mat <- seq2mat(gene.symbol$`Entrez Gene`, kegg.bio)
    message("KEGG information has been calculated")
  }

  if (react) {  # parallel # to run non parallel transform the %dopar% into %do%
    message("Calculating REACTOME information")
    react.bio <- foreach(i = seq_len(n.combin), .combine = c,
                         .verbose = F) %dopar% {
      comb <- .combinadic(gene.symbol$`Entrez Gene`, 2, i)
      corPathways(comb, genes, "Entrez Gene", "Reactome")
    }

    if (sum(!is.na(react.bio)) == length(genes_id)) {
      warning("REACTOME didn't found relevant information!\n")
    }
    react_mat <- seq2mat(gene.symbol$`Entrez Gene`, react.bio)
    message("REACTOME information has been calculated")
  }


  if (kegg & react) {
    cor_mat <- list(reactome = react_mat, kegg = kegg_mat)
  # } else if (kegg & go) {
  #   cor_mat <- list(kegg = kegg_mat, go = go_mat)
  # } else  if (go & react) {
  #   cor_mat <- list(reactome = react_mat, go = go_mat)
  # } else if (go) {
  #   cor_mat <- list(go = go_mat)
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
#' @param id is the ids we are looking for in column \code{type} of \code{genes}
#'  data.frame
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
    return(0)
  } else if (length(react_path) == 2) {
    if (nrow(react_path) == 0) {
      return(0)
    }
  } else if (is.na(react_path)) {
    return(0)
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
