# Functions required for the bio.cor or bio.cor2 computation.

".combinadic" <- function(n, r, i) {

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

# Function to estimate how much two graphs overlap by looking if the nodes
# are the same
# Function used in bio.cor2
compare_graphs <- function(g1, g2){

  # Check which case are we using
  if (is(g1, "graph") & is(g2, "graph")) {
    prot1 <- nodes(g1)
    prot2 <- nodes(g2)
    if (length(prot1) == 0 | length(prot2) == 0) {
      return(NA)
    }
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

  score <- (length(intersect(prot1, prot2)))*2/(
    length(prot2) + length(prot1))
  score
}

# function to  calculate the arrows from top to bottom
# ig a graph
# return all the paths from top to bottom
s.path <- function(ig) {
  lfi <- graph::leaves(ig, "in")
  degs <- graph::degree(ig)
  root <- names(degs$outDegree)[degs$outDegree == 0]
  paths <- RBGL::sp.between(ig, lfi, root)
  plens <- Biobase::subListExtract(paths, "length", simplify = TRUE)
  # out <- mean(plens)
  return(plens)
}

# Why just BP and note CC and MF ?? Because we want to join by function
# Calculates the degree of overlap of the GO BP ontologies of entrez ids.
# test genes
# 52 11342
# 52 80895
# 57654 58493
# 1164 1163
# 4150 2130
# 159 52

go_cor <- function(e_a, e_b, chip = "hgu133plus2.db", mapfun = NULL,
                   Ontology = "BP", ...) {
  # https://support.bioconductor.org/p/85702/#85732
  suppressPackageStartupMessages({
    library("GOstats")
    library("org.Hs.eg.db")
  })
  out <- 0.0

  if (is.na(e_a) | is.na(e_b)) {
    return(NA)
  }
  if (e_a == e_b) {
    return(1)
  }
  # Ensure proper format
  e_a <- as.character(e_a)
  e_b <- as.character(e_b)

  if (mapfun) {
    mapfunc <- function(z) {
      mget(z, revmap(org.Hs.egGO2EG), ifnotfound = NA)
    }

    LP <- GOstats::simLL(e_a, e_b, Ontology, measure = "LP", mapfun = mapfunc)
    UI <- GOstats::simLL(e_a, e_b, Ontology, measure = "UI", mapfun = mapfunc)
  } else {
    LP <- GOstats::simLL(e_a, e_b, Ontology, measure = "LP", chip = chip)
    UI <- GOstats::simLL(e_a, e_b, Ontology, measure = "UI", chip = chip)
  }

  if (length(LP) > 1 | length(UI) > 1) {
    if (!is.na(LP["sim"]) & !is.na(UI["sim"])) {

      # Calculates the score taking into account the size and the middle path
      # Taking advantage of the fact that in GO there is a root and leaves
      # UI: Union intersect, is the size of the intersection of the node
      #        sets divided by the size of the union of the node sets
      # LP: longest path, is the longest path in the intersection graph of
      #                the two supplied graph.
      # mean.gx mean number of steps from top of GO DAG to bottom

      # mean.g1 <- s.path(LP$g1)
      # mean.g2 <- s.path(LP$g2)
      # n.root.path <- length(c(mean.g1, mean.g2))
      # mean/LP attemps to normalize the mean length from top to bottom for each
      #  individual graph by the length of the longest path merging both graphs
      out <- UI$sim/LP$sim
      # out1 <- n.root.path/(mean(c(length(mean.g1), length(mean.g2))))/(
      #   UI$sim*LP$sim)
    }
  }
  # warning("score ", out)
  if (out > 1) {
    stop("go_cor is bigger than 1, for genes ", e_a, " and ", e_b)
  }
  # else if (out1 > 1) {
  #   stop("go_cor1 is bigger than 1 for genes ", e_a, " and ", e_b)
  # }
  return(out)
}

# function that given two kegg pathways calculates the similarity
# Just needed on bio.cor not in bio.cor2
kegg_cor <- function(react_a, react_b){
  # Function that correlates based on kegg ids
  # Basically calculates how many nodes do overlap between pathways

  if (is.na(react_a) | is.na(react_b)) {
    return(NA)
  } else if (react_a == react_b) {
    return(1)
  }

  # Retrieve the results from internet https site
  tmp <- tempfile("hsa")
  retrieveKGML(paste0("path:", react_a), "hsa", tmp, method = "wget",
               quiet = TRUE)
  g1 <- parseKGML2Graph(tmp, expandGenes = TRUE)

  tmp2 <- tempfile("hsa")
  retrieveKGML(paste0("path:", react_b), "hsa", tmp2, method = "wget",
               quiet = TRUE)
  g2 <- parseKGML2Graph(tmp2, expandGenes = TRUE)

  score <- compare_graphs(g1, g2)
  score
}

# Not useful because it doesn't have any real causation 02/08/2016
dist_cor <- function(a, b, info){
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

# Function that correlates based on reactome ids
# Just needed on bio.cor not in bio.cor2
react_cor <- function(react_a, react_b, hR){
  # Basically calculates how many nodes do overlap between pathways

  if (is.na(react_a) | is.na(react_b)) {
    return(NA)
  } else if (react_a == react_b) {
    return(1)
  }

  ids <- sapply(hR, function(x){x@id})
  react_name <- names(hR[ids %in% react_a])
  react_name2 <- names(hR[ids %in% react_b])

  if (length(react_name) != 1 | length(react_name2) != 1) {
    return(NA)
  }

  # Obtain each pathway
  g1 <- hR[[react_name]]
  g2 <- hR[[react_name2]]

  score <- compare_graphs(g1, g2)
  score
}

# Funcion to perform efficiently the conversion from combinations
# to symmetric matrix
comb2mat <- function(input, func, ...){
  # Perform all the combinations of 2 from the input
  cobs <- list()
  dots <- list(...)
  # parallel
  # cobs <- foreach(i = 1:length(input), .verbose = T) %dopar% {.combinadic(input, 2, i)}
  for (i in 1:length(input)) {
    cobs[[i]] <- .combinadic(input, 2, i)
  }
  # cobs <- combn(input, 2)
  func <- match.fun(func)
  # cobs <- lapply(seq_len(ncol(cobs)), function(i) func(i[1], i[2], ...)
  # simplify2array(bplapply(cobs, ))
  # p <- DoparParam()
  # N <- bplapply(cobs, function(x){func(x[1], x[2], dots)}, BPPARAM = p)
  # warning(length(N), " testing the N \n", head(N))
  N <- foreach(i = cobs, .combine = c, .verbose = F) %dopar% {
    func(i[1], i[2], mapfun = TRUE, Ontology = "BP")
  }
  # N <- sapply(cobs, function(x){func(x[1], x[2], ...)}) # maybe bplapply
  # Function that performs the calculus
  # N <- seq_len(ncol(combs))
  seq2mat(input, N)
}

# Transform a vector to a symetric matrix
#
# The matrix should be of ncol = x and nrow = x, so dat is at least
# choose(length(x), 2) of length.
# dat is the values or object to fill the matrix with
# x are the colnames and the rownames we want of the resulting matrix.
#
# testthat
# a <- seq2mat(LETTERS[1:5], 1:10)
# isSymmetric(a)
# func(LETTERS[3], LETTERS[2], a) == a[LETTERS[3], LETTERS[2]]
seq2mat <- function(x, dat) {
  if (length(dat) != choose(length(x), 2)) {
    stop("Data is not enough big to populate the matrix")
  }
  out <- matrix(ncol = length(x), nrow = length(x))
  out[upper.tri(out)] <- unlist(dat)
  out[lower.tri(out)] <- t(out)[lower.tri(t(out))]
  diag(out) <- 1
  if (ncol(out) != length(x)) {
    stop("Error, the symetrization doesn't work")
  }
  rownames(out) <- colnames(out) <- x
  return(out)
}

# Extract all the ids of biopath for each element on the combination
# comb ids of the column "by" to compare
# info is the data.fram with information of the by and biopath columns
# by Type of columns to subset
# biopath Column of info we want to compare
# return A matrix of combinations of all the ids selected from biopath
# and compare them all
# Used in bio.cor2
comb_biopath <- function(comb, info, by, biopath) {
  a <- unique(info[info[[by]] == comb[1], biopath])
  a <- a[a != ""]
  a <- a[!is.na(a)]

  b <- unique(info[info[[by]] == comb[2], biopath])
  b <- b[b != ""]
  b <- b[!is.na(b)]

  if (all(sapply(a, is.na)) | all(sapply(b, is.na))) {
    return(NA)
  } else if (length(a) == 0 | length(b) == 0) {
    return(NA)
  }
  expand.grid(a, b)
  # })
}

# Check if something is a matrix (internal use only)
# Not used in bio.cor2
check_na <- function(x){
  if (class(x) != "matrix") {
    if (length(x) == 0) {
      return(TRUE)
    } else if (length(x) == 1) {
      if (is.na(x)) {
        return(TRUE)
      }
    }
  }
  return(FALSE)
}

# Using data correlates biologically two genes
# not used in bio.cor2
bio.cor <- function(names_probes, ... ){
  # Using data correlates biologically two genes or probes
  # From the graphite package
  # x should be entrez id
  # or change the internals from "Entrez Gene" to "Symbols"
  humanReactome <- pathways("hsapiens", "reactome")
  # humanKegg <- pathways("hsapiens", "kegg")

  # mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # attri <- c("affy_hg_u133_plus_2","entrezgene",
  #            "gene_biotype", "start_position", "end_position", "chromosome_name",
  #            "strand", "reactome",
  #            "hgnc_symbol")
  # info <- getBM(attributes = attri,
  #               filters = "affy_hg_u133_plus_2", values = names_probes,
  #               mart = mart)

  # affy.id <- select(hgu133plus2.db, keys = names_probes,
                    # columns = c("PROBEID", "ENTREZID", "SYMBOL", "PATH"))
  gene.symbol <- unique(toTable(org.Hs.egSYMBOL2EG))
  gene.kegg <- unique(toTable(org.Hs.egPATH2EG))
  gene.reactome <- unique(toTable(reactomePATHID2EXTID))
  genes <- merge(gene.symbol, gene.kegg, by.x = "gene_id", by.y = "gene_id",
                 all = TRUE)
  genes <- unique(merge(genes, gene.reactome, all = TRUE))
  colnames(genes) <- c("Entrez Gene", "Symbol", "KEGG", "Reactome")
  go_mat <- comb2mat(names_probes, go_cor, mapfun = TRUE)
  # dist_mat <- comb2mat(x, dist_cor, info) # Not useful

  # TODO: extract a lab notebook
  # comb <- combn(names_probes, 2)
  react.bio <- rep(NA, choose(length(names_probes), 2))
  kegg.bio <- react.bio
  for (i in 1:choose(length(names_probes), 2)) {
    comb <- .combinadic(names_probes, 2, i)
    react_path <- comb_biopath(comb, genes, "Entrez Gene","Reactome")


    if (check_na(react_path)) {
      a <- NA
    } else {
      a <- apply(react_path, 1, function(y){
        react_cor(y[1], y[2], hR = humanReactome)
      })
    }

    # Calculate the max of the correlation of both ids
    if (length(a) != 0 | sum(!is.na(a)) >= 1) {
      react.bio[i] <- max(a, na.rm = TRUE)
    } else {
      react.bio[i] <- NA
    }
    # # })
    # print(head(comb))
    # print(head(affy.id))
    kegg_path <- comb_biopath(comb, genes, "Entrez Gene","KEGG")
    # Calculate the max of the correlation of both ids
    # N <- lapply(kegg_path, function(x){
    #   if (check_na(x)) {
    #     return(NA)
    #   }
    if (check_na(kegg_path)) {
      a <- NA
    } else {
      # print(head(kegg_path))
      a <- apply(kegg_path, 1, function(y) {
        result <- tryCatch({
          kegg_cor(y[1], y[2])
        }, warning = function(w) {
          message(paste("Warning downloading the data for ",
                        y[1], "or", y[2]))
          message(w)
        }, error = function(e) {
          message(paste("Couldn't download the data for ", y[1], "or", y[2]))
          # Choose a return value in case of error
          return(NA)
        })
      }
      )
    }

    if (length(a) != 0) {
      kegg.bio[i] <- max(a, na.rm = TRUE)
    } else {
      kegg.bio[i] <- NA
    }
    # }})

  }
  react_mat <- seq2mat(names_probes, react.bio)
  kegg_mat <- seq2mat(names_probes, kegg.bio)
  cor_mat <- list(reactome = react_mat, kegg = kegg_mat, go = go_mat)

}

# Not used in bio.cor2 but in cor.all
weighted <- function(x, w){
  if (length(x) != length(w)) {
    stop(paste("Weights and data don't match the length.\n",
               length(x), length(w)))
  }
  sum(x*w, na.rm = TRUE)
}

# adj Adjacency matrix
# bio_mat bio_mat object
rm.dup <- function(adj, bio_mat) {
  if (ncol(adj) != ncol(bio_mat[[1]])) {
    ids <- select(org.Hs.eg.db, keys = colnames(adj), keytype = "SYMBOL",
                  column = "ENTREZID")
    dup_ids <- ids[duplicated(ids$SYMBOL), "ENTREZID"]
    sapply(bio_mat, function(x){
      dupli <- colnames(x) %in% dup_ids
      colums_d <- x[,  dupli]
      if ((sum(is.na(colums_d)) | sum(colums_d == 0)) == (ncol(adj) - 1)) {
        merge(t(ids), colnames(x))
      }
    })
  }
}

# Function that used the previously calculated biological correlation to
# calculate the total correlation
# x is the datExpresssion
cor.all <- function(x, bio_mat, weights = c(0.5, 0.18, 0.10, 0.22), ...){
  # exp, reactome, kegg, go
  # cor_mat <- cor(x, use = "p")
  if (sum(weights) > 1) {
    stop("Weights are too big. The sum must be equal to 1")
  } else if (sum(weights) < 1) {
    warning("Weights are smaller than 1.")
  }
  cors <- c(list(exp = x), bio_mat)
  # Could use parallel
  # parApply(cl, simplify2array(cors), c(1, 2), weights, w = weights)
  apply(simplify2array(cors), c(1,2), weighted, w = weights)
}

# Builds a graph of the kegg pahtways known
# Not used in bio.cor and neigher in bio.cor2
kegg_build <- function(entrez_id){
  # We can build it without downloading from internet following this link
  # https://www.biostars.org/p/2449/#5887
  ## map -- KEGG id links the KEGG and org.Hs.eg packages
  name2id <- toTable(KEGGPATHNAME2ID)
  id2gene <- toTable(revmap(org.Hs.egPATH))
  paths <- unique(id2gene[id2gene$gene_id == entrez_id, "path_id"])
  # add gene SYMBOL for each Entrez id
  id2gene$symbol <- unlist(mget(id2gene$gene_id, org.Hs.egSYMBOL))
  path2gene <- merge(name2id, id2gene, by = "path_id") # 'join'

  met <- sapply(paths, function(x){
    ## create a graphBAM instance for each pathway
    pathway <- path2gene[path2gene$path_id == x, ]
    df <- with(pathway, data.frame(from = symbol, to = symbol,
                                     weight = 1,
                                     stringsAsFactors = FALSE))
    df <- unique(df)
    gr <- graphBAM(df, edgemode = "directed")
    gr <- as(gr,"graphNEL")
    nodes(gr)
  })
  met

}

# Finds the indices of the duplicated events of a vector
# The output is determined by the sapply function
indices.dup <- function(vec) {
  sapply(unique(vec[duplicated(vec)]), function(x){
     b <- 1:length(vec)
     b[vec == x]})
}


# Comparing information in databases
#
# Calculates the overlap or correlation of the information available in several
# systems.
#
# For Gene Ontologies, the DAG path structure is used to compute how similar two
# genes are. For metabolic pathways the max number of proteins involved in a
# pathway for each gene is calculated.
#
#
# genes_id is vector of ids to compare
# ids indicate if the id is eihter Entre Gene or Symbol
# go logical; indicates if the overlap in terms of GO is calculated
# react logical; indicates if the overlap in Reactome pathway is calculated
# kegg logical; indicates if the overlap in Kegg database is calculated
bio.cor2 <- function(genes_id, ids = "Entrez Gene",
                     go = FALSE, react = TRUE, kegg = FALSE, all = FALSE) {
  # Using data correlates biologically two genes or probes
  # From the graphite package
  # x should be entrez id
  # or change the internals from "Entrez Gene" to "Symbol"
  if (!ids %in% c("Entrez Gene", "Symbol")) {
    stop("Please check the input of genes in Symbol or Entrez Gene format")
  }
  if (all) {
    go <- kegg <- react <- all
  }

  # Obtain data from the annotation packages
  if (ids == "Symbol") {
    gene.symbol <- select(org.Hs.eg.db, keys = genes_id,
                          keytype = "SYMBOL", column = "ENTREZID")
    colnames(gene.symbol) <- c("Symbol", "Entrez Gene")
  } else {
    gene.symbol <- select(org.Hs.eg.db, keys = genes_id,
                          keytype = "ENTREZID", column = "SYMBOL")
    colnames(gene.symbol) <- c("Entrez Gene", "Symbol")
  }
  n.combin <- choose(length(gene.symbol$`Entrez Gene`), 2)
  # FIXME the length of input and for the calculus may differ!!!!
  if (sum(is.na(gene.symbol$`Entrez Gene`)) >= 1) {
    message("Some symbols are not mapped to Entrez Genes IDs")
  }
  dup_symb <- duplicated(gene.symbol$Symbol[!is.na(gene.symbol$`Entrez Gene`)])
  if (sum(dup_symb) >= 1) {
    message("Some symbols are mapped to several Entrez Genes IDs.")
  }

  if (kegg) {
    # Obtain data
    gene.kegg <- select(org.Hs.eg.db, keys = gene.symbol$`Entrez Gene`,
                        keytype = "ENTREZID", columns = "PATH")
    colnames(gene.kegg) <- c("Entrez Gene", "KEGG") # Always check it!
    # Merge data
    genes <- unique(merge(gene.symbol, gene.kegg, all = TRUE, sort = FALSE))
  }

  if (react) {
    if (!kegg) {
      genes <- gene.symbol
    }
    gene.reactome <- select(reactome.db, keys = gene.symbol$`Entrez Gene`,
                            keytype = "ENTREZID", columns = "REACTOMEID")
    colnames(gene.reactome) <- c("Entrez Gene", "Reactome")
    genes <- unique(merge(genes, gene.reactome, all = TRUE, sort = FALSE))
  }
  # detach(package:GOstats, unload = TRUE, character.only = TRUE)
  # detach(package:org.Hs.eg.db, unload = TRUE, character.only = TRUE)

  # Calculate the GO, or pathways overlap
  if (go) {  # parallel # to run non parallel transform the %dopar% into %do%
    message("Calculating GO information")
    go.mat <- foreach(i = seq_len(n.combin), .verbose = TRUE) %dopar% {
                        comb <- .combinadic(gene.symbol$`Entrez Gene`, 2, i)
                        # message("new comb ", paste(comb))
                        score <- go_cor(comb[1], comb[2], mapfun = TRUE,
                                        Ontology = "BP")
                        score
      # go.mat <- c(go.mat, score)
    }
    # registerDoParallel(4)
    # print(go.mat)
    go_mat <- seq2mat(gene.symbol$`Entrez Gene`, go.mat)
    # go_mat <- comb2mat(gene.symbol$`Entrez Gene`, func = go_cor, mapfun = TRUE,
                       # Ontology = "BP")
    # go_mat.mf <- comb2mat(genes_id, go_cor, mapfun = TRUE, Ontology = "MF")
    # go_mat.cc <- comb2mat(genes_id, go_cor, mapfun = TentrRUE, Ontology = "CC")
    if (sum(!is.na(go_mat))  == length(gene.symbol$`Entrez Gene`)) {
      warning("GO didn't found relevant information!")
    }
    message("GO information has been calculated")
  }

  if (kegg) {  # parallel # to run non parallel transform the %dopar% into %do%
    message("Calculating KEGG information")
    kegg.bio <- foreach(i = seq_len(n.combin), .combine = c,
                        .verbose = F) %dopar% {
      comb <- .combinadic(gene.symbol$`Entrez Gene`, 2, i)
      react_genes(comb, genes, "KEGG", "Entrez Gene")
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
      react_genes(comb, genes, "Reactome", "Entrez Gene")
    }

    if (sum(!is.na(react.bio)) == length(genes_id)) {
      warning("REACTOME didn't found relevant informationª\n")
    }
    react_mat <- seq2mat(gene.symbol$`Entrez Gene`, react.bio)
    message("REACTOME information has been calculated")
  }


  if (all) {
    cor_mat <- list(reactome = react_mat, kegg = kegg_mat, go = go_mat)
  } else if (kegg & react) {
    cor_mat <- list(reactome = react_mat, kegg = kegg_mat)
  } else if (kegg & go) {
    cor_mat <- list(kegg = kegg_mat, go = go_mat)
  } else  if (go & react) {
    cor_mat <- list(reactome = react_mat, go = go_mat)
  } else if (go) {
    cor_mat <- list(go = go_mat)
  } else if (react) {
    cor_mat <- list(react = react_mat)
  } else if (kegg) {
    cor_mat <- list(kegg = kegg_mat)
  }

  # Select duplicate ids
  dupli <- indices.dup(gene.symbol[, ids])
  if (is.matrix(dupli)) {
    dupli2 <- lapply(seq_len(ncol(dupli)), function(i) {
      dupli[,i]})
    names(dupli2) <- colnames(dupli)
    dupli <- dupli2
  }
  if (length(dupli) >= 1) {
    # Keep the interesting columns
    message("Removing duplicated columns")
    cor_mat <- Map(function(mat, x = dupli) {
      # Select the colums to delete
      rem.colum <- sapply(x, function(y, m) {
        mean.column <- apply(m[, y], 2, mean, na.rm = TRUE)
        i <- which.max(mean.column)
        # Select those who don't bring more information
        rem.colum <- setdiff(y, y[i])
      }, m = mat)

      mat[-rem.colum, -rem.colum]
    }, cor_mat)
  }

  return(cor_mat)
}

# Extract which genes are from which reactome
# genes is the data.frame with information
# colm is the colum where "id" is found
# id is the ids we are looking for in column Symbol of such data.frame
genes.info <- function(genes, colm, id) {
  # Genes is the df, colm is the column you want, id is the id of the pathway
  out <- unique(genes[genes[[colm]] == id, "Symbol"])
  out[!is.na(out)]
}

# Calculates for all the combinations of pathway the right score
# comb are the ids to compare
# genes is the matrix with the information about "id" and "react"
# id is the column of "genes" where comb are to be found
# react is the column where pathways should be found.
# returns the max score of comparing the ids of such reaction
react_genes <- function(comb, genes, react, id) {
  if (!react %in% colnames(genes)) {
    stop("Please check which type of reaction you want")
  }
  if (any(is.na(comb))) {
    return(NA)
  }
  react_path <- comb_biopath(comb, genes, id, react)

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

  react <- apply(react_path, 1, function(x){
    a <- genes.info(genes, react, x[1])
    b <- genes.info(genes, react, x[2])
    out <- compare_graphs(a, b)
    out
  })

  # If NA returns 0
  if (length(react) != sum(is.na(react))) {
    out <- max(react, na.rm = TRUE)
  } else {
    out <- 0
  }
  out
}

#
# Rprof(tmp <- tempfile())
# a <- bio.cor2(as.character(1:50), kegg = TRUE)
# a
# Rprof()
# summaryRprof(tmp)
