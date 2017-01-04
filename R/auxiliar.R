# combinadic ####
#' i-th combination of n elements taken from r
#'
#' Function similar to combn but for larger vectors
#' @param n Elements to extract the combination from
#' @param r Number of elements per combination
#' @param i ith combination
#' @return The combination ith of the elements
#' @rdname combinadic
#' @seealso \code{\link{combn}}
#' @export
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
# similarities ####
#' Similarities
#'
#' Function to join list of similarities by a function provided by the user.
#' @param sim list of similarities to be joined. All similarities must have the
#'  same dimensions, and in the same order.
#' @param func function to perform on those similarities: \code{prod},
#' \code{sum}... It should accept as many arguments as similarities matrices
#' are provided, and should use with numbers.
#' @param ... Other arguments passed to the function \code{func}. Usually na.rm
#' or similar.
#' @return A matrix of the size of the similarities
#' @export
similarities <- function(sim, func, ...) {
    # Check that all the matrices are of the same dimensions and squared
    if (length(unique(as.vector(sapply(sim, dim)))) >= 2) {
        stop("Dimensions of the similarities differ")
    }
    if (!all(sapply(sim, isSymmetric))) {
        stop("Similarities are not symmetric")
    }
    if (any(is.na(sim))) {
        warning("There are some NA values in the similarities provided")
    }
    FUN <- match.fun(func)
    # Apply weighted to each cell position of each similarity measure
    apply(simplify2array(sim), c(1,2), FUN, ...)
}

# addSimilarities ####
#' Additive integration of similarities
#'
#' Function that used the previously calculated similarities into a single
#' similarity matrix.
#'
#' The total weight sum can't be higher than 1 to prevent values above
#' 1 but can be below 1.
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
        stop("Expected a matrix, generally a similarity measure from ",
             "expression")
    }
    if (!all(dim(x) == dim(bio_mat[[1]]))) {
        stop("Dimensions of x and bio_mat matrices is different")
    }
    cors <- c(list(exp = x), bio_mat)

    # Apply weighted to each cell position of each similarity measure
    similarities(cors, weighted, w = weights)
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
    if (!all(weights >= 0)) {
        warning("There are negative weights, should they be positive?")
    }
    if (sum(weights) > 1) {
        warning("The sum of the weights is above 1")
    }

    sum(x*weights, na.rm = TRUE)
}

# Adding more genes to calculated similarities ####
#' Add a similarity measure for a gene to existing similarities
#'
#' Given a gene and a function to calculate the similarity with add it to the
#' previously calculated similarities in sim.
#' @param gene Entrez Identifier of the gene to add the similarity calculation
#' @param func Function to calculate the similarity
#' @param sim Matrix with similarities already calculated.
#' @param ... Other arguments passed to the function to calculate similarity
#' @return A matrix with the similarities with the new similarity
#' @export
gene2Sim <- function(gene, func, sim, ...) {
    FUNC <- match.fun(func)
    genes <- colnames(sim)
    dots <- list(...)
    new.sim <- sapply(genes, FUNC, gene, ... = dots)
    n.sim <- cbind(sim, new.sim)
    n.sim <- rbind(n.sim, c(new.sim, 1))
    colnames(n.sim) <- c(colnames(sim), gene)
    return(n.sim)

}

#' Add a similarity measure for a list of genes to existing similarities
#'
#' Given a gene and a function to calculate the similarity with add it to the
#' previously calculated similarities in sim.
#' @param gene A list of Entrez Identifier of the gene to add to the similarity
#' calculation
#' @param func Function to calculate the similarity
#' @param sim Matrix with similarities already calculated.
#' @param ... Other arguments passed to the functio to calculate the
#' similarities
#' @return A matrix with the similarities with the new similarity
mgene2Sim <- function(gene, func, sim, ...) {
    FUNC <- match.fun(func)
    dots <- list(...)
    for (g in gene) {
        m <- gene2Sim(g, FUNC, m, ... = dots)
    }
    return(m)
}

#' Parallel power
#'
#' A wrapper to registerDoParallel code
#' @param n Number of cores used
#' @return The result of registerDoParallel
#' @seealso \code{\link[doParallel]{registerDoParallel}}
#' @export
#' @import doParallel
setCores <- function(n = 2) {
    registerDoParallel(n)
}
