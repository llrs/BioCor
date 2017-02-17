# combinadic ####
#' i-th combination of n elements taken from r
#'
#' Function similar to combn but for larger vectors. To avoid allocating a big
#' vector with all the combinations each one can be computed with this
#' function.
#' @param n Elements to extract the combination from
#' @param r Number of elements per combination
#' @param i ith combination
#' @return The combination ith of the elements
#' @rdname combinadic
#' @seealso \code{\link{combn}}
#' @export
#' @examples
#' #Output of all combinations
#' combn(LETTERS[1:5], 2)
#' # Otuput of the second combination
#' combinadic(LETTERS[1:5], 2, 2)
#' @author Joshua Ulrich
#' @references
#' \href{http://stackoverflow.com/a/4494469/2886003}{StackOverflow answer
#' 4494469/2886003}
combinadic <- function(n, r, i) {

    # http://msdn.microsoft.com/en-us/library/aa289166(VS.71).aspx
    # http://en.wikipedia.org/wiki/Combinadic
    n0 <- length(n)
    if (i < 1L | i > choose(n0, r)) {
        stop("'i' must be 0 < i <= n0!/(n0-r)!")
    }
    largestV <- function(n, r, i) {
        v <- n # Adjusted for one-based indexing
        while (choose(v,r) >= i) { # Adjusted for one-based indexing
            v <- v - 1L
        }
        return(v)
    }

    res <- rep(NA,r)
    for (j in 1L:r) {
        res[j] <- largestV(n0, r, i)
        i <- i - choose(res[j], r)
        n0 <- res[j]
        r <- r - 1L
    }
    res <- res + 1L
    res <- n[res]
    return(res)
}

# seq2mat ####
#' Transforms a vector to a symmetric matrix
#'
#' Fils a matrix of \code{ncol = length(x)} and \code{nrow = length(x)} with
#' the values in \code{dat} and setting the diagnoal to 1.
#'
#' \code{dat} should be at least \code{choose(length(x), 2)} of length. It
#' assumes that the data provided comes from using the row and column id to
#' obtain it.
#' @param x names of columns and rows, used to define the size of the matrix
#' @param dat Data to fill with the matrix with except the diagonal.
#' @return A square matrix with the diagonal set to 1 and \code{dat} on the
#' upper and lower triangle with the columns ids and row ids from x.
#' @examples
#' seq2mat(LETTERS[1:5], 1:10)
#' seq2mat(LETTERS[1:5], seq(from = 0.1, to = 1, by = 0.1))
#' @export
#' @seealso \code{\link{upper.tri}} and \code{\link{lower.tri}}
#' @author Lluís Revilla
seq2mat <- function(x, dat) {
    if (length(dat) != choose(length(x), 2L)) {
        stop("Data is not enough big to populate the matrix")
    }
    out <- matrix(ncol = length(x), nrow = length(x))
    out[upper.tri(out)] <- unlist(dat)
    out[lower.tri(out)] <- t(out)[lower.tri(t(out))]
    diag(out) <- 1L
    rownames(out) <- colnames(out) <- x
    return(out)
}
# similarities ####
#' Apply a function to a list of similarities
#'
#' Function to join list of similarities by a function provided by the user.
#' @param sim list of similarities to be joined. All similarities must have the
#'  same dimensions. The genes are assumed to be in the same order for all the
#'  matrices.
#' @param func function to perform on those similarities: \code{prod},
#' \code{sum}... It should accept as many arguments as similarities matrices
#' are provided, and should use numbers.
#' @param ... Other arguments passed to the function \code{func}. Usually na.rm
#' or similar.
#' @return A matrix of the size of the similarities
#' @export
#' @seealso \code{\link{weighted}} for functions that can be used, and
#' \code{\link{addSimilarities}} for a wrapper to one of them
#' @author Lluís Revilla
#' @examples
#' set.seed(100)
#' a <- seq2mat(LETTERS[1:5], rnorm(10))
#' b <- seq2mat(LETTERS[1:5], seq(from = 0.1, to = 1, by = 0.1))
#' sim <- list(b, a)
#' similarities(sim, weighted.prod, c(0.5, 0.5))
#' # Note the differences in the sign of some values
#' similarities(sim, weighted.sum, c(0.5, 0.5))
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
    apply(simplify2array(sim), c(1L, 2L), FUN, ...)
}

# addSimilarities ####
#' Additive integration of similarities
#'
#' Function that use the previously calculated similarities into a single
#' similarity matrix.
#'
#' The total weight can't be higher than 1 to prevent values above
#' 1 but can be below 1.
#' @param x A matrix with the similarity of expression
#' @param bio_mat A list of matrices of the same dimension
#' @param weights A numeric vector of weight to multiply each similarity
#' @return A square matrix of the same dimensions as the input matrices.
#' @export
#' @seealso \code{\link{similarities}}, \code{\link{weighted}}.
#' @author Lluís Revilla
#' @examples
#' set.seed(100)
#' a <- seq2mat(LETTERS[1:5], rnorm(10))
#' b <- seq2mat(LETTERS[1:5], seq(from = 0.1, to = 1, by = 0.1))
#' sim <- list(b)
#' addSimilarities(a, sim, c(0.5, 0.5))
addSimilarities <- function(x, bio_mat, weights = c(0.5, 0.18, 0.10, 0.22)){
    # exp, reactome, kegg, go
    # cor_mat <- cor(x, use = "p")
    if (sum(weights) > 1L) {
        stop("Weights are too big. The sum must be equal to 1")
    } else if (sum(weights) < 1L) {
        warning("Weights are smaller than 1.")
    }

    if (!is.matrix(x)) {
        stop("Expected a matrix, generally a similarity measure from ",
             "expression")
    }
    if (!all(dim(x) == dim(bio_mat[[1L]]))) {
        stop("Dimensions of x and bio_mat matrices is different")
    }
    cors <- c(list(exp = x), bio_mat)

    # Apply weighted to each cell position of each similarity measure
    similarities(cors, weighted.sum, w = weights)
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
#' @author Lluís Revilla
#' @seealso \code{\link{removeDup}}
#' @examples
#' duplicateIndices(c("52", "52", "53", "55")) # One repeated element
#' duplicateIndices(c("52", "52", "53", "55", "55")) # Repeated elements
#' duplicateIndices(c("52", "55", "53", "55", "52")) # Mixed repeated elements
duplicateIndices <- function(vec) {
    if (!is.character(vec)) {
        stop("Expected a list of characters to find duplicates on it")
    }
    sapply(unique(vec[duplicated(vec)]), function(x){
        b <- 1:length(vec)
        b[vec == x]}, simplify = FALSE)
}

# weighted ####
#' Weighted operations
#'
#' Calculates the weighted sum or product of \code{x}. Each values should have
#' its weight, otherwise it will throw an error.
#'
#' This functions are thought to be used with \code{similarities}. As some
#' similarities might be positive and others negative the argument \code{abs} is
#' provided for \code{weighted.sum}, assuming that only one similarity will be
#' negative (usually the one comming from exprresion correlation).
#' @inheritParams stats::weighted.mean
#' @param x an object containing the values whose weighted operations is to be
#' computed
#' @param abs If any \code{x} is negative you want the result negative too?
#' @return \code{weighted.sum} returns the sum of the product of x*weights
#' removing all \code{NA} values. See parameter \code{abs} if there are any
#' negative values.
#' @aliases weighted
#' @rdname weighted
#' @name weighted
#' @aliases weighted
#' @author Lluís Revilla
#' @seealso
#' \code{\link{similarities}} and \code{\link{addSimilarities}}
#' @export
#' @examples
#' expr <- c(-0.2, 0.3, 0.5, 0.8, 0.1)
#' weighted.sum(expr, c(0.5, 0.2, 0.1, 0.1, 0.1))
#' weighted.sum(expr, c(0.5, 0.2, 0.1, 0.2, 0.1), FALSE)
#' weighted.sum(expr, c(0.4, 0.2, 0.1, 0.2, 0.1))
#' weighted.sum(expr, c(0.4, 0.2, 0.1, 0.2, 0.1), FALSE)
#' weighted.sum(expr, c(0.4, 0.2, 0, 0.2, 0.1))
#' weighted.sum(expr, c(0.5, 0.2, 0, 0.2, 0.1))
weighted.sum <- function(x, w, abs = TRUE) {
    if (length(x) != length(w)) {
        stop("Weights and data don't match the length: ", length(x), " != ",
             length(w))
    }

    if (!is.numeric(x) | !is.numeric(w)) {
        stop("weights and x should be numeric")
    }
    if (!all(w >= 0L)) {
        warning("There are negative weights, should they be positive?")
    }
    if (sum(w, na.rm = TRUE) > 1L) {
        warning("The sum of the weights is above 1")
    }
    if (abs) {
        if (any(sign(x) < 0)) {
            -sum(abs(x)*w, na.rm = TRUE)
        } else {
            sum(x*w, na.rm = TRUE)
        }
    } else {
        sum(x*w, na.rm = TRUE)
    }
}

#' @return \code{weighted.prod} returns the product of product of x*weights
#'  removing all \code{NA} values.
#' @rdname weighted
#' @export
#' @examples
#' # Compared to weighted.prod:
#' weighted.prod(expr, c(0.5, 0.2, 0.1, 0.1, 0.1))
#' weighted.prod(expr, c(0.4, 0.2, 0.1, 0.2, 0.1))
#' weighted.prod(expr, c(0.4, 0.2, 0, 0.2, 0.1))
#' weighted.prod(expr, c(0.5, 0.2, 0, 0.2, 0.1))
weighted.prod <- function(x, w) {
    if (length(x) != length(w)) {
        stop("Weights and data don't match the length: ", length(x), " != ",
             length(w))
    }

    if (!is.numeric(x) | !is.numeric(w)) {
        stop("weights and x should be numeric")
    }

    if (sum(w, na.rm = TRUE) > 1L) {
        warning("The sum of the weights is above 1")
    }

    prod(x*w, na.rm = TRUE)
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

# Conversions ####
#' Convert the similarities formats
#'
#' Functions to convert the similarity coefficients between Jaccard and Dice.
#' D2J is the opposite of J2D.
#' @param D Dice coefficient, as returned by pathSim, genesSim and bioCor
#' @param J Jaccard coefficient
#' @return A numeric value.
#' @author Lluís Revilla
#' @export
#' @rdname conversions
#' @name conversions
#' @examples
#' D2J(0.5)
#' J2D(0.5)
#' D2J(J2D(0.5))
D2J <- function(D) {
    if (all(D > 1)) {
        stop("Dice index can't be above 1")
    } else if (all(D < 0)) {
        stop("Dice index can't be below 0")
    }
    D/(2 - D)
}

#' @export
#' @rdname conversions
J2D <- function(J) {
    if (all(J > 1)) {
        stop("Jaccard index can't be above 1")
    } else if (all(J < 0)) {
        stop("Jaccard index can't be below 0")
    }
    2*J/(1 + J)
}

# combineScores ####
#' Combining values
#'
#' Combine several values into one by several methods.
#'
#' The methods return:
#' \describe{\item{avg}{The average or mean value}
#' \item{max}{The max value}
#' \item{rcmax}{The max of the column means or row means}
#' \item{rcmax.avg}{The sum of the max values by rows and columns divided by
#' the number of columns and rows}
#' \item{BMA}{The same as \code{rcmax.avg}}}
#' @param scores Scores to be combined, some methods accept a matrix and a
#' vector others only a matrix.
#' @param method one of \code{c("avg", "max", "rcmax", "rcmax.avg", "BMA")}
#' @param round Should the resulting value be rounded to the third digit?
#' @return A numeric value as described in details.
#' @note This is a version of combineScores inspired from GOSemSim
#' \code{\link[GOSemSim]{combineScores}} with optional rounding and some
#' internal differences.
#' @export
#' @examples
#' d <- structure(c(0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
#' 0.285714285714286), .Dim = c(3L, 3L), .Dimnames = list(c("a",
#' "b", "c"), c("d", "e", "f")))
#' d
#' sapply(c("avg", "max", "rcmax", "rcmax.avg", "BMA"), combineScores,
#'        scores = d)
combineScores <- function(scores, method, round = FALSE) {
    # Check input
    method <- match.arg(method, c("avg", "max", "rcmax", "rcmax.avg", "BMA"))
    if (is.list(scores) & length(scores) == 1) {
        scores <- scores[[1]]
    }
    if (!is.matrix(scores)) {
        stop("Please introduce a matrix")
    } else if (any(dim(scores) == 0L)) {
        return(NA)
    }

    # Remove NA
    if (length(dim(scores)) == 2) {
        row.na.idx <- apply(scores, 1, function(i){all(is.na(i))})
        col.na.idx <- apply(scores, 2, function(i){all(is.na(i))})
        if (any(row.na.idx)) {
            scores <- scores[-which(row.na.idx), ]
        }
        if (any(col.na.idx)) {
            scores <- scores[, -which(col.na.idx)]
        }
    }

    # Apply the methods
    if (method == "avg") {
        result <- mean(scores, na.rm = TRUE)
    } else if (method == "max") {
        result <- max(scores, na.rm = TRUE)
    } else if (method == "rcmax") {
        result <-  max(rowMeans(scores, na.rm = TRUE),
                       colMeans(scores, na.rm = TRUE))
    } else if (method == "rcmax.avg" || method == "BMA") {
        result <- sum(apply(scores, 1, max, na.rm = TRUE),
                      apply(scores, 2, max, na.rm = TRUE)) / sum(dim(scores))
    }

    # Return the value
    if (round) {
        return(round(result, digits = 3))
    } else {
        return(result)
    }

}
