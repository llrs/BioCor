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
#' @note It doesn't check that the columns and rows of the matrices are in the
#' same order or are the same.
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
    if (!is.list(sim)) {
        stop("Please introduce a list with the similarities")
    }
    # Check that all the matrices are of the same dimensions and squared
    dims <- sapply(sim, dim)
    if (is.list(dims)) {
        stop("There are arrays with different dimensions")
    } else if (!all(apply(dims, 2, function(x) length(unique(x)) == 1))) {
        stop("Matrices with different dimensions.")
    }
    if (!all(sapply(sim, isSymmetric))) {
        stop("Similarities are not symmetric")
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
#' 1 but can be below 1. It uses weighted.sum with abs = TRUE internally.
#' @param x A matrix with the similarity of expression
#' @param bio_mat A list of matrices of the same dimension as x.
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
addSimilarities <- function(x, bio_mat, weights = c(0.5, 0.18, 0.10, 0.22)) {
    if (sum(weights) > 1L) {
        stop("Weights are too big. The sum must be equal to 1")
    } else if (sum(weights) < 1L) {
        warning("Weights are smaller than 1.")
    }

    if (!is.matrix(x)) {
        stop(
            "Expected a matrix, generally a similarity measure from ",
            "expression"
        )
    }
    if (!all(dim(x) == dim(bio_mat[[1L]]))) {
        stop("Dimensions of x and bio_mat matrices is different")
    }
    cors <- c(list(exp = x), bio_mat)

    # Apply weighted to each cell position of each similarity measure
    similarities(cors, weighted.sum, w = weights)
}
