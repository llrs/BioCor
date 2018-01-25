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
#' Fills a matrix of \code{ncol = length(x)} and \code{nrow = length(x)} with
#' the values in \code{dat} and setting the diagonal to 1.
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
    out[upper.tri(out)] <- unlist(dat, use.names = TRUE)
    out[lower.tri(out)] <- t(out)[lower.tri(t(out))]
    diag(out) <- 1L
    rownames(out) <- colnames(out) <- x
    return(out)
}

# duplicateIndices ####
#' Finds the indices of the duplicated events of a vector
#'
#' Finds the indices of duplicated elements in the vector given.
#'
#' For each duplication it can return a list or if all the duplication events
#' are of the same length it returns a matrix, where each column is duplicated.
#' @param vec Vector of identifiers presumably duplicated
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

# removeDup ####
#' Remove duplicated rows and columns
#'
#' Given the indices of the duplicated entries remove the columns and rows
#' until just one is left, it keeps the duplicated with the highest absolute
#' mean value.
#'
#' @param cor_mat List of matrices
#' @param dupli List of indices with duplicated entries
#' @return A matrix with only one of the columns and rows duplicated
#' @export
#' @author Lluís Revilla
#' @seealso \code{\link{duplicateIndices}} to obtain the list of indices with
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

