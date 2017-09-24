# combineScores ####
#' Combining values
#'
#' Combine several values into one by several methods.
#'
#' The input matrix can be a base matrix or a matrix from package Matrix.
#' The methods return:
#' \describe{\item{avg}{The average or mean value}
#' \item{max}{The max value}
#' \item{rcmax}{The max of the column means or row means}
#' \item{rcmax.avg}{The sum of the max values by rows and columns divided by
#' the number of columns and rows}
#' \item{BMA}{The same as \code{rcmax.avg}}
#' \item{reciprocal}{The double of the sum of the reciprocal maximal
#' similarities (above  a threshold) divided by the number of elements.
#' See equation 3 of the Tao \emph{et al} 2007 article}}
#' @param scores Matrix of scores to be combined
#' @param method one of \code{c("avg", "max", "rcmax", "rcmax.avg", "BMA",
#' "reciprocal")}, see Details.
#' @param round Should the resulting value be rounded to the third digit?
#' @param t Numeric value to filter scores below this value. Only used in the
#' reciprocal method.
#' @return A numeric value as described in details.
#' @note \code{combineScores} is a version of the function of the same name in
#' package GOSemSim (\code{\link[GOSemSim]{combineScores}}) with optional
#' rounding and some internal differences.
#' @export
#' @author Lluís Revilla based on Guangchuang Yu
#' @references
#' Ying Tao, Lee Sam, Jianrong Li, Carol Friedman, Yves A. Lussier;
#' Information theory applied to the sparse gene ontology annotation network to
#' predict novel gene function. Bioinformatics 2007; 23 (13): i529-i538.
#' doi: 10.1093/bioinformatics/btm195
#' @examples
#' (d <- structure(c(0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
#'                   0.285714285714286), .Dim = c(3L, 3L),
#'                 .Dimnames = list(c("a", "b", "c"), c("d", "e", "f"))))
#' e <- d
#' sapply(c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal"),
#'        combineScores, scores = d)
#' d[1,2] <- NA
#' sapply(c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal"),
#'        combineScores, scores = d)
combineScores <- function(scores, method, round = FALSE, t = 0) {
    # Check input
    method <- match.arg(method,
                        c("avg", "max", "rcmax", "rcmax.avg", "BMA",
                          "reciprocal"))


    if (is.list(scores) && length(scores) == 1) {
        scores <- scores[[1]]
    }
    # Check Matrix classes and matrix
    is.Matrix <- function(m){
        is.matrix(m) | grepl("Matrix", as.character(class(m)),
                                   ignore.case = FALSE)
    }
    if (!is.Matrix(scores)) {
        stop("scores argument should be a matrix")
    }
    if (length(round) != 1 || !is.logical(round)) {
        stop("round argument is not logical")
    }

    if (any(dim(scores) == 0L)) {
        return(NA)
    }

    if (!is.numeric(t) | t < 0 | t > 1) {
        stop("t is for the reciprocal method, and should be between 0 and 1")
    }

    # Remove NA
    if (any(is.na(scores)) && is.Matrix(scores)) {
        row.na.idx <- apply(scores, 1, function(i){all(is.na(i))})
        if (any(row.na.idx)) {
            scores <- scores[-which(row.na.idx), ]
        }

    }
    if (any(is.na(scores)) && is.Matrix(scores)) {
        col.na.idx <- apply(scores, 2, function(i){all(is.na(i))})
        if (any(col.na.idx)) {
            scores <- scores[, -which(col.na.idx)]
        }
    }

    # Apply the methods
    if (method == "avg") {
        result <- mean(scores, na.rm = TRUE)
    } else if (method == "max") {
        result <- max(scores, na.rm = TRUE)
    } else if (is.Matrix(scores)) {
        if (method == "rcmax") {
            result <- rcmax(scores)
        } else if (method == "rcmax.avg" || method == "BMA") {
            result <- BMA(scores)
        } else if (method == "reciprocal") {
            result <- reciprocal(scores, t)
        }
    } else {
        warning("Using max method because after removing NAs ",
                "there isn't a matrix to use.")
        result <- max(scores, na.rm = TRUE)
    }

    # Return the value
    if (round) {
        return(round(result, digits = 3))
    } else {
        return(result)
    }

}

BMA <- function(scores) {
    if (length(dim(scores)) != 2) {
        stop("scores must be a matrix.")
    }
    sum(apply(scores, 1, max, na.rm = TRUE),
                  apply(scores, 2, max, na.rm = TRUE)) / sum(
                      dim(scores))
}
rcmax <- function(scores) {
    if (length(dim(scores)) != 2) {
        stop("scores must be a matrix.")
    }
    col <- mean(apply(scores, 1, max, na.rm = TRUE))
    row <- mean(apply(scores, 2, max, na.rm = TRUE))
    max(col, row, na.rm = TRUE)
}

reciprocal <- function(scores, t) {

    if (t < 0 | t > 1) {
        stop("t must be between 1 and 0")
    }
    if (length(dim(scores)) != 2) {
        stop("scores must be a matrix.")
    }
    # Find the ones that max row is also the max col
    fmax <- function(x){
        x == max(x, na.rm = TRUE)
    }
    # Max in rows
    # Transpose to correct for the fact that apply fill by column
    rowsMax <- which(t(apply(scores, 1, fmax)))
    # Max in columns
    colsMax <- which(apply(scores, 2, fmax))
    rScores <- scores[intersect(rowsMax, colsMax)]

    if (all(rScores < t)) {
        NA
    } else {
        2*sum(rScores[rScores >= t])/(sum(dim(scores)))
    }
}

combineScoresPrep <- function(x, y, prep, method, ...) {
    if (is.null(x) || is.null(y)) {
        NA
    } else if (all(is.na(x)) || all(is.na(y))) {
        NA
    } else {
        combineScores(prep[x[!is.na(x)], y[!is.na(y)], drop = FALSE],
                      method, ...)
    }
}
vcombineScoresPrep <- Vectorize(combineScoresPrep,
                                vectorize.args = c("x", "y"))

#' \code{cobmineScoresPar} performs multiple (parallel) combineScores based on
#' a list of elements to combine into one.
#' @param subSets List of combinations as info in other functions.
#' @param BPPARAM Determining the parallel back-end. By default a for loop is
#' used.
#' @param ... Other arguments passed to \code{\link{combineScoresPar}}
#' @seealso \code{\link[BiocParallel]{bpparam}}
#' @rdname combineScores
#' @export
#' @import Matrix
#' @import BiocParallel
#' @importFrom utils combn
#' @examples
#' colnames(e) <- rownames(e)
#' combineScoresPar(e, list(a= c("a", "b"), b = c("b", "c")),
#'                  method = "max")
combineScoresPar <- function(scores, method, subSets = NULL, BPPARAM = NULL, ...){

    # Check scores
    if (is.null(subSets) | sum(dim(scores)) == 0) {
        return(combineScores(scores, method = method, ...))
    } else { # To handle cases where subSets are not present in scores
        # Check the ids
        subId <- unique(unlist(subSets))
        if (!all(subId %in% unlist(dimnames(scores)))) {
            B <- matrix(NA, ncol = length(subSets), nrow = length(subSets),
                        dimnames = list(names(subSets), names(subSets)))
        }
        keep <- sapply(subSets, function(x) {
            all(x %in% unlist(dimnames(scores)))
            })
        subSets <- subSets[keep]
        if (length(subSets) == 0) {
            return(B)
        }
    }

    #all combinations of indices
    ij <- combn(seq_along(subSets), 2)

    #add all i = j to the combination of indices
    ij <- matrix(c(ij, rep(seq_along(subSets), each = 2)), nrow = 2)

    if (is.null(BPPARAM)) { # If not a parallel background is provided
        # only one loop
        res <- numeric(ncol(ij)) # preallocate
        for (k in seq_len(ncol(ij))) {
            rowIds <- subSets[[ij[1, k]]]
            colIds <- subSets[[ij[2, k]]]
            res[k] <- combineScores(scores[rowIds, colIds, drop = FALSE], method,
                ... = ...)
        }
    } else {
        # Use the parallel background provided
        res <- bplapply(seq_len(ncol(ij)), function(x){
            if (ncol(ij) < x) {
                message(print(ij))
            }
            rowIds <- subSets[[ij[1, x]]]
            colIds <- subSets[[ij[2, x]]]
            combineScores(scores[rowIds, colIds, drop = FALSE], method, ... = ...)
        }, BPPARAM = BPPARAM)
        res <- as.numeric(res)
    }

    A <- sparseMatrix(i = ij[1,], j = ij[2,],
                      x = res, dims = rep(length(subSets), 2),
                      symmetric = TRUE, index1 = TRUE,
                      dimnames = list(names(subSets), names(subSets)))
    if (exists("B")) {
        return(AintoB(as.matrix(A), B))
    } else {
        return(A)
    }
}
