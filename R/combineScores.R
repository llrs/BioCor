# combineScores ####
#' Combining values
#'
#' Combine several similarities into one using several methods.
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
combineScores <- function(scores, method = c("max", "avg", "rcmax", "rcmax.avg", "BMA",
                                             "reciprocal"), round = FALSE, t = 0) {
    # Check input
    method <- match.arg(method,
                        c("avg", "max", "rcmax", "rcmax.avg", "BMA",
                          "reciprocal"))


    if (is.list(scores) && length(scores) == 1) {
        scores <- scores[[1]]
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

    if (!is.numeric(t) || t < 0 || t > 1) {
        stop("t is for the reciprocal method, and should be between 0 and 1")
    }

    scores <- removeNA(scores)

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
    rounder(result, round)
}

rounder <- function(result, round) {
    if (round) {
        return(round(result, digits = 3))
    } else {
        return(result)
    }
}
removeNA <- function(m){
    # Remove NA
    if (any(is.na(m)) && is.Matrix(m)) {
        row.na.idx <- apply(m, 1, function(i){all(is.na(i))})
        col.na.idx <- apply(m, 2, function(i){all(is.na(i))})
        if (any(row.na.idx)) {
            m <- m[-which(row.na.idx), , drop = FALSE]
        }
        if (any(col.na.idx)) {
            m <- m[, -which(col.na.idx)]
        }
    }
    m
}

# Check Matrix classes and matrix
#' @importFrom methods is
is.Matrix <- function(m){
    is.matrix(m) || is(m, "Matrix")
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
    # In case of a draw select the first 4
    # (it wouldn't matter as it should be all 1)
    rScores <- rScores[seq_len(min(ncol(scores), nrow(scores)))]
    rScores <- rScores[!is.na(rScores)]

    if (all(rScores < t)) {
        NA
    } else {
        2*sum(rScores[rScores >= t])/(sum(dim(scores)))
    }
}

#' \code{cobmineScoresPar} performs multiple \code{combineScores}. It allows
#' parallel computing using BiocParallel back-end.
#' @param subSets List of combinations as info in other functions.
#' @param BPPARAM BiocParallel back-end parameters.
#' By default (\code{NULL}) a \code{for} loop is used.
#' @param ... Other arguments passed to \code{combineScores}
#' @seealso \link[BiocParallel]{register} in BiocParallel about the arguments
#' accepted by BPPARAM
#' @rdname combineScores
#' @export
#' @import Matrix
#' @import BiocParallel
#' @importFrom utils combn
#' @examples
#' colnames(e) <- rownames(e)
#' combineScoresPar(e, list(a= c("a", "b"), b = c("b", "c")),
#'                  method = "max")
combineScoresPar <- function(scores,
                             method,
                             subSets = NULL,
                             BPPARAM = NULL,
                             ...){

    # Check scores
    if (is.null(subSets) | sum(dim(scores)) == 0) {
        return(combineScores(scores, method = method, ...))
        # To handle cases when it is already simplified
    } else if (length(subSets) == 1) {
        return(combineScores(scores, method = method, ...))
    } else {
        # To handle cases where subSets are not present in scores
        B <- matrix(NA, ncol = length(subSets), nrow = length(subSets),
                    dimnames = list(names(subSets), names(subSets)))
        subSets <- keepSubSet(subSets, scores)
        if (length(subSets) == 0) {
            return(B)
        }

    }

    #all combinations of indices
    ij <- combn(seq_along(subSets), 2) # Use fun

    #add all i = j to the combination of indices
    ij <- matrix(c(ij, rep(seq_along(subSets), each = 2)), nrow = 2)

    if (is.null(BPPARAM)) { # If not a parallel background is provided
        # only one loop
        res <- numeric(ncol(ij)) # preallocate
        for (k in seq_len(ncol(ij))) {
            rowIds <- subSets[[ij[1, k]]]
            colIds <- subSets[[ij[2, k]]]
            if (anyNA(c(rowIds, colIds)) || is.null(colIds) || is.null(rowIds)) {
                res[k] <- NA
            } else {
                res[k] <- combineScores(scores[rowIds, colIds, drop = FALSE],
                                        method, ... = ...)
            }
        }
    } else {
        # Use the parallel background provided
        res <- bplapply(seq_len(ncol(ij)), function(x){
            if (ncol(ij) < x) {
                message(print(ij))
            }
            rowIds <- subSets[[ij[1, x]]]
            colIds <- subSets[[ij[2, x]]]
            if (anyNA(c(rowIds, colIds)) || is.null(rowIds) || is.null(colIds)) {
                NA
            } else {
                combineScores(scores[rowIds, colIds, drop = FALSE],
                                        method, ... = ...)
            }
        }, BPPARAM = BPPARAM)
        res <- as.numeric(res)
    }

    A <- sparseMatrix(i = ij[1,], j = ij[2,],
                      x = res, dims = rep(length(subSets), 2),
                      symmetric = TRUE, index1 = TRUE,
                      dimnames = list(names(subSets), names(subSets)))
    AintoB(as.matrix(A), B)
}


keepSubSet <- function(subSets, scores) {
    # Check the ids
    subId <- unique(unlist(subSets, use.names = FALSE))
    values <- unlist(dimnames(scores), use.names = FALSE)

        # browser()
    if (anyNA(subId)) {
        keep <- sapply(subSets, check_in, values = values)
        subSets[keep]
    } else {
        subSets
    }
}

check_in <- function(x, values) {
    if (all(is.na(x) | is.null(x))) {
        FALSE
    } else {
        all(x %in% values)
    }
}

