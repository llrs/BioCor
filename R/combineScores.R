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
#' @note This is a version of combineScores from
#' \code{\link[GOSemSim]{combineScores}} with optional rounding and some
#' internal differences.
#' @export
#' @author Lluís Revilla based on Guangchuang Yu
#' @references
#' Ying Tao, Lee Sam, Jianrong Li, Carol Friedman, Yves A. Lussier;
#' Information theory applied to the sparse gene ontology annotation network to
#' predict novel gene function. Bioinformatics 2007; 23 (13): i529-i538.
#' doi: 10.1093/bioinformatics/btm195
#' @examples
#' d <- structure(c(0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
#' 0.285714285714286), .Dim = c(3L, 3L), .Dimnames = list(c("a",
#' "b", "c"), c("d", "e", "f")))
#' d
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
    if (!is.matrix(scores)) {
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
    if (any(is.na(scores)) && is.matrix(scores)) {
        row.na.idx <- apply(scores, 1, function(i){all(is.na(i))})
        if (any(row.na.idx)) {
            scores <- scores[-which(row.na.idx), ]
        }

    }
    if (any(is.na(scores)) && is.matrix(scores)) {
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
    } else if (is.matrix(scores)) {
        if (method == "rcmax") {
            result <-  max(rowMeans(scores, na.rm = TRUE),
                           colMeans(scores, na.rm = TRUE))
        } else if (method == "rcmax.avg" || method == "BMA") {
            result <- sum(apply(scores, 1, max, na.rm = TRUE),
                          apply(scores, 2, max, na.rm = TRUE)) / sum(
                              dim(scores))
        } else if (method == "reciprocal") {
            rows <- apply(scores, 2, which.max)
            columns <- apply(scores, 1, which.max)
            reciprocal <- rows[columns[rows] == seq_along(rows)]
            rScores <- sapply(seq_along(reciprocal), function(x) {
                scores[reciprocal[x], names(reciprocal)[x]]
            })

            if (all(rScores <= t)) {
                result <- NA
            } else {
                result <- 2*sum(rScores[rScores >= t])/(sum(dim(scores)))
            }

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
