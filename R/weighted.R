# weighted ####
#' Weighted operations
#'
#' Calculates the weighted sum or product of \code{x}. Each values should have
#' its weight, otherwise it will throw an error.
#'
#' This functions are thought to be used with \code{similarities}. As some
#' similarities might be positive and others negative the argument \code{abs}
#' is provided for \code{weighted.sum}, assuming that only one similarity will
#' be negative (usually the one coming from expression correlation).
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
#' \code{\link{weighted.mean}}, \code{\link{similarities}} and
#' \code{\link{addSimilarities}}
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
        stop(
            "Weights and data don't match the length: ", length(x), " != ",
            length(w)
        )
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
        if (any(sign(x[!is.na(x)]) < 0)) {
            -sum(abs(x) * w, na.rm = TRUE)
        } else {
            sum(x * w, na.rm = TRUE)
        }
    } else {
        sum(x * w, na.rm = TRUE)
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
        stop(
            "Weights and data don't match the length: ", length(x), " != ",
            length(w)
        )
    }

    if (!is.numeric(x) | !is.numeric(w)) {
        stop("weights and x should be numeric")
    }

    if (sum(w, na.rm = TRUE) > 1L) {
        warning("The sum of the weights is above 1")
    }

    prod(x * w, na.rm = TRUE)
}
