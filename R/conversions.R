# Conversions ####
#' Convert the similarities formats
#'
#' Functions to convert the similarity coefficients between Jaccard and Dice.
#' D2J is the opposite of J2D.
#' @param D Dice coefficient, as returned by [diceSim()],
#' [geneSim()], [clusterSim()] and
#' [clusterGeneSim()]
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
    D / (2 - D)
}

#' @export
#' @rdname conversions
J2D <- function(J) {
    if (all(J > 1)) {
        stop("Jaccard index can't be above 1")
    } else if (all(J < 0)) {
        stop("Jaccard index can't be below 0")
    }
    2 * J / (1 + J)
}
