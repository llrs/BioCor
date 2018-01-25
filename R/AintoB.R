#' Insert a matrix into another
#'
#' Insert values from a matrix into another matrix based on the rownames and
#' colnames replacing the values.
#'
#' If all the genes with pathway information are already calculated but you
#' would like to use more genes when performing analysis. insert the once you
#' have calculated on the matrix of genes.
#' @param A A matrix to be inserted.
#' @param B A matrix to insert in.
#' @return A matrix with the values of A in the matrix B.
#' @author Lluís Revilla
#' @examples
#' B <- matrix(ncol = 10, nrow = 10,
#'     dimnames = list(letters[1:10], letters[1:10]))
#' A <- matrix(c(1:15), byrow=TRUE, nrow=5,
#'     dimnames = list(letters[1:5], letters[1:3]))
#' AintoB(A, B)
#'
#' # Mixed orders
#' colnames(A) <- c("c", "h", "e")
#' rownames(A) <- c("b", "a", "f", "c", "j")
#' AintoB(A, B)
#'
#' # Missing colums or rows
#' colnames(A) <- c("d", "f", "k")
#' AintoB(A, B)
#' @export
AintoB <- function(A, B) {
    if (!is.matrix(A) | !is.matrix(B)) {
        stop("Input should be matrices")
    }
    # Select the order for columns
    if (ncol(A) <= ncol(B) & nrow(A) <= nrow(B)) {
        mc <- match(colnames(A), colnames(B))
        mr <- match(rownames(A), rownames(B))
    } else {
        stop("Impossible to insert matrix A into matrix B\n",
            "Matrix A is bigger than matrix B.")
    }
    # Omit those with NA
    nar <- is.na(mr)
    nac <- is.na(mc)
    mr <- mr[!nar]
    mc <- mc[!nac]

    B[mr, mc] <- A[!nar, !nac]
    B
}
