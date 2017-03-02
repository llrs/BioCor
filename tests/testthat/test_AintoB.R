library("BioCor")
context("Testing AintoB functions")

test_that("AintoB", {
    B <- matrix(ncol = 10, nrow = 10,
                dimnames = list(letters[1:10], letters[1:10]))
    A <- matrix(c(1:15), byrow = T, nrow = 5,
                dimnames = list(letters[1:5], letters[1:3]))
    test <- AintoB(A, B)
    expect_equal(test[1:5, 1:3], A)

    # Mixed orders
    colnames(A) <- c("c", "h", "e")
    rownames(A) <- c("b", "a", "f", "c", "j")
    test <- AintoB(A, B)
    expect_equal(test[1, 3], 4L)
    expect_equal(colnames(test)[4L], "d")

    # Missing colums or rows
    colnames(A) <- c("d", "f", "k")
    test <- AintoB(A, B)
    expect_equal(test[2L, 4L], 1L)
})
