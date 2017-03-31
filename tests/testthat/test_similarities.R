library("BioCor")
context("Testing Similarities function")

set.seed(1)
a <- seq2mat(LETTERS[1:5], runif(10))
b <- seq2mat(LETTERS[1:5], runif(10))
x <- seq2mat(LETTERS[1:5], runif(10))
mat <- list("kegg" = a , "react" = b)

test_that("addSimilarities", {
    weights <- c(0.1, 0.2, 0.8)
    expect_error(addSimilarities(x, mat, weights), "Weights are too big")
    weights <- c(0.1, 0.2, 0.3)
    expect_warning(addSimilarities(x, mat, weights), "smaller than 1")
    weights <- c(0.3, 0.3, 0.4)
    expect_error(addSimilarities("a", mat, weights), "Expected a matrix")
    test <- addSimilarities(x, mat, weights)

    expect_true(isSymmetric(test))
    expect_true(all(diag(test) == 1))
    expect_true(all(test <= 1))
    expect_true(all(test >= -1))
    expect_true(test["B", "D"] <= 0.6129278)

    x <- matrix(runif(10), ncol = 5)
    expect_error(addSimilarities(x, mat, weights), "is different")
})

test_that("similarities", {

    test <- similarities(mat, mean)
    expect_error(similarities(matrix(, ncol = 2), mean), "introduce")
    expect_error(similarities(list(mat = matrix(, ncol = 2, nrow = 2),
                                   met = NA), mean), "differ")
    expect_error(similarities(list(mat = matrix(, ncol = 2, nrow = 2),
                                   met = matrix(ncol = 2, nrow = 3)), mean),
                 "Dimensions")

    expect_equal(test["A", "B"], (b["A", "B"] + a["A", "B"])/2)
    expect_equal(test["A", "B"], 0.2357416190207)

    test <- similarities(mat, weighted.sum, w = c(0.8, 0.1))
    expect_equal(test["A", "B"], 0.23300438800361)

    test <- similarities(mat, prod)
    expect_equal(test["A", "B"], 0.0546880340227757)

    mat2 <- mat
    mat2$react[1, 2] <- 1
    expect_error(similarities(mat2, prod), "symmetric")
})
