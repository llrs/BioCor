library("BioCor")
context("Testing geneSim")


test_that("geneSim", {

    test <- geneSim("2", "9", info)
    expect_equal(test, 0.4)
    test <- geneSim("2", "9", info, NULL)
    expect_false(is.null(colnames(test)))
    expect_false(is.null(rownames(test)))
    expect_equal(ncol(test), 4L)
})

test_that("mgeneSim", {
    test <- mgeneSim(c("1", "10", "2", "3", "4", "5", "6", "7", "8", "9"), info)
    expect_true(is.na(test["1", "1"]))
    expect_true(isSymmetric(test))
    expect_equal(test["2", "10"], 0.4)
})
