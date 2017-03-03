library("BioCor")
context("Testing clusterGeneSim")

test_that("clusterGeneSim", {

    expect_warning(test <- clusterGeneSim(c("2", "1"), c("9", "4"), info))
    expect_warning(test2 <- clusterGeneSim(c("9", "4"), c("2", "1"), info))
    expect_equal(test, 0.4)
    expect_equal(test, test2)
    test <- clusterGeneSim(c("2", "1"), c("9", "4"), info, "max")
    expect_true(is.matrix(test))
    test <- clusterGeneSim(c("2", "1"), c("9", "4"), info, c("avg", "max"))
    expect_equal(test, 0.00909090909090909)
    test <- clusterGeneSim(c("2", "1"), c("9", "4"), info, c("avg", "max"),
                           round = TRUE)
    expect_equal(test, 0.009)

})

test_that("mclusterGeneSim", {
    test <- mclusterGeneSim(clusters, info)
    expect_equal(test[1L, 1L], 1)

    test <- mclusterGeneSim(clusters, info, c("max", "avg"))
    expect_equal(test[1L, 1L], 1)
    expect_true(isSymmetric(test))
    expect_error(mclusterGeneSim(clusters, info, NULL))
})