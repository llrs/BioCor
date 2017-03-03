library("BioCor")
context("Testing clusterSim")

test_that("clusterSim", {

    test <- clusterSim(c("2", "1"), c("9", "4"), info)
    test2 <- clusterSim(c("9", "4"), c("2", "1"), info)
    expect_equal(test, 0.4)
    expect_equal(test2, test)
    test <- clusterSim(c("2", "1"), c("9", "4"), info, "avg")
    expect_equal(test, 0.00909090909090909)
    test <- clusterSim(c("2", "1"), c("9", "4"), info, "avg", round = TRUE)
    expect_equal(test, 0.009)
})

test_that("clusterSim", {
    clusters <- list(cluster1 = c("10", "2", "3"),
                     cluster2 = c("10", "2", "9"))
    test <- mclusterSim(clusters, info)
    expect_true(isSymmetric(test))
    expect_equal(test[1, 1], 1)
    expect_equal(colnames(test), names(clusters))
    expect_equal(rownames(test), names(clusters))
    clusters <- list(cluster1 = c("10", "3"),
                     cluster2 = c("10", "2", "9"),
                     cluster3 = c("2", "9", "3", "4"))
    test <- mclusterSim(clusters, info)
    expect_equal(test[1, 1], 1)
    expect_equal(colnames(test), names(clusters))
    expect_equal(rownames(test), names(clusters))
    test <- mclusterSim(clusters, info, "avg", round = TRUE)
    expect_equal(test[1, 1], 0.664)

    expect_warning(mclusterSim(clusters, info, method = NULL, round = TRUE))
})
