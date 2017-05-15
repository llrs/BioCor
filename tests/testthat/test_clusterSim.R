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
    expect_error(clusterSim(c("9", "9"), c("2", "2"), info), "several")
    expect_error(clusterSim(c("9", "9"), c(1, 1, 2), info), "character")
    expect_error(clusterSim(c("9", "9"), c(1, 1), info), "several")
    expect_warning(clusterSim(c("11", "13"), c("12", "15"), info), "provided")
})

test_that("clusterSim", {
    clusters <- list(cluster1 = c("10", "2", "3"),
                     cluster2 = c("10", "2", "9"))
    test <- mclusterSim(clusters, info)
    expect_true(isSymmetric(test))
    expect_equal(test[1L, 1L], 1L)
    expect_equal(colnames(test), names(clusters))
    expect_equal(rownames(test), names(clusters))
    clusters <- list(cluster1 = c("10", "3"),
                     cluster2 = c("10", "2", "9"),
                     cluster3 = c("2", "9", "3", "4"))
    test <- mclusterSim(clusters, info)
    expect_equal(test[1L, 1L], 1L)
    expect_equal(colnames(test), names(clusters))
    expect_equal(rownames(test), names(clusters))
    test <- mclusterSim(clusters, info, "avg", round = TRUE)
    expect_equal(test[1L, 1L], 0.664)
    clusters3 <- list(cluster1 = c("10", "2", "3"),
                     cluster2 = c(10, 2, 9))
    expect_error(mclusterSim(clusters3, info))
    clusters2 <- list(cluster1 = c("10", "2", "3"))
    expect_error(mclusterSim(clusters2, info), "several")
    expect_error(mclusterSim(clusters, as.environment(info)), "list")
    expect_warning(mclusterSim(list(cluster1 = c("10", "2", "3"),
                                  cluster2 = c("13", "2", "9")), info),
                   "in the list")
    expect_warning(mclusterSim(clusters, info, method = NULL, round = TRUE))
})

