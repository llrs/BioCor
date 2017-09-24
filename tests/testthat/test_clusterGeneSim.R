library("BioCor")
context("Testing clusterGeneSim")

test_that("clusterGeneSim", {

    expect_warning(test <- clusterGeneSim(c("2", "1"), c("9", "4"), info),
                   "Using max method because after removing NAs")
    expect_error(clusterGeneSim(c("2", "2"), c("4", "4"), info), "several")
    expect_error(clusterGeneSim(c("2", "2"), c(9, 4), info), "character")
    expect_error(clusterGeneSim(c("2", "2"), c("9", "4"),
                                as.environment(info)),
                 "list")
    expect_warning(clusterGeneSim(c("13", "12"), c("9", "4"), info), "list")
    expect_error(clusterGeneSim(c("2", "1"), c("9", "3"), info,
                                method = NULL), "method")
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
    expect_error(mclusterGeneSim(c("a", "b"), info), "list")
    expect_error(mclusterGeneSim(list(a = c("a", "b"), b = c(1, 2)), info),
                 "character")
    expect_error(mclusterGeneSim(list(a = c("a", "b")), info),
                 "several clusters")
    expect_error(mclusterGeneSim(clusters, c("a", "b")), "list")
    cluster2 <- clusters
    cluster2$cluster1 <- "199"
    expect_warning(mclusterGeneSim(cluster2, info), "not in")
    test <- mclusterGeneSim(clusters, info)
    expect_equal(test[1L, 1L], 1)

    test <- mclusterGeneSim(clusters, info, c("max", "avg"))
    expect_equal(test[1L, 1L], 1)
    expect_true(isSymmetric(test))
    expect_error(mclusterGeneSim(clusters, info, NULL))
})
