library("BioCor")
context("Testing clusterGeneSim")

test_that("clusterGeneSim", {

    expect_warning(test0 <- clusterGeneSim(c("2", "1"), c("9", "4"), info),
                   "Using max method because after removing NAs")
    expect_error(clusterGeneSim(c("2", "2"), c("4", "4"), info), "several")
    expect_error(clusterGeneSim(c("2", "2"), c(9, 4), info), "character")
    expect_error(clusterGeneSim(c("2", "2"), c("9", "4"),
                                as.environment(info)),
                 "list")
    expect_warning(test <- clusterGeneSim(c("13", "12"), c("9", "4"), info), "list")
    expect_true(is.na(test))
    expect_error(clusterGeneSim(c("2", "1"), c("9", "3"), info,
                                method = NULL), "method")
    expect_warning(test2 <- clusterGeneSim(c("9", "4"), c("2", "1"), info))
    expect_equal(test0, 0.4)
    expect_equal(test0, test2)
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



test_that("clusterGeneSim GeneSetCollection", {

    expect_warning(test0 <- clusterGeneSim(c("2", "1"), c("9", "4"), Info),
                   "Some genes are not in the GeneSetCollection")
    expect_error(clusterGeneSim(c("2", "2"), c("4", "4"), Info), "several")
    expect_error(clusterGeneSim(c("2", "2"), c(9, 4), Info), "character")
    expect_warning(test <- clusterGeneSim(c("13", "12"), c("9", "4"), Info), "list")
    expect_true(is.na(test))
    expect_warning(expect_error(clusterGeneSim(c("2", "1"), c("9", "3"), Info,
                                method = NULL), "method"))
    expect_warning(test2 <- clusterGeneSim(c("9", "4"), c("2", "1"), Info))
    expect_equal(test0, 0.4)
    expect_equal(test0, test2)
    expect_warning(test <- clusterGeneSim(c("2", "1"), c("9", "4"), Info, "max"))
    expect_true(is.matrix(test))
    expect_warning(test <- clusterGeneSim(c("2", "1"), c("9", "4"), Info, c("avg", "max")))
    expect_equal(test, 0.1)
    expect_warning(test <- clusterGeneSim(c("2", "1"), c("9", "4"), Info, c("avg", "max"),
                           round = TRUE))
    expect_equal(test, 0.1)

})

test_that("mclusterGeneSim GeneSetCollection", {
    expect_error(mclusterGeneSim(c("a", "b"), Info), "list")
    expect_error(mclusterGeneSim(list(a = c("a", "b"), b = c(1, 2)), Info),
                 "character")
    expect_error(mclusterGeneSim(list(a = c("a", "b")), Info),
                 "several clusters")
    expect_error(mclusterGeneSim(clusters, c("a", "b")), "list")
    cluster2 <- clusters
    cluster2$cluster1 <- "199"
    expect_warning(mclusterGeneSim(cluster2, Info), "not in")
    expect_warning(test <- mclusterGeneSim(clusters, Info))
    expect_equal(test[1L, 1L], 1)

    expect_warning(test <- mclusterGeneSim(clusters, Info, c("max", "avg")))
    expect_equal(test[1L, 1L], 1)
    expect_true(isSymmetric(test))
    expect_error(expect_warning(mclusterGeneSim(clusters, Info, NULL)))
})
