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
    expect_error(clusterSim(c("2", "1"), c("9", "4"), "a"), "list")
    expect_warning(clusterSim(c("2", "11"), c("9", "4"), info), "not in")
    expect_true(is.na(clusterSim(c("4", "5"), c("6", "7"), info)))
    test <- clusterSim(c("2", "1"), c("9", "4"), info, NULL)
    expect_true(is.matrix(test))
    expect_true(is.na(clusterSim(c("1", "5"), c("7", "9"), info)))
})

test_that("mclusterSim", {
    clusters <- list(cluster1 = c("10", "2", "3"),
                     cluster2 = c("10", "2", "9"))
    test <- mclusterSim(clusters, info)
    expect_equal(test[1L, 1L], 1L)
    expect_equal(colnames(test), names(clusters))
    expect_equal(rownames(test), names(clusters))


    expect_error(mclusterSim(c("10", "2"), info), "list")
    test <- mclusterSim(list(a = c("4", "5"), b = c("6", "7")), info)
    expect_true(all(is.na(unlist(test))))
    expect_true(isSymmetric(test))
    clusters <- list(cluster1 = c("10", "3"),
                     cluster2 = c("10", "2", "9"),
                     cluster3 = c("2", "9", "3", "4"))
    test <- mclusterSim(clusters, info)
    expect_equal(test[1L, 1L], 1L)
    expect_equal(colnames(test), names(clusters))
    expect_equal(rownames(test), names(clusters))
    test <- mclusterSim(clusters, info, "avg", round = TRUE)
    expect_equal(test[1L, 1L], 0.7) # 0.664
    clusters3 <- list(cluster1 = c("10", "2", "3"),
                     cluster2 = c(10, 2, 9))
    expect_error(mclusterSim(clusters3, info))
    clusters2 <- list(cluster1 = c("10", "2", "3"))
    expect_error(mclusterSim(clusters2, info), "several")
    expect_error(mclusterSim(clusters, as.environment(info)), "list")
    expect_warning(mclusterSim(list(cluster1 = c("10", "2", "3"),
                                  cluster2 = c("13", "2", "9")), info),
                   "in the list")

    clusters <- list(cluster1 = c("4", "2"),
                     cluster2 = c("10", "2", "9"),
                     cluster3 = c("4", "9", "10"))
    test4 <- mclusterSim(clusters, info, "rcmax.avg")
    expect_equal(test4["cluster1", "cluster2"], 0.76)
    expect_warning(mclusterSim(clusters, info, method = NULL, round = TRUE))
})

