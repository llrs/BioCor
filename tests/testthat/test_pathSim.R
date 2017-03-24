library("BioCor")
context("Testing pathSim")

test_that("pathSim", {
    test <- pathSim("194840", "156580", info)
    expect_equal(test, 0.4)
    expect_error(pathSim(NULL, "156580", info))
    expect_error(pathSim(1, "156580", info))
    expect_error(pathSim("1", "156580", "a"))
})

test_that("mpathSim", {
    test <- mpathSim(c("194840", "156580"), info, method = "max")
    expect_equal(test, 1L)
    test <- mpathSim(c("194840", "156580"), info, NULL)
    expect_equal(test[2L, 1L], 0.4)
    expect_error(mpathSim("156580", info), "Introduce several")
    expect_error(mpathSim(1, "156580", info))
    expect_error(mpathSim("a", "156580", "a"))
    expect_error(mpathSim("a", "156580", info))
})

test_that("pathSims_matrix", {
    test <- pathSims_matrix(info)
    pathways <- unique(unlist(info))
    nPathways <- length(pathways[!is.na(pathways)])
    expect_equal(dim(test), c(nPathways, nPathways))
    expect_true(all(diag(test) == 1L))
    eq <- pathSim("1430728", "156580", info)
    expect_equal(test["1430728", "156580"], eq)
})
