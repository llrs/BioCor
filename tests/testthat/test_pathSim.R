library("BioCor")
context("Testing pathSim")

test_that("pathSim", {
    test <- pathSim("194840", "156580", info)
    expect_equal(test, 0.4)
    expect_error(pathSim(NULL, "156580", info), "several")
    expect_error(pathSim(1, "156580", info), "characters")
    expect_error(pathSim("1", "156580", "a"), "list")
    expect_true(is.na(pathSim("1", "156580", info)))

})

test_that("mpathSim", {
    test <- mpathSim(c("194840", "156580"), info, method = "max")
    expect_equal(test, 1L)
    test <- mpathSim(c("194840", "156580"), info, NULL)
    expect_equal(test[2L, 1L], 0.4)
    test <- mpathSim(c("a" = "194840","b" = "156580"), info, method = NULL)
    expect_equal(test["a", "b"], 0.4)
    expect_error(mpathSim("156580", info), "Introduce several")
    expect_error(mpathSim(c(1, 156580), info), "character")
    expect_error(mpathSim(c("a", "156580"), "a"), "list")
    expect_warning(mpathSim(c("a", "156580"), info), "not in")
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


test_that("mpathSim for GeneSetCollections", {
    fl <- system.file("extdata", "Broad.xml", package="GSEABase")
    gss <- getBroadSets(fl) # GeneSetCollection of 2 sets

    a <- mpathSim(names(gss), info = gss)
    expect_true(all(diag(a) == 1))
    expect_equal(colnames(a), rownames(a))

    a <- expect_warning(mpathSim(c(names(gss), "A"), info = gss),
                   "not in the GeneSetCollection")
    expect_true(is.na(all(diag(a) == 1)))
    expect_equal(colnames(a), rownames(a))
})
