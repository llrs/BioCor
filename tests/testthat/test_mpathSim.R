test_that("mpathSim", {
    test <- mpathSim(c("194840", "156580"), info, method = "max")
    expect_equal(test, 1L)
    test <- mpathSim(c("194840", "156580"), info, NULL)
    expect_equal(test[2L, 1L], 0.4)
    test <- mpathSim(c("a" = "194840", "b" = "156580"), info, method = NULL)
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

test_that("mpathSim for GeneSetCollections and list is equal", {

    # With missing all
    a <- mpathSim(info = inverseList(geneIds(gss)))
    b <- mpathSim(info = gss)
    expect_equal(a, b)
    pathways <- colnames(a)

    # Without missing pathways
    a <- mpathSim(names(geneIds(gss)), info = inverseList(geneIds(gss)))
    b <- mpathSim(names(geneIds(gss)), info = gss)
    expect_equal(a, b)


    # Without missing pathways
    expect_warning(a <- mpathSim(c(pathways, "A"),
                                 info = inverseList(geneIds(gss))
    ),
    "not in the list"
    )
    expect_warning(b <- mpathSim(c(pathways, "A"), info = gss),
                   "not in the GeneSetCollection"
    )
    expect_equal(a, b)

    # Without missing method
    a <- mpathSim(info = inverseList(geneIds(gss)), method = "avg")
    b <- mpathSim(info = gss, method = "avg")
    expect_equal(a, b)

    # Without missing pathways and method
    a <- mpathSim(names(geneIds(gss)),
        info = inverseList(geneIds(gss)),
        method = "avg"
    )
    b <- mpathSim(names(geneIds(gss)), info = gss, method = "avg")
    expect_equal(a, b)

    # Without missing pathways and methods
    expect_warning(a <- mpathSim(c(pathways, "A"),
            info = inverseList(geneIds(gss))
        ),
        "not in the list"
    )
    expect_warning(b <-
        mpathSim(c(pathways, "A"), info = gss),
        "not in the GeneSetCollection"
    )
    expect_equal(a, b)
})

test_that("mpathSim for GeneSetCollections", {
    a <- mpathSim(info = gss)
    expect_true(all(diag(a) == 1))
    expect_equal(colnames(a), rownames(a))
    expect_warning(
      a <- mpathSim(c(colnames(a), "A"), info = gss),
      "not in the GeneSetCollection"
    )
    expect_true(!all(is.na(diag(a) == 1)))
    expect_equal(colnames(a), rownames(a))
})
