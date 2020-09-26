library("BioCor")
context("Testing combineScores")
methods <- c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")

test_that("combineScores", {
    test <- sapply(methods, combineScores, scores = d)
    expect_equal(as.numeric(test["rcmax.avg"]), as.numeric(test["BMA"]))
    expect_true(all(test["max"] > test[-2]))
    expect_equal(as.numeric(test["avg"]), 0.339770723104056)
    expect_equal(as.numeric(test["max"]), 0.6)
    expect_equal(as.numeric(test["rcmax"]), 0.5)
    expect_equal(as.numeric(test["rcmax.avg"]), 0.464285714285714)
    expect_equal(as.numeric(test["reciprocal"]), 1 / 3)

    expect_error(combineScores(d, "reciprocal", t = "a"), "between")
    expect_error(combineScores(d, "reciprocal", t = 2), "between")
    expect_error(combineScores(d, "reciprocal", t = -2), "between")
    expect_true(is.na(combineScores(d, "reciprocal", t = 0.61)))
    d2 <- list("d" = d)
    expect_equal(combineScores(d2, "max"), 0.6)

    # Rounding
    test2 <- sapply(methods, combineScores, scores = d, round = TRUE)
    expect_equal(as.numeric(test2["rcmax.avg"]), as.numeric(test2["BMA"]))
    expect_true(all(test2["max"] > test2[-2]))
    expect_equal(as.numeric(test2["avg"]), 0.340)
    expect_equal(as.numeric(test2["max"]), 0.6)
    expect_equal(as.numeric(test2["rcmax"]), 0.5)
    expect_equal(as.numeric(test2["rcmax.avg"]), 0.464)
    expect_equal(as.numeric(test2["reciprocal"]), 0.333)
    expect_error(combineScores(as.vector(d), "max"), "matrix")
    expect_error(combineScores(d, "max", round = "a"), "logical")

    # Repeated "max values"
    d[1, 1] <- 0.5
    test3 <- sapply(methods, combineScores, scores = d)
    expect_equal(as.numeric(test3["rcmax.avg"]), as.numeric(test3["BMA"]))
    expect_true(all(test3["max"] > test3[-2]))
    expect_equal(as.numeric(test3["avg"]), 0.350881834215168)
    expect_equal(as.numeric(test3["max"]), 0.6)
    expect_equal(as.numeric(test3["rcmax"]), 0.5)
    expect_equal(as.numeric(test3["rcmax.avg"]), 0.480952380952381)
    expect_equal(as.numeric(test3["reciprocal"]), 0.2)
    # Nas in the matrix
    d[1, 1] <- NA
    test4 <- sapply(methods, combineScores, scores = d)
    expect_equal(as.numeric(test4["rcmax.avg"]), as.numeric(test4["BMA"]))
    expect_true(all(test4["max"] > test4[-2]))
    expect_equal(as.numeric(test4["avg"]), 0.332242063492063)
    expect_equal(as.numeric(test4["max"]), 0.6)
    expect_equal(as.numeric(test4["rcmax"]), 0.5)
    expect_equal(as.numeric(test4["rcmax.avg"]), 0.464285714285714)
    expect_equal(as.numeric(test4["reciprocal"]), 1 / 3)
    # Not square matrices
    scores2 <- structure(c(
        0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
        0.285714285714286, 0.13, 0.2, 0.6
    ),
    .Dim = c(4L, 3L),
    .Dimnames = list(
        c("a", "b", "c", "d"),
        c("e", "f", "g")
    )
    )
    test5 <- sapply(methods, combineScores, scores = scores2)
    expect_equal(as.numeric(test5["rcmax.avg"]), as.numeric(test5["BMA"]))
    expect_true(all(test5["max"] > test5[-2]))
    expect_equal(as.numeric(test5["avg"]), 0.332328042328042)
    expect_equal(as.numeric(test5["max"]), 0.6)
    expect_equal(as.numeric(test5["rcmax"]), 0.566666666666667)
    expect_equal(as.numeric(test5["rcmax.avg"]), 0.507142857142857)
    expect_equal(as.numeric(test5["reciprocal"]), 0.342857142857143)
})

test_that("combineScoresPar", {
    methods <- c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")
    e <- structure(c(
        0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
        0.285714285714286
    ),
    .Dim = c(3L, 3L),
    .Dimnames = list(c("a", "b", "c"), c("a", "b", "c"))
    )
    subSet <- list(a = c("a", "b"), b = c("b", "c"))
    test3 <- combineScoresPar(e, method = "max", subSet)
    expect_that(test3, is_a("matrix"))
})

test_that("combineScoresPar equivalent to combineScores", {
    methods <- c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")
    test <- sapply(methods, combineScoresPar, scores = d)
    test_s <- sapply(methods, combineScores, scores = d)
    expect_equal(test, test_s)

    expect_error(combineScoresPar(d, "reciprocal", t = "a"), "between")
    expect_error(combineScoresPar(d, "reciprocal", t = 2), "between")
    expect_error(combineScoresPar(d, "reciprocal", t = -2), "between")
    expect_true(is.na(combineScoresPar(d, "reciprocal", t = 0.61)))
    d2 <- list("d" = d)
    expect_equal(combineScoresPar(d2, "max"), combineScores(d2, "max"))

    # Rounding
    test2 <- sapply(methods, combineScoresPar, scores = d, round = TRUE)
    test2_s <- sapply(methods, combineScores, scores = d, round = TRUE)
    expect_equal(test2, test2_s)

    expect_error(combineScoresPar(as.vector(d), "max"), "matrix")
    expect_error(combineScoresPar(d, "max", round = "a"), "logical")

    # Repeated "max values"
    d[1, 1] <- 0.5
    test3 <- sapply(methods, combineScoresPar, scores = d)
    test3_s <- sapply(methods, combineScores, scores = d)
    expect_equal(test3, test3_s)

    # Nas in the matrix
    d[1, 1] <- NA
    test4 <- sapply(methods, combineScoresPar, scores = d)
    test4_s <- sapply(methods, combineScores, scores = d)
    expect_equal(test4, test4_s)

    # Not square matrices
    scores2 <- structure(c(
        0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
        0.285714285714286, 0.13, 0.2, 0.6
    ),
    .Dim = c(4L, 3L),
    .Dimnames = list(
        c("a", "b", "c", "d"),
        c("e", "f", "g")
    )
    )
    test5 <- sapply(methods, combineScoresPar, scores = scores2)
    test5_s <- sapply(methods, combineScores, scores = scores2)
    expect_equal(test5, test5_s)
})

test_that("reciprocal", {
    expect_error(reciprocal(d, 2), "between")
    expect_error(reciprocal(as.vector(d), 0.5), "matrix")

    p2 <- structure(c(
        1, 0.8, 1, 1, 0.8, 1, 0.8, 0.8, 1, 0.8, 1, 1, 1,
        0.8, 1, 1
    ),
    .Dim = c(4L, 4L),
    .Dimnames = list(
        c("1430728", "156580", "156582", "211859"),
        c("1430728", "156580", "156582", "211859")
    )
    )
    scores <- structure(c(
        0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
        0.285714285714286
    ),
    .Dim = c(3L, 3L),
    .Dimnames = list(c("a", "b", "c"), c("d", "e", "f"))
    )
    scores2 <- scores
    scores2[1, 1] <- 0.5
    expect_equal(reciprocal(p2, 0L), 1L) # Test many draws
    expect_equal(reciprocal(scores, 0), 1 / 3)
    expect_equal(reciprocal(scores2, 0), 0.2)
})

test_that("BMA", {
    expect_error(BMA(as.vector(d)), "matrix")
})

test_that("rcmax", {
    expect_error(rcmax(as.vector(d)), "matrix")
})

test_that("rowIds  &colIds are not null in combineScorePar", {
    set.seed(456)
    library("reactome.db")
    genes2Pathways <- as.list(reactomeEXTID2PATHID)
    pathways <- unlist(genes2Pathways, use.names = FALSE)
    genes <- rep(names(genes2Pathways), lengths(genes2Pathways))
    paths2genes <- split(genes, pathways)
    human <- grep("R-HSA-", names(paths2genes))
    paths2genes <- paths2genes[human]
    paths2genes <- lapply(paths2genes, unique)
    paths2genes <- paths2genes[lengths(paths2genes) >= 2]

    a <- unlist(paths2genes, use.names = FALSE)
    b <- rep(names(paths2genes), lengths(paths2genes))
    genes2paths <- split(b, a)

    # clusters
    clusters <- list(a = sample(genes, 50), b = sample(genes, 25))
    expect_warning(
        mclusterGeneSim(clusters,
            info = genes2paths,
            method = c("max", "BMA")
        ),
        "Some genes are not in the list provided."
    )
})
