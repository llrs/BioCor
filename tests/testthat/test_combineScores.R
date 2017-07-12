library("BioCor")
context("Testing combineScores")

test_that("combineScores", {
    methods <- c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")
    test <- sapply(methods, combineScores, scores = d)
    expect_equal(as.numeric(test["rcmax.avg"]), as.numeric(test["BMA"]))
    expect_true(all(test["max"] > test[-2]))
    expect_equal(as.numeric(test["avg"]), 0.339770723104056)
    expect_equal(as.numeric(test["max"]), 0.6)
    expect_equal(as.numeric(test["rcmax"]), 0.5)
    expect_equal(as.numeric(test["rcmax.avg"]), 0.464285714285714)
    expect_equal(as.numeric(test["reciprocal"]), 1/3)

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
    d[1,1] <- 0.5
    test3 <- sapply(methods, combineScores, scores = d)
    expect_equal(as.numeric(test3["rcmax.avg"]), as.numeric(test3["BMA"]))
    expect_true(all(test3["max"] > test3[-2]))
    expect_equal(as.numeric(test3["avg"]), 0.350881834215168)
    expect_equal(as.numeric(test3["max"]), 0.6)
    expect_equal(as.numeric(test3["rcmax"]), 0.5)
    expect_equal(as.numeric(test3["rcmax.avg"]), 0.480952380952381)
    expect_equal(as.numeric(test3["reciprocal"]), 0.2)
    # Nas in the matrix
    d[1,1] <- NA
    test4 <- sapply(methods, combineScores, scores = d)
    expect_equal(as.numeric(test4["rcmax.avg"]), as.numeric(test4["BMA"]))
    expect_true(all(test4["max"] > test4[-2]))
    expect_equal(as.numeric(test4["avg"]), 0.332242063492063)
    expect_equal(as.numeric(test4["max"]), 0.6)
    expect_equal(as.numeric(test4["rcmax"]), 0.5)
    expect_equal(as.numeric(test4["rcmax.avg"]), 0.464285714285714)
    expect_equal(as.numeric(test4["reciprocal"]), 1/3)
    # Not square matrices
    scores2 <- structure(c(0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
                           0.285714285714286, 0.13, 0.2, 0.6),
                         .Dim = c(4L, 3L),
                         .Dimnames = list(c("a", "b", "c", "d"),
                                          c("e", "f", "g")))
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
    e <- structure(c(0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
                      0.285714285714286), .Dim = c(3L, 3L),
                    .Dimnames = list(c("a", "b", "c"), c("a", "b", "c")))
    subSet <- list(a = c("a", "b"), b = c("b", "c"))
    expect_warning(combineScoresPar(e, subSet, NULL,
                     method = "max"), "symmetric")
    expect_warning(test3 <- combineScoresPar(e, subSet, bpparam(1),
                                             method = "max"),
                   "symmetric")
    class(test3)
})
