library("BioCor")
context("Testing combineScores")

test_that("combineScores", {
    test <- sapply(c("avg", "max", "rcmax", "rcmax.avg", "BMA"), combineScores,
                   scores = d)
    expect_equal(as.numeric(test["rcmax.avg"]), as.numeric(test["BMA"]))
    expect_true(all(test["max"] > test[-2]))
    expect_equal(as.numeric(test["avg"]), 0.339770723104056)
    expect_equal(as.numeric(test["max"]), 0.6)
    expect_equal(as.numeric(test["rcmax"]), 0.5)
    expect_equal(as.numeric(test["rcmax.avg"]), 0.464285714285714)

    test2 <- sapply(c("avg", "max", "rcmax", "rcmax.avg", "BMA"), combineScores,
                    scores = d, round = TRUE)
    expect_equal(as.numeric(test2["rcmax.avg"]), as.numeric(test2["BMA"]))
    expect_true(all(test2["max"] > test2[-2]))
    expect_equal(as.numeric(test2["avg"]), 0.340)
    expect_equal(as.numeric(test2["max"]), 0.6)
    expect_equal(as.numeric(test2["rcmax"]), 0.5)
    expect_equal(as.numeric(test2["rcmax.avg"]), 0.464)
    expect_error(combineScores(as.vector(d), "max"), "matrix")
    expect_error(combineScores(d, "max", round = "a"), "logical")
})
