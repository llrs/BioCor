library("BioCor")
context("Testing geneSim")


test_that("geneSim", {

    test <- geneSim("2", "9", info)
    test2 <- geneSim("9", "2", info)
    expect_equal(test, 0.4)
    expect_equal(test, test2)
    test <- geneSim("2", "9", info, NULL)
    expect_true(is.matrix(test))
    expect_false(is.null(colnames(test)))
    expect_false(is.null(rownames(test)))
    expect_equal(ncol(test), 4L)
    expect_error(geneSim(c("2", "9"), "9", info), "just one gene")
    expect_true(is.na(geneSim("11", "12", info)))
    expect_error(geneSim("2", "9", as.environment(info)))
    expect_true(is.na(geneSim("7", "8", info)))
})

test_that("mgeneSim", {
    test <- mgeneSim(c("1", "10", "2", "3", "4", "5", "6", "7", "8", "9"),
                     info)
    expect_true(is.na(test["1", "1"]))
    expect_true(isSymmetric(test))
    expect_equal(test["2", "10"], 0.4)
    expect_error(mgeneSim(c("2", "2"), info))
    expect_error(geneSim(c("2", "9"), as.environment(info)))
    test2 <- mgeneSim(c("7", "8"), info)
    expect_true(is.na(test2[1L, 1L]))
    expect_error(mgeneSim(c("11", "12"), info), "provided")
    test3 <- mgeneSim(c("a" = "2", "b" = "9"), info)
    expect_equal(colnames(test3), letters[1:2])
    expect_equal(rownames(test3), letters[1:2])
    expect_error(mgeneSim(c(2, 3), info), "character")
})
