library("BioCor")
context("Testing geneSim")


test_that("geneSim", {
    test <- geneSim("2", "9", info)
    test2 <- geneSim("9", "2", info)
    expect_error(geneSim(2, "9", info), "character")
    expect_equal(test, 0.4)
    expect_equal(test, test2)
    test <- geneSim("2", "9", info, NULL)
    expect_true(is.matrix(test))
    expect_false(is.null(colnames(test)))
    expect_false(is.null(rownames(test)))
    expect_equal(ncol(test), 4L)
    expect_error(geneSim(c("2", "9"), "9", info), "just one gene")
    expect_true(is.na(geneSim("11", "12", info)))
    expect_error(geneSim("2", "9", as.environment(info)), "list")
    expect_true(is.na(geneSim("7", "8", info)))
})

test_that("mgeneSim", {
    test <- mgeneSim(
        c("1", "10", "2", "3", "4", "5", "6", "7", "8", "9"),
        info
    )
    expect_true(is.na(test["1", "1"]))
    expect_true(isSymmetric(test))
    expect_equal(test["2", "10"], 0.4)
    expect_error(mgeneSim(c("2", "2"), info), "geneSim")
    expect_error(mgeneSim(c("2", "9"), as.environment(info)), "list")
    test2 <- mgeneSim(c("7", "8"), info)
    expect_true(is.na(test2[1L, 1L]))
    expect_error(mgeneSim(c("11", "12"), info), "provided")
    expect_warning(mgeneSim(c("11", "10"), info), "not in")
    test3 <- mgeneSim(c("a" = "2", "b" = "9"), info)
    expect_equal(colnames(test3), letters[1:2])
    expect_equal(rownames(test3), letters[1:2])
    expect_error(mgeneSim(c(2, 3), info), "character")
    expect_warning(mgeneSim(as.character(1:10), info, NULL), "be null")
    expect_error
})

test_that("geneSim with GSC", {
    test <- geneSim("2", "9", Info)
    test2 <- geneSim("9", "2", Info)
    expect_error(geneSim(2, "9", Info), "character")
    expect_equal(test, 0.4)
    expect_equal(test, test2)
    test3 <- geneSim("9", "2", info)
    expect_equal(test2, test3)

    test <- geneSim("2", "9", Info, NULL)
    expect_true(is.matrix(test))
    expect_false(is.null(colnames(test)))
    expect_false(is.null(rownames(test)))
    expect_equal(ncol(test), 4L)
    expect_error(geneSim(c("2", "9"), "9", Info), "just one gene")
    expect_true(is.na(geneSim("11", "12", Info)))
    expect_error(geneSim("2", "9", as.environment(Info)), " character vector")
    expect_true(is.na(geneSim("7", "8", Info)))
})

test_that("mgeneSim with GSC", {
    test <- mgeneSim(c("10", "2", "3", "9", "1"), Info)
    expect_true(is.na(test))

    test <- mgeneSim(c("10", "2", "3", "9"), Info)
    expect_true(isSymmetric(test))
    expect_equal(test["2", "10"], 0.4)
    expect_error(mgeneSim(c("2", "2"), Info), "geneSim")
    test2 <- mgeneSim(c("7", "8"), Info)
    expect_true(is.na(test2))
    test3 <- mgeneSim(c("a" = "2", "b" = "9"), Info)
    expect_equal(colnames(test3), letters[1:2])
    expect_equal(rownames(test3), letters[1:2])
    expect_error(mgeneSim(c(2, 3), Info), "character")
})
