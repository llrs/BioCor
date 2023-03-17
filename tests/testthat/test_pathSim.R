library("BioCor")
context("Testing pathSim")

test_that("pathSim", {
    test <- pathSim("194840", "156580", info)
    expect_equal(test, 0.4)
    expect_error(pathSim(NULL, "156580", info), "several")
    expect_error(pathSim(1, "156580", info), "characters")
    expect_error(pathSim("1", "156580", "a"), "list")
    expect_true(is.na(pathSim("1", "156580", info)))
    expect_equal(pathSim("109582", "114608", info), 1)
})

test_that("pathSim with GSC", {
    test <- pathSim("194840", "156580", Info)
    expect_equal(test, 0.4)
    expect_error(pathSim(NULL, "156580", Info), "several")
    expect_error(pathSim(1, "156580", Info), "characters")
    expect_error(pathSim("1", "156580", "a"), "list")
    expect_true(is.na(pathSim("1", "156580", Info)))
    expect_true(is.na(pathSim("109582", "114608", Info)))
})
