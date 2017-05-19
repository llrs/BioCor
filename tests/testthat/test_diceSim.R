library("BioCor")
context("Testing diceSim")

test_that("diceSim", {
    genes.id2 <- c("52", "11342", "80895", "57654", "548953", "11586", "45985")
    genes.id <- c("52", "11342", "80895", "57654", "58493", "1164", "1163",
                  "4150", "2130", "159")
    test <- diceSim(genes.id, genes.id2)
    expect_equal(test, 0.47058823)
    expect_equal(diceSim(genes.id, genes.id), 1L)
    expect_true(is.na(diceSim(c(), c())))
    expect_true(is.na(diceSim(c(), genes.id)))
})
