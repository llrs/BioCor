library("BioCor")
context("Testing auxiliar functions")

test_that("seq2mat puts a combination to the right place", {
  test <- seq2mat(LETTERS[1:5], 1:10)

  expect_true(isSymmetric(test))

  expect_error(seq2mat(LETTERS[1:5], 1:11))
  expect_error(seq2mat(LETTERS[1:5], 1:9))

  com <- combn(as.character(1:5), 2)
  dat <- apply(com, 2, function(x){sum(as.numeric(x))})

  expect_equal()
})
test_that("comb_biopath", {})
test_that("react_genes", {})
test_that("indices.dup", {})
test_that(".combinadic", {})
