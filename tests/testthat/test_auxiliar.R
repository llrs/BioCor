library("BioCor")
context("Testing auxiliar functions")

test_that("seq2mat puts a combination to the right place", {
  test <- seq2mat(LETTERS[1:5], 1:10)

  expect_true(isSymmetric(test))

  expect_error(seq2mat(LETTERS[1:5], 1:11))
  expect_error(seq2mat(LETTERS[1:5], 1:9))

  com <- combn(as.character(1:5), 2)
  dat <- apply(com, 2, function(x){sum(as.numeric(x))})

  expect_equal(seq2mat(as.character(1:5), dat)["3", "5"], 8L)
  expect_equal(seq2mat(as.character(1:5), dat)["5", "3"], 8L)
})

test_that("comb_biopath", {})
test_that("react_genes", {})
test_that("indices.dup", {})
test_that(".combinadic", {
  bcombin <- combn(LETTERS[1:5], 2)
  .combin <- .combinadic(LETTERS[1:5], 2, 1)
  expect_true(all.equal(bcombin[, 1], rev(.combin)))
  .combin <- .combinadic(LETTERS[1:5], 2, 2)
  expect_true(all.equal(bcombin[, 2], rev(.combin)))

})
test_that("rem.dup", {})
test_that("indices.dup", {
  test <- indices.dup(as.character(c(vec, 3)))
  expect_true(is.list(test))
  expect_length(test$`3`, 3L)
})

test_that("removeDup", {
  a <- seq2mat(c("52", "52", "53", "55"), rnorm(choose(4, 2)))
  mat <- list("kegg" = a, "react" = a)
  dupli <- indices.dup(rownames(a))
  remat <- removeDup(mat, dupli)

  expect_equal(ncol(remat[[1]]), ncol(mat[[1]]) - 1)
  expect_equal(nrow(remat[[1]]), nrow(mat[[1]]) - 1)
  expect_equal(nrow(remat[[1]]), ncol(remat[[1]]))
})
