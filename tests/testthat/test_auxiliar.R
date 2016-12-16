library("BioCor")
context("Testing auxiliar functions")

# seq2mat ####
test_that("seq2mat puts a combination to the right place", {
  test <- seq2mat(LETTERS[1:5], 1:10)

  expect_true(isSymmetric(test)) # already checks if it is squared.

  expect_error(seq2mat(LETTERS[1:5], 1:11))
  expect_error(seq2mat(LETTERS[1:5], 1:9))


  com <- combn(as.character(1:5), 2)
  dat <- apply(com, 2, function(x){sum(as.numeric(x))})

  expect_equal(seq2mat(as.character(1:5), dat)["3", "5"], 8L)
  expect_equal(seq2mat(as.character(1:5), dat)["5", "3"], 8L)
})

# combinadic ####
test_that("combinadic", {
  bcombin <- combn(LETTERS[1:5], 2)
  .combin <- combinadic(LETTERS[1:5], 2, 1)
  expect_true(all.equal(bcombin[, 1], rev(.combin)))
  .combin <- combinadic(LETTERS[1:5], 2, 2)
  expect_true(all.equal(bcombin[, 2], rev(.combin)))
  expect_error(combinadic(LETTERS[1:5], 0, 2), "must be")

})

# weighted ####
test_that("weighted", {
  set.seed(5)
  x <- rnorm(5)
  weights <-  c(0.1, 0.2, 0.3, 0.4, 0.0)
  test <- weighted(x, weights)
  expect_equal(test, -0.15580413)
  expect_error(weighted(x, c(0.1, 0.2)), "match the length")
  expect_error(weighted("a", c(0.1)), "should be numeric")
  expect_error(weighted(0.1, "a"), "should be numeric")
  weights <-  c(0.1, 0.2, 0.3, 0.4, 0.2)
  expect_warning(weighted(x, weights), "the weights is above 1")
})

# duplicateIndices ####
test_that("duplicateIndices", {
  vec <- c("52", "52", "52", "53", "55")
  test <- duplicateIndices(as.character(c(vec, 3)))
  expect_true(is.list(test))
  expect_length(test$`52`, 3L)

  b <- matrix(ncol = 2, nrow = 3)
  expect_error(duplicateIndices(b), "Expected")
})

# removeDup ####
test_that("removeDup", {
  a <- seq2mat(c("52", "52", "53", "55"), rnorm(choose(4, 2)))

  mat <- list("kegg" = a, "react" = a)
  dupli <- duplicateIndices(rownames(a))
  remat <- removeDup(mat, dupli)

  expect_equal(ncol(remat[[1]]), ncol(mat[[1]]) - 1)
  expect_equal(nrow(remat[[1]]), nrow(mat[[1]]) - 1)
  expect_equal(nrow(remat[[1]]), ncol(remat[[1]]))

  b <- matrix(ncol = 2, nrow = 3)
  mat <- list("kegg" = a, "react" = b)
  expect_error(removeDup(mat, dupli), "should be symmetric")
})

# genesInfo ####
test_that("genesInfo", {
  info <- structure(
    list(ENTREZID = c("1", "2", "2", "2", "2", "2", "2", "2",
                      "2", "2", "2", "2", "3", "4", "5", "6",
                      "7", "8", "9", "9", "9", "9", "10", "10",
                      "10", "10"),
         REACTOMEID = c(NA, "109582", "114608", "140837",
                        "140877", "1474228", "1474244",
                        "162582", "194315", "194840", "76002",
                        "76005", NA, NA, NA, NA, NA, NA,
                        "1430728", "156580", "156582", "211859",
                        "1430728", "156580","156582", "211859")
    ),
    .Names = c("ENTREZID", "REACTOMEID"),
    class = "data.frame", row.names = c(NA, -26L))

  test <- genesInfo(info, "REACTOMEID", c("9", "10"), "REACTOMEID")
  expect_equal(test, as.character())

})
