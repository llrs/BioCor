library("BioCor")
context("Testing auxiliar functions")

# seq2mat ####
test_that("seq2mat puts a combination to the right place", {
    test <- seq2mat(LETTERS[1:5], 1:10)

    expect_true(isSymmetric(test)) # already checks if it is squared.

    expect_error(seq2mat(LETTERS[1:5], 1:11))
    expect_error(seq2mat(LETTERS[1:5], 1:9))

    com <- combn(as.character(1:5), 2)
    dat <- apply(com, 2, function(x) {
        sum(as.numeric(x))
    })

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
    a <- seq2mat(c("52", "52", "53", "55"), runif(choose(4, 2)))

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


test_that("inverseList", {
    l <- list("a" = NA, "b" = letters[1:4])
    expect_equal(inverseList(l), list(a = "b", b = "b", c = "b", d = "b"),
        .Names = c("a", "b", "c", "d")
    )
})
