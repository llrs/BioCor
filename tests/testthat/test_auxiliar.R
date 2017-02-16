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
test_that("weighted.sum", {
    set.seed(5)
    x <- runif(5)
    weights <-  c(0.1, 0.2, 0.3, 0.4, 0.0)
    test <- weighted.sum(c(-0.5, 0.5), c(0.6, 0.2))
    expect_equal(test, -0.4)
    test <- weighted.sum(c(-0.5, 0.5), c(0.6, 0.2), abs = FALSE)
    expect_equal(test, -0.2)
    test <- weighted.sum(x, weights)
    expect_equal(test, 0.54588768)
    expect_error(weighted.sum(x, c(0.1, 0.2)), "match the length")
    expect_error(weighted.sum("a", c(0.1)), "should be numeric")
    expect_error(weighted.sum(0.1, "a"), "should be numeric")
    weights <-  c(0.1, 0.2, 0.3, 0.4, 0.2)
    expect_warning(weighted.sum(x, weights), "the weights is above 1")
})

test_that("weighted.prod", {
    set.seed(5)
    x <- runif(5)
    weights <-  c(0.1, 0.2, 0.3, 0.4, NA)
    test <- weighted.prod(x, weights)
    expect_equal(weighted.prod(x, c(0.1, 0.2, 0.3, 0.4, 0)), 0L)
    expect_equal(test, 8.58568731153938e-05)
    expect_error(weighted.prod(x, c(0.1, 0.2)), "match the length")
    expect_error(weighted.prod("a", c(0.1)), "should be numeric")
    expect_error(weighted.prod(0.1, "a"), "should be numeric")
    weights <-  c(0.1, 0.2, 0.3, 0.4, 0.2)
    expect_warning(weighted.prod(x, weights), "the weights is above 1")
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
    expect_is(test, "list")
    expect_is(test[[1L]], "character")
    expect_equal(names(test), c("9", "10"))
    expect_length(test[[1L]], 0L)
})

# D2J and J2D ####
test_that("Conversions", {
    expect_error(D2J(1.1))
    expect_error(D2J(-1.1))
    expect_error(J2D(1.1))
    expect_error(J2D(-1.1))

    expect_equal(J2D(D2J(0.5)), 0.5)

    expect_equal(D2J(0.5), 1/3)
    expect_equal(J2D(0.5), 2/3)
})
# combineScores ####
test_that("combineScores", {
    d <- structure(c(0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
                     0.285714285714286), .Dim = c(3L, 3L),
                   .Dimnames = list(c("a","b", "c"), c("d", "e", "f")))
    test <- sapply(c("avg", "max", "rcmax", "rcmax.avg", "BMA"), combineScores,
           scores = d)
    expect_equal(as.numeric(test["rcmax.avg"]), as.numeric(test["BMA"]))
    expect_true(all(test["max"] > test[-2]))
    expect_equal(as.numeric(test["avg"]), 0.339770723104056)
    expect_equal(as.numeric(test["max"]), 0.6)
    expect_equal(as.numeric(test["rcmax"]), 0.5)
    expect_equal(as.numeric(test["rcmax.avg"]), 0.464285714285714)

})
