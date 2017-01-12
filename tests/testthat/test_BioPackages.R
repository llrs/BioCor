library("BioCor")
library("doParallel")
context("Testing BioCor fundamental functions")

# corPathways ####
test_that("corPathways", {
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
    test <- corPathways(c("2", "9"), info, "ENTREZID", "REACTOMEID")
    expect_equal(test, 0L)
    expect_error(corPathways(c("2", "9"), info, "ENTREZID", "REACTOME"),
                 "pathway")
})

# goCor ####
genes.id <- as.character(c(52, 11342, 80895, 57654, 58493, 1164, 1163, 4150,
                           2130, 159))

test_that("go_cor", {

})

# comparePathways ####
test_that("comparePathways", {
    genes.id2 <- c("52", "11342", "80895", "57654", "548953", "11586", "45985")
    test <- comparePathways(genes.id, genes.id2)
    expect_equal(test, 0.47058823)
    graph1 <- graph::graphNEL(nodes = genes.id)
    graph2 <- graph::graphNEL(nodes = genes.id2)

    test2 <- comparePathways(graph1, genes.id2)
    expect_equal(test, test2)
    test3 <- comparePathways(graph1, graph2)
    expect_equal(test2, test3)
    expect_equal(comparePathways(graph1, graph1), 1L)
    expect_equal(comparePathways(c(), c()), 0L)
})

# addSimilarities ####
test_that("addSimilarities", {
    set.seed(1)
    a <- seq2mat(LETTERS[1:5], runif(10))
    b <- seq2mat(LETTERS[1:5], runif(10))
    x <- seq2mat(LETTERS[1:5], runif(10))
    mat <- list("kegg" = a , "react" = b)
    weights <- c(0.1, 0.2, 0.8)
    expect_error(addSimilarities(x, mat, weights), "Weights are too big")
    weights <- c(0.1, 0.2, 0.3)
    expect_warning(addSimilarities(x, mat, weights), "smaller than 1")
    weights <- c(0.3, 0.3, 0.4)
    expect_error(addSimilarities("a", mat, weights), "Expected a matrix")
    test <- addSimilarities(x, mat, weights)

    expect_true(isSymmetric(test))
    expect_true(all(diag(test) == 1))
    expect_true(all(test <= 1))
    expect_true(all(test >= -1))
    expect_true(test["B", "D"] <= 0.6129278)

    x <- matrix(runif(10), ncol = 5)
    expect_error(addSimilarities(x, mat, weights), "is different")
})

# bioCor ####
test_that("bioCor react", {
    similarities <- bioCor(genes.id, kegg = FALSE, react = TRUE)
    expect_is(similarities, "list")
    expect_length(similarities, 1L)
    expect_equal(ncol(similarities$react), nrow(similarities$react))
    expect_equal(ncol(similarities$react), 10L)
})

test_that("bioCor kegg", {
    similarities <- bioCor(genes.id, kegg = TRUE, react = FALSE)
    expect_message(bioCor(genes.id, kegg = TRUE, react = FALSE))
    expect_is(similarities, "list")
    expect_length(similarities, 1L)
    expect_equal(ncol(similarities$kegg), nrow(similarities$kegg))
    expect_equal(ncol(similarities$kegg), 10L)
})

test_that("bioCor is replicable", {
    registerDoParallel(2)
    sim1 <- bioCor(genes.id, kegg = FALSE, react = TRUE)
    stopImplicitCluster()
    registerDoParallel(2)
    sim2 <- bioCor(genes.id, kegg = FALSE, react = TRUE)
    stopImplicitCluster()
    expect_equal(sim1[[1]]["1163", "159"], sim2[[1]]["1163", "159"])
})

test_that("bioCor works whith non existing keys", {
    sim <- bioCor(c("3855", "3880"), react = TRUE)
    expect_warning(bioCor(c("3855", "3880"), react = TRUE))
    expect_message(bioCor(c("3855", "3880"), react = TRUE))
    expect_true(is.na(sim[[1]][1,2]))
})

# test_that("bioCor go", {
#   similarities <- bioCor(genes, kegg = FALSE, react = FALSE, go = TRUE)
# })
