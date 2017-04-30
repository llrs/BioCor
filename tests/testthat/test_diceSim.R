library("BioCor")
context("Testing diceSim")

test_that("diceSim", {
    genes.id2 <- c("52", "11342", "80895", "57654", "548953", "11586", "45985")
    genes.id <- c("52", "11342", "80895", "57654", "58493", "1164", "1163",
                  "4150", "2130", "159")
    test <- diceSim(genes.id, genes.id2)
    expect_equal(test, 0.47058823)
    graph1 <- graph::graphNEL(nodes = genes.id)
    graph2 <- graph::graphNEL(nodes = genes.id2)

    test2 <- diceSim(graph1, genes.id2)
    expect_equal(test, test2)
    test3 <- diceSim(graph1, graph2)
    expect_equal(test2, test3)
    test4 <- diceSim(genes.id, graph2)
    expect_equal(test2, test4)
    expect_equal(diceSim(graph1, graph1), 1L)
    expect_warning(diceSim(c(), c()), "is not")
})
