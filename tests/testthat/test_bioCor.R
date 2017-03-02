library("BioCor")
context("Testing BioCor function")

test_that("bioCor react entrez", {
    similarities <- bioCor(genes.id, kegg = FALSE, react = TRUE)
    expect_is(similarities, "list")
    expect_length(similarities, 1L)
    expect_equal(ncol(similarities$react), nrow(similarities$react))
    expect_equal(ncol(similarities$react), 10L)
    expect_true(all(diag(as.matrix(similarities$react)) == 1L))

})

test_that("bioCor react symbols", {
    similarities <- bioCor(genes.symbol,  ids = "SYMBOL", kegg = FALSE,
                           react = TRUE)
    expect_is(similarities, "list")
    expect_length(similarities, 1L)
    expect_equal(ncol(similarities$react), nrow(similarities$react))
    expect_equal(ncol(similarities$react), 10L)
    expect_true(all(diag(as.matrix(similarities$react)) == 1L))

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
    expect_warning(sim1 <- bioCor(genes.id[1:3], kegg = FALSE, react = TRUE))
    sim2 <- bioCor(genes.id, kegg = FALSE, react = TRUE)
    expect_equal(sim1[[1]]["52", "11342"], sim2[[1]]["52", "11342"])
})

test_that("bioCor works whith non existing keys in a table", {
    expect_warning(sim <- bioCor(c("3855", "3880"), react = TRUE))
    expect_warning(bioCor(c("3855", "3880"), react = TRUE))
    expect_warning(expect_message(bioCor(c("3855", "3880"), react = TRUE)))
    expect_true(is.na(sim[[1]][1,2]))
})
