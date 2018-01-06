library("BioCor")
context("Testing list to GeneSetCollection method")

test_that("info", {

    test <- as.GeneSetCollection(info)
    sum <- info[sapply(info, function(x)all(!is.na(x)))]
    expect_true(class(test) == "GeneSetCollection")
    expect_equal(lengths(inverseList(geneIds(test))),
                 lengths(sum))


})
