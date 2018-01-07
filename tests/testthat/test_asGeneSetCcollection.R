library("BioCor")
context("Testing list to GeneSetCollection method")

test_that("info", {

    test <- as.GeneSetCollection(info)
    sum <- info[sapply(info, function(x)all(!is.na(x)))]
    expect_true(class(test) == "GeneSetCollection")
    expect_equal(lengths(inverseList(geneIds(test))),
                 lengths(sum))


})


fl <- system.file("extdata", "Broad.xml", package="GSEABase")
gss <- getBroadSets(fl) # GeneSetCollection of 2 sets

test_that("as.GeneSetCollection", {
    l <-  inverseList(geneIds(gss))
    gss2 <- as.GeneSetCollection(l)
    expect_equal(length(geneIds(gss)), length(geneIds(gss2)))
    expect_true(all(names(gss) %in% names(gss2)))
    expect_equal(length(gss), length(gss2))
})
