test_that("combineSources", {
    DB1 <- list(g1 = letters[6:8], g2 = letters[1:5], g3 = letters[4:7])
    DB2 <- list(
        g1 = c("one", "two"), g2 = c("three", "four"),
        g3 = c("another", "two")
    )
    test <- combineSources(DB1, DB2)
    expect_length(test, 3)
    expect_equal(lengths(test), c(g1 = 5L, g2 = 7L, g3 = 6L))
    expect_equal(combineSources(DB1, DB1), DB1)
    DB3 <- list(
        g1 = c("one", "two"), g2 = c("three", "four"),
        g4 = c("five", "six", "seven"), g5 = c("another", "two")
    )
    expect_warning(combineSources(DB1, DB3), "25%")
    DB4 <- list(
        g1 = c("one", "two"), g2 = c("three", "four"),
        g4 = c("five", "six", "seven"), g5 = c("another", "two"),
        g6 = c("a", "b"), g7 = c("test", "bla")
    )
    expect_warning(combineSources(DB1, DB4), "50%")
})
