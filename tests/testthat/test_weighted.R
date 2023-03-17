library("BioCor")
context("Testing weighted functions")

test_that("weighted.sum", {
    set.seed(5)
    x <- runif(5)
    weights <- c(0.1, 0.2, 0.3, 0.4, 0.0)
    test <- weighted.sum(c(-0.5, 0.5), c(0.6, 0.2))
    expect_equal(test, -0.4)
    test <- weighted.sum(c(-0.5, 0.5), c(0.6, 0.2), abs = FALSE)
    expect_equal(test, -0.2)
    test <- weighted.sum(x, weights)
    expect_equal(test, 0.54588768)
    expect_error(weighted.sum(x, c(0.1, 0.2)), "match the length")
    expect_error(weighted.sum("a", c(0.1)), "should be numeric")
    expect_error(weighted.sum(0.1, "a"), "should be numeric")
    weights <- c(0.1, 0.2, 0.3, 0.4, 0.2)
    expect_warning(weighted.sum(x, weights), "the weights is above 1")
    expect_warning(
        weighted.sum(c(-0.5, 0.5), c(0.6, -0.2), abs = FALSE),
        "negative"
    )
})

test_that("weighted.prod", {
    set.seed(5)
    x <- runif(5)
    weights <- c(0.1, 0.2, 0.3, 0.4, NA)
    test <- weighted.prod(x, weights)
    expect_equal(weighted.prod(x, c(0.1, 0.2, 0.3, 0.4, 0)), 0L)
    expect_equal(test, 8.58568731153938e-05)
    expect_error(weighted.prod(x, c(0.1, 0.2)), "match the length")
    expect_error(weighted.prod("a", c(0.1)), "should be numeric")
    expect_error(weighted.prod(0.1, "a"), "should be numeric")
    weights <- c(0.1, 0.2, 0.3, 0.4, 0.2)
    expect_warning(weighted.prod(x, weights), "the weights is above 1")
})
