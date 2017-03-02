library("BioCor")
context("Testing conversions functions")

test_that("Conversions", {
    expect_error(D2J(1.1))
    expect_error(D2J(-1.1))
    expect_error(J2D(1.1))
    expect_error(J2D(-1.1))

    expect_equal(J2D(D2J(0.5)), 0.5)

    expect_equal(D2J(0.5), 1/3)
    expect_equal(J2D(0.5), 2/3)
})
