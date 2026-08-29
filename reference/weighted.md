# Weighted operations

Calculates the weighted sum or product of `x`. Each values should have
its weight, otherwise it will throw an error.

## Usage

``` r
weighted.sum(x, w, abs = TRUE)

weighted.prod(x, w)
```

## Arguments

- x:

  an object containing the values whose weighted operations is to be
  computed

- w:

  a numerical vector of weights the same length as `x` giving the
  weights to use for elements of `x`.

- abs:

  If any `x` is negative you want the result negative too?

## Value

`weighted.sum` returns the sum of the product of x\*weights removing all
`NA` values. See parameter `abs` if there are any negative values.

`weighted.prod` returns the product of product of x\*weights removing
all `NA` values.

## Details

This functions are thought to be used with `similarities`. As some
similarities might be positive and others negative the argument `abs` is
provided for `weighted.sum`, assuming that only one similarity will be
negative (usually the one coming from expression correlation).

## See also

[`weighted.mean()`](https://rdrr.io/r/stats/weighted.mean.html),
[`similarities()`](https://biocor.llrs.dev/reference/similarities.md)
and
[`addSimilarities()`](https://biocor.llrs.dev/reference/addSimilarities.md)

## Author

Lluís Revilla

## Examples

``` r
expr <- c(-0.2, 0.3, 0.5, 0.8, 0.1)
weighted.sum(expr, c(0.5, 0.2, 0.1, 0.1, 0.1))
#> [1] -0.3
weighted.sum(expr, c(0.5, 0.2, 0.1, 0.2, 0.1), FALSE)
#> Warning: The sum of the weights is above 1
#> [1] 0.18
weighted.sum(expr, c(0.4, 0.2, 0.1, 0.2, 0.1))
#> [1] -0.36
weighted.sum(expr, c(0.4, 0.2, 0.1, 0.2, 0.1), FALSE)
#> [1] 0.2
weighted.sum(expr, c(0.4, 0.2, 0, 0.2, 0.1))
#> [1] -0.31
weighted.sum(expr, c(0.5, 0.2, 0, 0.2, 0.1))
#> [1] -0.33
# Compared to weighted.prod:
weighted.prod(expr, c(0.5, 0.2, 0.1, 0.1, 0.1))
#> [1] -2.4e-07
weighted.prod(expr, c(0.4, 0.2, 0.1, 0.2, 0.1))
#> [1] -3.84e-07
weighted.prod(expr, c(0.4, 0.2, 0, 0.2, 0.1))
#> [1] 0
weighted.prod(expr, c(0.5, 0.2, 0, 0.2, 0.1))
#> [1] 0
```
