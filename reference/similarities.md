# Apply a function to a list of similarities

Function to join list of similarities by a function provided by the
user.

## Usage

``` r
similarities(sim, func, ...)
```

## Arguments

- sim:

  list of similarities to be joined. All similarities must have the same
  dimensions. The genes are assumed to be in the same order for all the
  matrices.

- func:

  function to perform on those similarities: `prod`, `sum`... It should
  accept as many arguments as similarities matrices are provided, and
  should use numbers.

- ...:

  Other arguments passed to the function `func`. Usually `na.rm` or
  similar.

## Value

A matrix of the size of the similarities

## Note

It doesn't check that the columns and rows of the matrices are in the
same order or are the same.

## See also

[`weighted()`](https://biocor.llrs.dev/reference/weighted.md) for
functions that can be used, and
[`addSimilarities()`](https://biocor.llrs.dev/reference/addSimilarities.md)
for a wrapper to one of them

## Author

Lluís Revilla

## Examples

``` r
set.seed(100)
a <- seq2mat(LETTERS[1:5], rnorm(10))
b <- seq2mat(LETTERS[1:5], seq(from = 0.1, to = 1, by = 0.1))
sim <- list(b, a)
similarities(sim, weighted.prod, c(0.5, 0.5))
#>              A            B            C           D           E
#> A  0.250000000 -0.012554809  0.006576558  0.08867848 -0.10181337
#> B -0.012554809  0.250000000 -0.005918782  0.01462141  0.14290654
#> C  0.006576558 -0.005918782  0.250000000  0.04779451 -0.18568337
#> D  0.088678481  0.014621409  0.047794513  0.25000000 -0.08996553
#> E -0.101813370  0.142906542 -0.185683371 -0.08996553  0.25000000
# Note the differences in the sign of some values
similarities(sim, weighted.sum, c(0.5, 0.5))
#>            A          B          C          D          E
#> A  1.0000000 -0.3010962  0.1657656  0.6433924 -0.6408953
#> B -0.3010962  1.0000000 -0.1894585  0.3084856  0.7572664
#> C  0.1657656 -0.1894585  1.0000000  0.4593150 -0.8626297
#> D  0.6433924  0.3084856  0.4593150  1.0000000 -0.6799311
#> E -0.6408953  0.7572664 -0.8626297 -0.6799311  1.0000000
```
