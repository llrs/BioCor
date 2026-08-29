# Additive integration of similarities

Function that use the previously calculated similarities into a single
similarity matrix.

## Usage

``` r
addSimilarities(x, bio_mat, weights = c(0.5, 0.18, 0.1, 0.22))
```

## Arguments

- x:

  A matrix with the similarity of expression

- bio_mat:

  A list of matrices of the same dimension as x.

- weights:

  A numeric vector of weight to multiply each similarity

## Value

A square matrix of the same dimensions as the input matrices.

## Details

The total weight can't be higher than 1 to prevent values above 1 but
can be below 1. It uses weighted.sum with abs = TRUE internally.

## See also

[`similarities()`](https://biocor.llrs.dev/reference/similarities.md),
[`weighted()`](https://biocor.llrs.dev/reference/weighted.md).

## Author

Lluís Revilla

## Examples

``` r
set.seed(100)
a <- seq2mat(LETTERS[1:5], rnorm(10))
b <- seq2mat(LETTERS[1:5], seq(from = 0.1, to = 1, by = 0.1))
sim <- list(b)
addSimilarities(a, sim, c(0.5, 0.5))
#>            A          B          C          D          E
#> A  1.0000000 -0.3010962  0.1657656  0.6433924 -0.6408953
#> B -0.3010962  1.0000000 -0.1894585  0.3084856  0.7572664
#> C  0.1657656 -0.1894585  1.0000000  0.4593150 -0.8626297
#> D  0.6433924  0.3084856  0.4593150  1.0000000 -0.6799311
#> E -0.6408953  0.7572664 -0.8626297 -0.6799311  1.0000000
```
