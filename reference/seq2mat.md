# Transforms a vector to a symmetric matrix

Fills a matrix of `ncol = length(x)` and `nrow = length(x)` with the
values in `dat` and setting the diagonal to 1.

## Usage

``` r
seq2mat(x, dat)
```

## Arguments

- x:

  names of columns and rows, used to define the size of the matrix

- dat:

  Data to fill with the matrix with except the diagonal.

## Value

A square matrix with the diagonal set to 1 and `dat` on the upper and
lower triangle with the columns ids and row ids from x.

## Details

`dat` should be at least `choose(length(x), 2)` of length. It assumes
that the data provided comes from using the row and column id to obtain
it.

## See also

[`upper.tri()`](https://rdrr.io/r/base/lower.tri.html) and
[`lower.tri()`](https://rdrr.io/r/base/lower.tri.html)

## Author

Lluís Revilla

## Examples

``` r
seq2mat(LETTERS[1:5], 1:10)
#>   A B C  D  E
#> A 1 1 2  4  7
#> B 1 1 3  5  8
#> C 2 3 1  6  9
#> D 4 5 6  1 10
#> E 7 8 9 10  1
seq2mat(LETTERS[1:5], seq(from = 0.1, to = 1, by = 0.1))
#>     A   B   C   D   E
#> A 1.0 0.1 0.2 0.4 0.7
#> B 0.1 1.0 0.3 0.5 0.8
#> C 0.2 0.3 1.0 0.6 0.9
#> D 0.4 0.5 0.6 1.0 1.0
#> E 0.7 0.8 0.9 1.0 1.0
```
