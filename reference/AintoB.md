# Insert a matrix into another

Insert values from a matrix into another matrix based on the rownames
and colnames replacing the values.

## Usage

``` r
AintoB(A, B)
```

## Arguments

- A:

  A matrix to be inserted.

- B:

  A matrix to insert in.

## Value

A matrix with the values of A in the matrix B.

## Details

If all the genes with pathway information are already calculated but you
would like to use more genes when performing analysis. insert the once
you have calculated on the matrix of genes.

## Author

Lluís Revilla

## Examples

``` r
B <- matrix(
  ncol = 10, nrow = 10,
  dimnames = list(letters[1:10], letters[1:10])
)
A <- matrix(c(1:15),
  byrow = TRUE, nrow = 5,
  dimnames = list(letters[1:5], letters[1:3])
)
AintoB(A, B)
#>    a  b  c  d  e  f  g  h  i  j
#> a  1  2  3 NA NA NA NA NA NA NA
#> b  4  5  6 NA NA NA NA NA NA NA
#> c  7  8  9 NA NA NA NA NA NA NA
#> d 10 11 12 NA NA NA NA NA NA NA
#> e 13 14 15 NA NA NA NA NA NA NA
#> f NA NA NA NA NA NA NA NA NA NA
#> g NA NA NA NA NA NA NA NA NA NA
#> h NA NA NA NA NA NA NA NA NA NA
#> i NA NA NA NA NA NA NA NA NA NA
#> j NA NA NA NA NA NA NA NA NA NA

# Mixed orders
colnames(A) <- c("c", "h", "e")
rownames(A) <- c("b", "a", "f", "c", "j")
AintoB(A, B)
#>    a  b  c  d  e  f  g  h  i  j
#> a NA NA  4 NA  6 NA NA  5 NA NA
#> b NA NA  1 NA  3 NA NA  2 NA NA
#> c NA NA 10 NA 12 NA NA 11 NA NA
#> d NA NA NA NA NA NA NA NA NA NA
#> e NA NA NA NA NA NA NA NA NA NA
#> f NA NA  7 NA  9 NA NA  8 NA NA
#> g NA NA NA NA NA NA NA NA NA NA
#> h NA NA NA NA NA NA NA NA NA NA
#> i NA NA NA NA NA NA NA NA NA NA
#> j NA NA 13 NA 15 NA NA 14 NA NA

# Missing colums or rows
colnames(A) <- c("d", "f", "k")
AintoB(A, B)
#>    a  b  c  d  e  f  g  h  i  j
#> a NA NA NA  4 NA  5 NA NA NA NA
#> b NA NA NA  1 NA  2 NA NA NA NA
#> c NA NA NA 10 NA 11 NA NA NA NA
#> d NA NA NA NA NA NA NA NA NA NA
#> e NA NA NA NA NA NA NA NA NA NA
#> f NA NA NA  7 NA  8 NA NA NA NA
#> g NA NA NA NA NA NA NA NA NA NA
#> h NA NA NA NA NA NA NA NA NA NA
#> i NA NA NA NA NA NA NA NA NA NA
#> j NA NA NA 13 NA 14 NA NA NA NA
```
