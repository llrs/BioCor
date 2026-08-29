# Remove duplicated rows and columns

Given the indices of the duplicated entries remove the columns and rows
until just one is left, it keeps the duplicated with the highest
absolute mean value.

## Usage

``` r
removeDup(cor_mat, dupli)
```

## Arguments

- cor_mat:

  List of matrices

- dupli:

  List of indices with duplicated entries

## Value

A matrix with only one of the columns and rows duplicated

## See also

[`duplicateIndices()`](https://biocor.llrs.dev/reference/duplicateIndices.md)
to obtain the list of indices with duplicated entries.

## Author

Lluís Revilla

## Examples

``` r
a <- seq2mat(c("52", "52", "53", "55"), runif(choose(4, 2)))
b <- seq2mat(c("52", "52", "53", "55"), runif(choose(4, 2)))
mat <- list("kegg" = a, "react" = b)
mat
#> $kegg
#>           52        52        53        55
#> 52 1.0000000 0.1325060 0.6171091 0.3380487
#> 52 0.1325060 1.0000000 0.7913045 0.9054476
#> 53 0.6171091 0.7913045 1.0000000 0.1975567
#> 55 0.3380487 0.9054476 0.1975567 1.0000000
#> 
#> $react
#>           52         52        53         55
#> 52 1.0000000 0.79408518 0.7546029 0.32268770
#> 52 0.7940852 1.00000000 0.9113948 0.08616947
#> 53 0.7546029 0.91139485 1.0000000 0.91120045
#> 55 0.3226877 0.08616947 0.9112005 1.00000000
#> 
dupli <- duplicateIndices(rownames(a))
remat <- removeDup(mat, dupli)
remat
#> $kegg
#>           52        53        55
#> 52 1.0000000 0.7913045 0.9054476
#> 53 0.7913045 1.0000000 0.1975567
#> 55 0.9054476 0.1975567 1.0000000
#> 
#> $react
#>           52        53        55
#> 52 1.0000000 0.7546029 0.3226877
#> 53 0.7546029 1.0000000 0.9112005
#> 55 0.3226877 0.9112005 1.0000000
#> 
```
