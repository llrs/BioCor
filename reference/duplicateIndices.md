# Finds the indices of the duplicated events of a vector

Finds the indices of duplicated elements in the vector given.

## Usage

``` r
duplicateIndices(vec)
```

## Arguments

- vec:

  Vector of identifiers presumably duplicated

## Value

The format is determined by the simplify2array

## Details

For each duplication it can return a list or if all the duplication
events are of the same length it returns a matrix, where each column is
duplicated.

## See also

[`removeDup()`](https://biocor.llrs.dev/reference/removeDup.md)

## Author

Lluís Revilla

## Examples

``` r
duplicateIndices(c("52", "52", "53", "55")) # One repeated element
#> $`52`
#> [1] 1 2
#> 
duplicateIndices(c("52", "52", "53", "55", "55")) # Repeated elements
#> $`52`
#> [1] 1 2
#> 
#> $`55`
#> [1] 4 5
#> 
duplicateIndices(c("52", "55", "53", "55", "52")) # Mixed repeated elements
#> $`55`
#> [1] 2 4
#> 
#> $`52`
#> [1] 1 5
#> 
```
