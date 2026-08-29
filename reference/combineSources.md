# Combine different sources of pathways

Given several sources of pathways with the same for the same id of the
genes it merge them.

## Usage

``` r
combineSources(...)
```

## Arguments

- ...:

  Lists of genes and their pathways.

## Value

A single list with the pathways of each source on the same gene.

## Details

It assumes that the identifier of the genes are the same for both
sources but if many aren't equal it issues a warning. Only unique
pathways identifiers are returned.

## Examples

``` r
DB1 <- list(g1 = letters[6:8], g2 = letters[1:5], g3 = letters[4:7])
DB2 <- list(
  g1 = c("one", "two"), g2 = c("three", "four"),
  g3 = c("another", "two")
)
combineSources(DB1, DB2)
#> $g1
#> [1] "f"   "g"   "h"   "one" "two"
#> 
#> $g2
#> [1] "a"     "b"     "c"     "d"     "e"     "three" "four" 
#> 
#> $g3
#> [1] "d"       "e"       "f"       "g"       "another" "two"    
#> 
combineSources(DB1, DB1)
#> $g1
#> [1] "f" "g" "h"
#> 
#> $g2
#> [1] "a" "b" "c" "d" "e"
#> 
#> $g3
#> [1] "d" "e" "f" "g"
#> 
DB3 <- list(
  g1 = c("one", "two"), g2 = c("three", "four"),
  g4 = c("five", "six", "seven"), g5 = c("another", "two")
)
combineSources(DB1, DB3) # A warning is expected
#> Warning: More than 25% of genes identifiers of a source are unique
#> Check the identifiers of the genes
#> $g1
#> [1] "f"   "g"   "h"   "one" "two"
#> 
#> $g2
#> [1] "a"     "b"     "c"     "d"     "e"     "three" "four" 
#> 
#> $g3
#> [1] "d" "e" "f" "g"
#> 
#> $g4
#> [1] "five"  "six"   "seven"
#> 
#> $g5
#> [1] "another" "two"    
#> 
```
