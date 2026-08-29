# Compare pathways

Function to estimate how much two list of genes overlap by looking how
much of the nodes are shared. Calculates the Dice similarity

## Usage

``` r
diceSim(g1, g2)
```

## Arguments

- g1, g2:

  A character list with the names of the proteins in each pathway.

## Value

A score between 0 and 1 calculated as the double of the proteins shared
by g1 and g2 divided by the number of genes in both groups.

## Details

It requires a vector of characters otherwise will return an `NA`.

## See also

Used for [`geneSim()`](https://biocor.llrs.dev/reference/geneSim.md),
see [`conversions()`](https://biocor.llrs.dev/reference/conversions.md)
help page to transform Dice score to Jaccard score.

## Author

Lluís Revilla

## Examples

``` r
genes.id2 <- c("52", "11342", "80895", "57654", "548953", "11586", "45985")
genes.id1 <- c(
  "52", "11342", "80895", "57654", "58493", "1164", "1163",
  "4150", "2130", "159"
)
diceSim(genes.id1, genes.id2)
#> [1] 0.4705882
diceSim(genes.id2, genes.id2)
#> [1] 1
```
