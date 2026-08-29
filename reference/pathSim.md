# Calculates the Dice similarity between pathways

Calculates the similarity between pathways using dice similarity score.
`diceSim` is used to calculate similarities between the two pathways.

## Usage

``` r
pathSim(pathway1, pathway2, info)

# S4 method for class 'character,character,GeneSetCollection'
pathSim(pathway1, pathway2, info)
```

## Arguments

- pathway1, pathway2:

  A single pathway to calculate the similarity

- info:

  A GeneSetCollection or a list of genes and the pathways they are
  involved.

## Value

The similarity between those pathways or all the similarities between
each comparison.

## Methods (by class)

- `pathSim(pathway1 = character, pathway2 = character, info = GeneSetCollection)`:
  Calculates all the similarities of a GeneSetCollection and combine
  them using `combineScoresPar`

## See also

[`conversions()`](https://biocor.llrs.dev/reference/conversions.md) help
page to transform Dice score to Jaccard score.
[`mpathSim()`](https://biocor.llrs.dev/reference/mpathSim.md) for
multiple pairwise comparison of pathways.

## Author

Lluís Revilla

## Examples

``` r
if (require("reactome.db")) {
  # Extracts the paths of all genes of org.Hs.eg.db from reactome
  genes.react <- as.list(reactomeEXTID2PATHID)
  (paths <- sample(unique(unlist(genes.react)), 2))
  pathSim(paths[1], paths[2], genes.react)
} else {
  warning("You need reactome.db package for this example")
}
#> [1] 0
```
