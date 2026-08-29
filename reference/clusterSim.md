# Similarity score between clusters of genes based on pathways similarity

Looks for the similarity between genes in groups

## Usage

``` r
clusterSim(cluster1, cluster2, info, method = "max", ...)

# S4 method for class 'character,character,GeneSetCollection'
clusterSim(cluster1, cluster2, info, method = "max", ...)
```

## Arguments

- cluster1, cluster2:

  A vector with genes.

- info:

  A GeneSetCollection or a list of genes and the pathways they are
  involved.

- method:

  one of `c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")`,
  see Details.

- ...:

  Other arguments passed to `combineScores`

## Value

`clusterSim` returns a similarity score of the two clusters

## Details

Once the pathways for each cluster are found they are combined using
[`combineScores()`](https://biocor.llrs.dev/reference/combineScores.md).

## Methods (by class)

- `clusterSim( cluster1 = character, cluster2 = character, info = GeneSetCollection )`:
  Calculates all the similarities of the GeneSetCollection and combine
  them using
  [`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)

## See also

For a different approach see
[`clusterGeneSim()`](https://biocor.llrs.dev/reference/clusterGeneSim.md),
[`combineScores()`](https://biocor.llrs.dev/reference/combineScores.md)
and [`conversions()`](https://biocor.llrs.dev/reference/conversions.md)

## Author

Lluís Revilla

## Examples

``` r
if (require("org.Hs.eg.db")) {
  # Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
  # data of June 31st 2011)
  genes.kegg <- as.list(org.Hs.egPATH)
  clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg)
  clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg, NULL)
  clusterSim(c("9", "15", "10"), c("33", "19", "20"), genes.kegg, "avg")
} else {
  warning("You need org.Hs.eg.db package for this example")
}
#> [1] 0.07083146
```
