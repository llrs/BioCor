# Similarity score between clusters of genes based on pathways similarity

Looks for the similarity between genes in groups. Once the pathways for
each cluster are found they are combined using
[`combineScores()`](https://biocor.llrs.dev/reference/combineScores.md).

## Usage

``` r
mclusterSim(clusters, info, method = "max", ...)

# S4 method for class 'list,GeneSetCollection'
mclusterSim(clusters, info, method = "max", ...)
```

## Arguments

- clusters:

  A list of clusters of genes to be found in `id`.

- info:

  A GeneSetCollection or a list of genes and the pathways they are
  involved.

- method:

  one of `c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")`,
  see Details.

- ...:

  Other arguments passed to `combineScores`

## Value

`mclusterSim` returns a matrix with the similarity scores for each
cluster comparison.

## Methods (by class)

- `mclusterSim(clusters = list, info = GeneSetCollection)`: Calculates
  all the similarities of the GeneSetCollection and combine them using
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

  clusters <- list(
    cluster1 = c("18", "81", "10"),
    cluster2 = c("100", "10", "1"),
    cluster3 = c("18", "10", "83")
  )
  mclusterSim(clusters, genes.kegg)
  mclusterSim(clusters, genes.kegg, "avg")
} else {
  warning("You need org.Hs.eg.db package for this example")
}
#>            cluster1   cluster2  cluster3
#> cluster1 0.11837329 0.07739749 0.1158909
#> cluster2 0.07739749 0.26653339 0.1448643
#> cluster3 0.11589087 0.14486431 0.2183986
```
