# Similarity score between clusters of genes based on genes similarity

Looks for the similarity between genes of a group and then between each
group's genes.

## Usage

``` r
mclusterGeneSim(clusters, info, method = c("max", "rcmax.avg"), ...)

# S4 method for class 'list,GeneSetCollection'
mclusterGeneSim(clusters, info, method = c("max", "rcmax.avg"), ...)
```

## Arguments

- clusters:

  A list of clusters of genes to be found in `id`.

- info:

  A GeneSetCollection or a list of genes and the pathways they are
  involved.

- method:

  A vector with two or one argument to be passed to combineScores the
  first one is used to summarize the similarities of genes, the second
  one for clusters.

- ...:

  Other arguments passed to `combineScores`

## Value

Returns a matrix with the similarity scores for each cluster comparison.

## Methods (by class)

- `mclusterGeneSim(clusters = list, info = GeneSetCollection)`:
  Calculates all the similarities of the GeneSetCollection and combine
  them using
  [`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)

## See also

[`clusterGeneSim()`](https://biocor.llrs.dev/reference/clusterGeneSim.md),
[`clusterSim()`](https://biocor.llrs.dev/reference/clusterSim.md) and
[`combineScores()`](https://biocor.llrs.dev/reference/combineScores.md)

## Author

Lluís Revilla

## Examples

``` r
if (require("org.Hs.eg.db")) {
  genes.kegg <- as.list(org.Hs.egPATH)
  clusters <- list(
    cluster1 = c("18", "81", "10"),
    cluster2 = c("100", "594", "836"),
    cluster3 = c("18", "10", "83")
  )
  mclusterGeneSim(clusters, genes.kegg)
  mclusterGeneSim(clusters, genes.kegg, c("max", "avg"))
  mclusterGeneSim(clusters, genes.kegg, c("max", "BMA"))
} else {
  warning("You need org.Hs.eg.db package for this example")
}
#>           cluster1  cluster2  cluster3
#> cluster1 1.0000000 1.0000000 0.8022654
#> cluster2 1.0000000 1.0000000 0.8304646
#> cluster3 0.8022654 0.8304646 1.0000000
```
