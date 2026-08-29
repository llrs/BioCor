# Similarity score between clusters of genes based on genes similarity

Looks for the similarity between genes of a group and then between each
group's genes.

## Usage

``` r
clusterGeneSim(cluster1, cluster2, info, method = c("max", "rcmax.avg"), ...)

# S4 method for class 'character,character,GeneSetCollection'
clusterGeneSim(cluster1, cluster2, info, method = c("max", "rcmax.avg"), ...)
```

## Arguments

- cluster1, cluster2:

  A vector with genes.

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

Returns a similarity score between the genes of the two clusters.

## Details

Differs with clusterSim that first each combination between genes is
calculated, and with this values then the comparison between the two
clusters is done. Thus applying combineScores twice, one at gene level
and another one at cluster level.

## Methods (by class)

- `clusterGeneSim( cluster1 = character, cluster2 = character, info = GeneSetCollection )`:
  Calculates the gene similarities in a GeneSetCollection and combine
  them using
  [`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)

## See also

[`mclusterGeneSim()`](https://biocor.llrs.dev/reference/mclusterGeneSim.md),
[`combineScores()`](https://biocor.llrs.dev/reference/combineScores.md)
and [`clusterSim()`](https://biocor.llrs.dev/reference/clusterSim.md)

## Author

Lluís Revilla

## Examples

``` r
if (require("org.Hs.eg.db")) {
  # Extract the paths of all genes of org.Hs.eg.db from KEGG (last update in
  # data of June 31st 2011)
  genes.kegg <- as.list(org.Hs.egPATH)
  clusterGeneSim(c("18", "81", "10"), c("100", "10", "1"), genes.kegg)
  clusterGeneSim(
    c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
    c("avg", "avg")
  )
  clusterGeneSim(
    c("18", "81", "10"), c("100", "10", "1"), genes.kegg,
    c("avg", "rcmax.avg")
  )
  (clus <- clusterGeneSim(
    c("18", "81", "10"), c("100", "10", "1"),
    genes.kegg, "avg"
  ))
  combineScores(clus, "rcmax.avg")
} else {
  warning("You need org.Hs.eg.db package for this example")
}
#> Loading required package: org.Hs.eg.db
#> Loading required package: AnnotationDbi
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: IRanges
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> 
#> [1] 0.2061464
```
