# Similarity score genes based on pathways similarity

Given two genes, calculates the Dice similarity between each pathway
which is combined to obtain a similarity between the genes.

## Usage

``` r
geneSim(gene1, gene2, info, method = "max", ...)

# S4 method for class 'character,character,GeneSetCollection'
geneSim(gene1, gene2, info, method = "max", ...)
```

## Arguments

- gene1, gene2:

  Ids of the genes to calculate the similarity, to be found in genes.

- info:

  A GeneSetCollection or a list of genes and the pathways they are
  involved.

- method:

  one of `c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")`,
  see Details.

- ...:

  Other arguments passed to `combineScores`

## Value

The highest Dice score of all the combinations of pathways between the
two ids compared if a method to combine scores is provided or NA if
there isn't information for one gene. If an `NA` is returned this means
that there isn't information available for any pathways for one of the
genes. Otherwise a number between 0 and 1 (both included) is returned.
Note that there isn't a negative value of similarity.

## Details

Given the information about the genes and their pathways, uses the ids
of the genes to find the Dice similarity score for each pathway
comparison between the genes. Later this similarities are combined using
[`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md).

## Methods (by class)

- `geneSim(gene1 = character, gene2 = character, info = GeneSetCollection)`:
  Calculates all the similarities of the GeneSetCollection and combine
  them using
  [`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)

## See also

[`mgeneSim()`](https://biocor.llrs.dev/reference/mgeneSim.md),
[`conversions()`](https://biocor.llrs.dev/reference/conversions.md) help
page to transform Dice score to Jaccard score. For the method to combine
the scores see
[`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md).

## Author

Lluís Revilla

## Examples

``` r
if (require("org.Hs.eg.db") & require("reactome.db")) {
  # Extract the paths of all genes of org.Hs.eg.db from KEGG
  # (last update in data of June 31st 2011)
  genes.kegg <- as.list(org.Hs.egPATH)
  # Extracts the paths of all genes of org.Hs.eg.db from reactome
  genes.react <- as.list(reactomeEXTID2PATHID)
  geneSim("81", "18", genes.react)
  geneSim("81", "18", genes.kegg)
  geneSim("81", "18", genes.react, NULL)
  geneSim("81", "18", genes.kegg, NULL)
} else {
  warning("You need reactome.db and org.Hs.eg.db package for this example")
}
#> Loading required package: reactome.db
#>       00250 00280 00410 00640 00650       01100
#> 04510     0     0     0     0     0 0.001503759
#> 04520     0     0     0     0     0 0.000000000
#> 04530     0     0     0     0     0 0.000000000
#> 04670     0     0     0     0     0 0.003210273
#> 04810     0     0     0     0     0 0.004467610
#> 05146     0     0     0     0     0 0.011326861
#> 05322     0     0     0     0     0 0.000000000
#> 05412     0     0     0     0     0 0.000000000
```
