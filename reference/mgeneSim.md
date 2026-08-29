# Similarity score genes based on pathways similarity

Given two genes, calculates the Dice similarity between each pathway
which is combined to obtain a similarity between the genes.

## Usage

``` r
mgeneSim(genes, info, method = "max", ...)

# S4 method for class 'character,GeneSetCollection'
mgeneSim(genes, info, method = "max", ...)

# S4 method for class 'missing,GeneSetCollection'
mgeneSim(genes, info, method = "max", ...)
```

## Arguments

- genes:

  A vector of genes.

- info:

  A GeneSetCollection or a list of genes and the pathways they are
  involved.

- method:

  one of `c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")`,
  see Details.

- ...:

  Other arguments passed to `combineScores`

## Value

`mgeneSim` returns the matrix of similarities between the genes in the
vector

## Details

Given the information about the genes and their pathways, uses the ids
of the genes to find the Dice similarity score for each pathway
comparison between the genes. Later this similarities are combined using
[`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md).

## Methods (by class)

- `mgeneSim(genes = character, info = GeneSetCollection)`: Calculates
  all the similarities of the list and combine them using
  [`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)

- `mgeneSim(genes = missing, info = GeneSetCollection)`: Calculates all
  the similarities of the list and combine them using
  [`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)

## Note

genes accept named characters and the output will use the names of the
genes.

## See also

[`geneSim()`](https://biocor.llrs.dev/reference/geneSim.md),
[`conversions()`](https://biocor.llrs.dev/reference/conversions.md) help
page to transform Dice score to Jaccard score. For the method to combine
the scores see
[`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md).

## Examples

``` r
if (require("org.Hs.eg.db") & require("reactome.db")) {
  # Extract the paths of all genes of org.Hs.eg.db from KEGG
  # (last update in data of June 31st 2011)
  genes.kegg <- as.list(org.Hs.egPATH)
  # Extracts the paths of all genes of org.Hs.eg.db from reactome
  genes.react <- as.list(reactomeEXTID2PATHID)
  mgeneSim(c("81", "18", "10"), genes.react)
  mgeneSim(c("81", "18", "10"), genes.react, "avg")
  named_genes <- structure(c("81", "18", "10"),
    .Names = c("ACTN4", "ABAT", "NAT2")
  )
  mgeneSim(named_genes, genes.react, "max")
} else {
  warning("You need reactome.db and org.Hs.eg.db package for this example")
}
#>            ACTN4       ABAT       NAT2
#> ACTN4 1.00000000 0.12379110 0.06050641
#> ABAT  0.12379110 1.00000000 0.05112045
#> NAT2  0.06050641 0.05112045 1.00000000
```
