# Calculates the Dice similarity between pathways

Calculates the similarity between several pathways using dice similarity
score. If one needs the matrix of similarities between pathways set the
argument methods to `NULL`.

## Usage

``` r
mpathSim(pathways, info, method = NULL, ...)

# S4 method for class 'character,GeneSetCollection,ANY'
mpathSim(pathways, info, method = NULL, ...)

# S4 method for class 'missing,GeneSetCollection,ANY'
mpathSim(pathways, info, method = NULL, ...)

# S4 method for class 'missing,list,ANY'
mpathSim(pathways, info, method = NULL, ...)

# S4 method for class 'missing,list,missing'
mpathSim(pathways, info, method = NULL, ...)
```

## Arguments

- pathways:

  Pathways to calculate the similarity for

- info:

  A list of genes and the pathways they are involved or a
  GeneSetCollection object

- method:

  To combine the scores of each pathway, one of
  `c("avg", "max", "rcmax", "rcmax.avg", "BMA")`, if NULL returns the
  matrix of similarities.

- ...:

  Other arguments passed to
  [`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)

## Value

The similarity between those pathways or all the similarities between
each comparison.

## Methods (by class)

- `mpathSim(pathways = character, info = GeneSetCollection, method = ANY)`:
  Calculates the similarity between the provided pathways of the
  GeneSetCollection using `combineScoresPar`

- `mpathSim(pathways = missing, info = GeneSetCollection, method = ANY)`:
  Calculates all the similarities of the GeneSetCollection and combine
  them using `combineScoresPar`

- `mpathSim(pathways = missing, info = list, method = ANY)`: Calculates
  all the similarities of the list and combine them using
  `combineScoresPar`

- `mpathSim(pathways = missing, info = list, method = missing)`:
  Calculates all the similarities of the list

## Note

`pathways` accept named characters, and then the output will have the
names

## See also

[`pathSim()`](https://biocor.llrs.dev/reference/pathSim.md) For single
pairwise comparison.
[`conversions()`](https://biocor.llrs.dev/reference/conversions.md) To
convert the Dice similarity to Jaccard similarity

## Examples

``` r
if (require("reactome.db")) {
  genes.react <- as.list(reactomeEXTID2PATHID)
  (pathways <- sample(unique(unlist(genes.react)), 10))
  mpathSim(pathways, genes.react, NULL)
  named_paths <- structure(
    c("R-HSA-112310", "R-HSA-112316", "R-HSA-112315"),
    .Names = c(
      "Neurotransmitter Release Cycle",
      "Neuronal System",
      "Transmission across Chemical Synapses"
    )
  )
  mpathSim(named_paths, genes.react, NULL)
  many_pathways <- sample(unique(unlist(genes.react)), 152)
  mpathSim(many_pathways, genes.react, "avg")
} else {
  warning("You need reactome.db package for this example")
}
#> [1] 0.007304779
```
