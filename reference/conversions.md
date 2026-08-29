# Convert the similarities formats

Functions to convert the similarity coefficients between Jaccard and
Dice. D2J is the opposite of J2D.

## Usage

``` r
D2J(D)

J2D(J)
```

## Arguments

- D:

  Dice coefficient, as returned by
  [`diceSim()`](https://biocor.llrs.dev/reference/diceSim.md),
  [`geneSim()`](https://biocor.llrs.dev/reference/geneSim.md),
  [`clusterSim()`](https://biocor.llrs.dev/reference/clusterSim.md) and
  [`clusterGeneSim()`](https://biocor.llrs.dev/reference/clusterGeneSim.md)

- J:

  Jaccard coefficient

## Value

A numeric value.

## Author

Lluís Revilla

## Examples

``` r
D2J(0.5)
#> [1] 0.3333333
J2D(0.5)
#> [1] 0.6666667
D2J(J2D(0.5))
#> [1] 0.5
```
