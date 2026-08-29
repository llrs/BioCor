# Combining values

Combine several similarities into one using several methods.

## Usage

``` r
combineScores(
  scores,
  method = c("max", "avg", "rcmax", "rcmax.avg", "BMA", "reciprocal"),
  round = FALSE,
  t = 0
)

combineScoresPar(scores, method, subSets = NULL, BPPARAM = NULL, ...)
```

## Arguments

- scores:

  Matrix of scores to be combined

- method:

  one of `c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal")`,
  see Details.

- round:

  Should the resulting value be rounded to the third digit?

- t:

  Numeric value to filter scores below this value. Only used in the
  reciprocal method.

- subSets:

  List of combinations as info in other functions.

- BPPARAM:

  BiocParallel back-end parameters. By default (`NULL`) a `for` loop is
  used.

- ...:

  Other arguments passed to `combineScores`

## Value

A numeric value as described in details.

## Details

The input matrix can be a base matrix or a matrix from package Matrix.
The methods return:

- **avg**: The average or mean value.

- **max**: The max value.

- **rcmax**: The max of the column means or row means.

- **rcmax.avg**: The sum of the max values by rows and columns divided
  by the number of columns and rows.

- **BMA**: The same as `rcmax.avg`.

- **reciprocal**: The double of the sum of the reciprocal maximal
  similarities (above a threshold) divided by the number of elements.
  See equation 3 of the Tao *et al* 2007 article.

## Note

`combineScores` is a version of the function of the same name in package
GOSemSim
([`GOSemSim::combineScores()`](https://rdrr.io/pkg/GOSemSim/man/combineScores.html))
with optional rounding and some internal differences.

## References

Ying Tao, Lee Sam, Jianrong Li, Carol Friedman, Yves A. Lussier;
Information theory applied to the sparse gene ontology annotation
network to predict novel gene function. Bioinformatics 2007; 23 (13):
i529-i538. doi: 10.1093/bioinformatics/btm195

## See also

[register](https://rdrr.io/pkg/BiocParallel/man/register.html) in
BiocParallel about the arguments accepted by BPPARAM.

## Author

Lluís Revilla based on Guangchuang Yu.

## Examples

``` r
(d <- structure(
  c(
    0.4, 0.6, 0.222222222222222, 0.4, 0.4, 0, 0.25, 0.5,
    0.285714285714286
  ),
  .Dim = c(3L, 3L),
  .Dimnames = list(c("a", "b", "c"), c("d", "e", "f"))
))
#>           d   e         f
#> a 0.4000000 0.4 0.2500000
#> b 0.6000000 0.4 0.5000000
#> c 0.2222222 0.0 0.2857143
e <- d
sapply(c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal"),
  combineScores,
  scores = d
)
#>        avg        max      rcmax  rcmax.avg        BMA reciprocal 
#>  0.3397707  0.6000000  0.5000000  0.4642857  0.4642857  0.3333333 
d[1, 2] <- NA
sapply(c("avg", "max", "rcmax", "rcmax.avg", "BMA", "reciprocal"),
  combineScores,
  scores = d
)
#>        avg        max      rcmax  rcmax.avg        BMA reciprocal 
#>  0.3322421  0.6000000  0.5000000  0.4642857  0.4642857  0.2000000 
colnames(e) <- rownames(e)
combineScoresPar(e, list(a = c("a", "b"), b = c("b", "c")),
  method = "max"
)
#>     a   b
#> a 0.6 0.5
#> b 0.5 0.5
```
