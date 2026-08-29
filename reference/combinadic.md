# i-th combination of n elements taken from r

Function similar to combn but for larger vectors. To avoid allocating a
big vector with all the combinations each one can be computed with this
function.

## Usage

``` r
combinadic(n, r, i)
```

## Arguments

- n:

  Elements to extract the combination from

- r:

  Number of elements per combination

- i:

  ith combination

## Value

The combination ith of the elements

## References

[StackOverflow answer
4494469/2886003](https://stackoverflow.com/a/4494469/2886003)

## See also

[`combn()`](https://rdrr.io/r/utils/combn.html)

## Author

Joshua Ulrich

## Examples

``` r
# Output of all combinations
combn(LETTERS[1:5], 2)
#>      [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
#> [1,] "A"  "A"  "A"  "A"  "B"  "B"  "B"  "C"  "C"  "D"  
#> [2,] "B"  "C"  "D"  "E"  "C"  "D"  "E"  "D"  "E"  "E"  
# Otuput of the second combination
combinadic(LETTERS[1:5], 2, 2)
#> [1] "C" "A"
```
