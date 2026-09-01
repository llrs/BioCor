# BioCor

This project wants to allow the user to calculate functional
similarities (or biological correlation as it was named originally hence
the name) and use them for network building or other purposes.

# Installation

It is an R package you can install it from the Bioconductor project
with:

``` r

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install("BioCor")
```

You can install this version of *BioCor* with:

``` r

if (!requireNamespace("devtools", quietly = TRUE)) {
    install.packages("devtools")
}
devtools::install_github("llrs/BioCor")
```

# How to use BioCor?

See the
[vignette](https://bioconductor.org/packages/release/bioc/vignettes/BioCor/inst/doc/BioCor_1_basics.html)
in Bioconductor site and the [advanced
vignette](https://bioconductor.org/packages/release/bioc/vignettes/BioCor/inst/doc/BioCor_2_advanced.html).  
Here is a minimum example:

``` r

set.seed(24)
# The data must be provided, see the vignette for more details.
# Get some pathways from the pathway data
(pathways <- sample(unlist(genesReact, use.names = FALSE), 5))
#> [1] "R-HSA-5689603" "R-HSA-168256"  "R-HSA-5696397" "R-HSA-9665230"
#> [5] "R-HSA-9658195"
# Calculate the pathway similarity of them
mpathSim(pathways, genesReact, NULL)
#>               R-HSA-5689603 R-HSA-168256 R-HSA-5696397 R-HSA-9665230
#> R-HSA-5689603   1.000000000 0.0415194346   0.070796460  0.0000000000
#> R-HSA-168256    0.041519435 1.0000000000   0.003634711  0.0009174312
#> R-HSA-5696397   0.070796460 0.0036347115   1.000000000  0.0000000000
#> R-HSA-9665230   0.000000000 0.0009174312   0.000000000  1.0000000000
#> R-HSA-9658195   0.007936508 0.0871794872   0.000000000  0.0000000000
#>               R-HSA-9658195
#> R-HSA-5689603   0.007936508
#> R-HSA-168256    0.087179487
#> R-HSA-5696397   0.000000000
#> R-HSA-9665230   0.000000000
#> R-HSA-9658195   1.000000000
```

# Who might use this package?

It is intended for bioinformaticians, both people interested in
*knowing* the functionally *similarity of some genes* or clusters and
people *developing* some other analysis at the top of it.

# What is the goal of this project?

The goal of this project is to provide methods to calculate functional
similarities based on pathways.

# What can be BioCor used for?

Here is a non-comprehensive list:

- Diseases or drug:  
  By observing which genes with the same pathways are more affected
- Gene/protein functional analysis:  
  By testing how new pathways are similar to existing pathways
- Protein-protein interaction:  
  By testing if they are involved in the same pathways
- miRNA-mRNA interaction:  
  By comparing clusters they affect
- sRNA regulation:  
  By observing the relationship between sRNA and genes
- Evolution:  
  By comparing similarities of genes between species
- Networks improvement:  
  By adding information about the known relationship between genes
- Evaluate pathways databases:  
  By comparing scores of the same entities

See the [advanced
vignette](https://bioconductor.org/packages/release/bioc/vignettes/BioCor/inst/doc/BioCor_2_advanced.html)

# Contributing

Please read [how to contribute](https://biocor.llrs.dev/CONTRIBUTING.md)
for details on the code of conduct, and the process for submitting pull
requests.

# Acknowledgments

Anyone that has contributed to make this package be as is, specially my
advisor.
