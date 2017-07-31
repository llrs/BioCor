[![Build Status](https://travis-ci.org/llrs/BioCor.svg?branch=master)](
https://travis-ci.org/llrs/BioCor) [![Build Status](http://www.bioconductor.org/shields/build/devel/bioc/BioCor.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/BioCor/)

[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/BioCor.svg)](https://www.bioconductor.org/packages/devel/bioc/html/BioCor.html#since)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![codecov](https://codecov.io/gh/llrs/BioCor/branch/master/graph/badge.svg)](https://codecov.io/gh/llrs/BioCor/)  [![commit](http://www.bioconductor.org/shields/commits/bioc/BioCor.svg)](https://www.bioconductor.org/packages/devel/bioc/html/BioCor.html#svn_source)

# BioCor package

This project wants to allow the user to calculate functional similarities (or biological correlation as it was named originally hence the name) and 
use them for network building or other purposes.

# How does it work?

It is an R package you can install it from the Bioconductor project with:

```r
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("BioCor")
```
You can install this version of *BioCor* with:
```r
library("devtools")
install_github("llrs/BioCor")
```

# Who will use this repo or project?

It is intended for bioinformaticians, both people interested in *knowing* the functionally similarity of some genes or clusters and people *developing* some other analysis at the top of it.

# How to use BioCor?

See the [vignette](http://bioconductor.org/packages/release/bioc/vignettes/BioCor/inst/doc/BioCor.html) in Bioconductor site and the [advanced vignette](https://llrs.github.io/vignette2.html).

# What is the goal of this project?

The goal of this project is to provide methods to calculate functional similarities based on pathways. 

# What can be BioCor used for?

 - Diseases or drug:  
  By observing which genes with the same pathways are more affected
 - Gene/protein functional analysis:  
  By testing how new pathways
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

See [advanced vignette](https://llrs.github.io/vignette2.html)

# Contributing

Please read [how to contribute](.github/CONTRIBUTING.md) for details on the code of conduct, and the process for submitting pull requests.

# Acknowledgments

Anyone that has contributed to make this package be as is, including, my advisor.
