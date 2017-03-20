[![Build Status](https://travis-ci.org/llrs/BioCor.svg?branch=master)](
https://travis-ci.org/llrs/BioCor)  


[![Bioc](http://www.bioconductor.org/shields/years-in-bioc/BioCor.svg)](https://www.bioconductor.org/packages/devel/bioc/html/BioCor.html#since) [![total](https://img.shields.io/badge/downloads-43487/total-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/BioCor) [![month](https://img.shields.io/badge/downloads-1663/month-blue.svg?style=flat)](https://bioconductor.org/packages/stats/bioc/BioCor)

[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![codecov](https://codecov.io/gh/llrs/BioCor/branch/master/graph/badge.svg)](https://codecov.io/gh/llrs/BioCor/)  [![commit](http://www.bioconductor.org/shields/commits/bioc/BioCor.svg)](https://www.bioconductor.org/packages/devel/bioc/html/BioCor.html#svn_source)

# BioCor package

This project wants to allow the user to calculate functional similarities and 
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

# What is the goal of this project?

The goal of this project is to provide methods to calculate functional similarities based on pathways. 

# Contributing

Please read [CONTRIBUTING.md](.github/CONTRIBUTING.md) for details on our code of conduct, and the process for submitting pull requests.

# Acknowledgments

Anyone that has contributed to make this package be as is, including, my advisor.
