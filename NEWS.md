# BioCor 1.28

* Update tests to latest testthat change.

# BioCor 1.26

* Added support for plots of similarities with `plot_data()` and 
  `plot_similarity()`.
* Removed internal data accidentally exposed (It was only used for testing).
* Fix problem with documenting the package via roxygen2.
* Remove the `\describe{}` macro used in some pages.

# Changes in version 1.19.1

* Added alternative text to plots on vignettes.
* Reactivated again the comparison with GOSemSim

# Changes in version 1.12

* Move to NEWS.md.
* Corrected "value of length: 2 type: logical".
* Avoid cropping the images.
* Updated documentation.

# Changes in version 1.10.1

* Added a pkgdown website

# Changes in version 1.8.1

* Exported inverseList function
* Improved spelling
* Corrected a bug related to logical coercion of length greater than 1.

# Changes in version 1.6.1

* Corrected an empty condition in mclusterGeneSim
* Corrected a logical check in mclusterSim
* Improved coherence in the output between lists and GeneSetCollections
* Added links between vignettes

# Changes in version 1.5.4

* Improved spelling
* Reduced the complexity fo combineScoresPar and combineScores.
* Improve efficiency of code in vignettes.
* Change the news section
* Adding a section about GeneOverlap package.
* Changed the License

# Changes in version 1.3.3

* Add methods for GeneSetCollections

# Changes in version 1.1.12

* Improve speed on geneSim and clusterSim
* Add a vignette describing more uses of BioCor
* Adding more tests
* Calculating similarities correctly when non-existent genes are provided.

# Changes in version 1.1.4

* Correct incorrect implementation of rcmax method

# Changes in version 1.1.3


* Improving reciprocal method
* Increased testind of combineScores

# Changes in version 1.1.0

* Adding combineScoresPar for high number of comparisons
* Improving the speed of mpathSim
* Removing the dependency of graph

# Changes in version 0.99.33

* Add more test to increase coverage
* Document better the BiocViews.

# Changes in version 0.99.32

* Accept named pathways in mpathSim
* Accept named genes in mgeneSim
* diceSim returns NA when not able to calculate similarities
* Increase test coverage

# Changes in version 0.99.31

* Refer to weighted.mean in weighted help page
* Change the default behaviour of mpathSim: now methods is NULL
* Correct the misleading use of geneSetsCollection.
* Add test for when large number of pathways are needed (pathSims_matrix)

# Changes in version 0.99.30

* Add reciprocal method to combineScores

# Changes in version 0.99.30

* Update documentation of combineScores

# Changes in version 0.99.29

* Update documentation of combineScores

# Changes in version 0.99.28

* Update documentation of combineScores

# Changes in version 0.99.27

* Correct a test for clusterSim

# Changes in version 0.99.26

* Update the required R version from 3.3 to 3.4's
* Adding more tests to increase the coverage
* Using mpathSim to calculate path similarities

# Changes in version 0.99.25

* Update Vignette to just output html

# Changes in version 0.99.24

* Submit the package to the Bioconductor project

# Changes in version 0.99.23

* Add the Bioconductor webhook
