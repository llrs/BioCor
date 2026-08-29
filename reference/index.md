# Package index

## General introduction

- [`BioCor`](https://biocor.llrs.dev/reference/BioCor-package.md)
  [`BioCor-package`](https://biocor.llrs.dev/reference/BioCor-package.md)
  : BioCor: A package to calculate functional similarities

## Similarity of pathways

Compare pairs of pathways, just one (pathSim) or as many as you want
(mpathSim)

- [`pathSim()`](https://biocor.llrs.dev/reference/pathSim.md) :
  Calculates the Dice similarity between pathways
- [`mpathSim()`](https://biocor.llrs.dev/reference/mpathSim.md) :
  Calculates the Dice similarity between pathways

## Similarity of genes

Compare pairs of genes, just one (geneSim) or as many as you want
(mgeneSim)

- [`geneSim()`](https://biocor.llrs.dev/reference/geneSim.md) :
  Similarity score genes based on pathways similarity
- [`mgeneSim()`](https://biocor.llrs.dev/reference/mgeneSim.md) :
  Similarity score genes based on pathways similarity

## Pooling similarity of groups of genes

Compare groups of genes by pooling all the pathways they are in

- [`clusterSim()`](https://biocor.llrs.dev/reference/clusterSim.md) :
  Similarity score between clusters of genes based on pathways
  similarity
- [`mclusterSim()`](https://biocor.llrs.dev/reference/mclusterSim.md) :
  Similarity score between clusters of genes based on pathways
  similarity

## Similarity of groups of genes

Compare groups of genes by comparing the similarity of the genes in each
group

- [`clusterGeneSim()`](https://biocor.llrs.dev/reference/clusterGeneSim.md)
  : Similarity score between clusters of genes based on genes similarity
- [`mclusterGeneSim()`](https://biocor.llrs.dev/reference/mclusterGeneSim.md)
  : Similarity score between clusters of genes based on genes similarity

## Convert between similarities

Convert from Jaccard to Dice or vice versa

- [`D2J()`](https://biocor.llrs.dev/reference/conversions.md)
  [`J2D()`](https://biocor.llrs.dev/reference/conversions.md) : Convert
  the similarities formats

## Combine similarities

- [`combineScores()`](https://biocor.llrs.dev/reference/combineScores.md)
  [`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)
  : Combining values

## Helper functions

- [`plot_data()`](https://biocor.llrs.dev/reference/plot_similarity.md)
  [`plot_similarity()`](https://biocor.llrs.dev/reference/plot_similarity.md)
  : The position of the nodes is based on the similarity between them.
- [`combineSources()`](https://biocor.llrs.dev/reference/combineSources.md)
  : Combine different sources of pathways
- [`AintoB()`](https://biocor.llrs.dev/reference/AintoB.md) : Insert a
  matrix into another
- [`seq2mat()`](https://biocor.llrs.dev/reference/seq2mat.md) :
  Transforms a vector to a symmetric matrix
- [`combinadic()`](https://biocor.llrs.dev/reference/combinadic.md) :
  i-th combination of n elements taken from r
- [`diceSim()`](https://biocor.llrs.dev/reference/diceSim.md) : Compare
  pathways
- [`duplicateIndices()`](https://biocor.llrs.dev/reference/duplicateIndices.md)
  : Finds the indices of the duplicated events of a vector
- [`removeDup()`](https://biocor.llrs.dev/reference/removeDup.md) :
  Remove duplicated rows and columns
- [`inverseList()`](https://biocor.llrs.dev/reference/inverseList.md) :
  Invert a list
- [`similarities()`](https://biocor.llrs.dev/reference/similarities.md)
  : Apply a function to a list of similarities
- [`addSimilarities()`](https://biocor.llrs.dev/reference/addSimilarities.md)
  : Additive integration of similarities
- [`weighted.sum()`](https://biocor.llrs.dev/reference/weighted.md)
  [`weighted.prod()`](https://biocor.llrs.dev/reference/weighted.md) :
  Weighted operations
