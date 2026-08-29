# About BioCor

Abstract

Describes the background of the package, important functions defined in
the package and some of the applications and usages, including the
integration with other packages and analysis, and comparisons with
related packages. Some frequent or important questions about the package
are answered at the end of the document. For more advanced usage you can
look at the other vignette.

## Introduction

Methods to find similarities have been developed for several purposes,
being Jaccard and Dice similarities the most known. In bioinformatics
much of the research on the topic is centered around [Gene
Ontologies](https://www.geneontology.org/) because they provide
controlled vocabularies, as part of their mission:

> The mission of the GO Consortium is to develop an up-to-date,
> comprehensive, computational model of biological systems, from the
> molecular level to larger pathways, cellular and organism-level
> systems.

However, there is another resource of similarities between genes:
metabolic pathways. Metabolic pathways describe the relationship between
genes, proteins, lipids and other elements of the cells. A pathway
describes, to some extent, the function in which it is involved in the
cell. There exists several databases about which gene belong to which
pathway. Together with pathways, gene sets related to a function or to a
phenotype are a source of information of the genes function. With this
package we provide the methods to calculate functional similarities
based on this information.

Here we provides functions to calculate *functional similarities
distances* for pathways, gene sets, genes and clusters of genes. The
name BioCor stands from biological correlation, shortened to BioCor,
because as said we look if some genes are in the same pathways or gene
sets as other genes.

BioCor is different from
*[GeneOverlap](https://bioconductor.org/packages/3.23/GeneOverlap)*
because here we use the Dice index instead of the Jaccard index
(although we provide a function to change from one to the other, see
[this section](#D2J))and that package only allows to compare pathways
but not genes or groups of genes. But
*[GeneOverlap](https://bioconductor.org/packages/3.23/GeneOverlap)*
provides some functionalities to plot the similarity scores and provides
the associated p-value to the comparison of pathways.

The development of this package aimed initially to improve clustering of
genes by functionality in weighted gene co-expression networks using
*[WGCNA](https://CRAN.R-project.org/package=WGCNA)*. The package has
some functions to combine similarities in order to integrate with
`WGCNA`. For other uses you can check the [advanced
vignette.](https://bioconductor.org/packages/3.23/BioCor/vignettes/BioCor_2_advanced.html).

## Citation

You can cite the package as:

``` r

citation("BioCor")
```

## Installation

The BioCor package is available at
[Bioconductor](https://bioconductor.org) and can be downloaded and
installed via BiocManager:

``` r

install.packages("BiocManager")
BiocManager::install("BioCor")
```

You can install the latest version of
*[BioCor](https://github.com/llrs/BioCor)* from Github with:

``` r

library("devtools")
install_github("llrs/BioCor")
```

## Using BioCor

### Preparation

We can load the package and prepare the data for which we want to
calculate the similarities:

``` r

library("BioCor")
## Load libraries with the data of the pathways
library("org.Hs.eg.db")
library("reactome.db")
genesKegg <- as.list(org.Hs.egPATH)
genesReact <- as.list(reactomeEXTID2PATHID)
# Remove genes and pathways which are not from human pathways 
genesReact <- lapply(genesReact, function(x){
    unique(grep("R-HSA-", x, value = TRUE))
    })
genesReact <- genesReact[lengths(genesReact) >= 1] 
```

To avoid having biased data it is important to have all the data about
the pathways and genes associated to all pathways for the organism under
study. Here we assume that we are interested in human pathways. We use
this two databases KEGG and Reactome as they are easy to obtain the
data. However KEGG database is no longer free for large retrievals
therefore it is not longer updated in the Bioconductor annotation
packages.

However, one can use any list where the names of the list are the genes
and the elements of the list the pathways or groups where the gene
belong. One could also read from a GMT file or use GeneSetCollections in
addition or instead of those associations from a pathway database and
convert it to list using:

``` r

library("GSEABase")
paths2Genes <- geneIds(getGmt("/path/to/file.symbol.gmt",
                 geneIdType=SymbolIdentifier()))

genes <- unlist(paths2Genes, use.names = FALSE)
pathways <- rep(names(paths2Genes), lengths(paths2Genes))
genes2paths <- split(pathways, genes) # List of genes and the gene sets
```

With `genes2paths` we have the information ready to use.

### Pathway similarities

We can compute similarities (Dice similarity, see [question 1](#FAQ1) of
FAQ) between two pathways or between several pathways and combine them,
or not:

``` r

(paths <- sample(unique(unlist(genesReact)), 2))
## [1] "R-HSA-5663020" "R-HSA-5654716"
pathSim(paths[1], paths[2], genesReact)
## [1] 0

(pathways <- sample(unique(unlist(genesReact)), 10))
##  [1] "R-HSA-9732724" "R-HSA-9856532" "R-HSA-9031525" "R-HSA-5339700"
##  [5] "R-HSA-9635644" "R-HSA-964739"  "R-HSA-9702600" "R-HSA-9694614"
##  [9] "R-HSA-8877330" "R-HSA-211945"
mpathSim(pathways, genesReact)
##               R-HSA-9732724 R-HSA-9856532 R-HSA-9031525 R-HSA-5339700
## R-HSA-9732724     1.0000000    0.00000000    0.00000000             0
## R-HSA-9856532     0.0000000    1.00000000    0.00000000             0
## R-HSA-9031525     0.0000000    0.00000000    1.00000000             0
## R-HSA-5339700     0.0000000    0.00000000    0.00000000             1
## R-HSA-9635644     0.0000000    0.00000000    0.00000000             0
## R-HSA-964739      0.0000000    0.00000000    0.00000000             0
## R-HSA-9702600     0.0000000    0.00000000    0.00000000             0
## R-HSA-9694614     0.0000000    0.04878049    0.00000000             0
## R-HSA-8877330     0.1111111    0.00000000    0.00000000             0
## R-HSA-211945      0.0000000    0.00000000    0.01801802             0
##               R-HSA-9635644 R-HSA-964739 R-HSA-9702600 R-HSA-9694614
## R-HSA-9732724             0            0             0    0.00000000
## R-HSA-9856532             0            0             0    0.04878049
## R-HSA-9031525             0            0             0    0.00000000
## R-HSA-5339700             0            0             0    0.00000000
## R-HSA-9635644             1            0             0    0.00000000
## R-HSA-964739              0            1             0    0.00000000
## R-HSA-9702600             0            0             1    0.00000000
## R-HSA-9694614             0            0             0    1.00000000
## R-HSA-8877330             0            0             0    0.00000000
## R-HSA-211945              0            0             0    0.00000000
##               R-HSA-8877330 R-HSA-211945
## R-HSA-9732724     0.1111111   0.00000000
## R-HSA-9856532     0.0000000   0.00000000
## R-HSA-9031525     0.0000000   0.01801802
## R-HSA-5339700     0.0000000   0.00000000
## R-HSA-9635644     0.0000000   0.00000000
## R-HSA-964739      0.0000000   0.00000000
## R-HSA-9702600     0.0000000   0.00000000
## R-HSA-9694614     0.0000000   0.00000000
## R-HSA-8877330     1.0000000   0.00000000
## R-HSA-211945      0.0000000   1.00000000
```

When the method to combine the similarities is set to `NULL`
[`mpathSim()`](https://biocor.llrs.dev/reference/mpathSim.md) returns a
matrix of pathway similarities, otherwise it combines the values. In the
next section we can see the methods to combine pathway similarities.

#### Combining values

To combine values we provide a function with several methods:

``` r

sim <- mpathSim(pathways, genesReact)
methodsCombineScores <- c("avg", "max", "rcmax", "rcmax.avg", "BMA",
                          "reciprocal")
sapply(methodsCombineScores, BioCor::combineScores, scores = sim)
##        avg        max      rcmax  rcmax.avg        BMA reciprocal 
##  0.1035582  1.0000000  1.0000000  1.0000000  1.0000000  1.0000000
```

We can also specify the method to combine the similarities in
[`mpathSim()`](https://biocor.llrs.dev/reference/mpathSim.md),
[`geneSim()`](https://biocor.llrs.dev/reference/geneSim.md),
[`mgeneSim()`](https://biocor.llrs.dev/reference/mgeneSim.md),
[`clusterSim()`](https://biocor.llrs.dev/reference/clusterSim.md),
[`mclusterSim()`](https://biocor.llrs.dev/reference/mclusterSim.md),
[`clusterGeneSim()`](https://biocor.llrs.dev/reference/clusterGeneSim.md)
and
[`mclusterGeneSim()`](https://biocor.llrs.dev/reference/mclusterGeneSim.md),
argument method. By default the method is set to “max” to combine
pathways (except in mpathSim where the default is to show all the
pathway similarities) and “BMA” to combine similarities of genes or for
cluster analysis. This function is adapted from
*[GOSemSim](https://bioconductor.org/packages/3.23/GOSemSim)* package.

The function
[`combineScoresPar()`](https://biocor.llrs.dev/reference/combineScores.md)
allows to use a parallel background (using
*[BiocParallel](https://bioconductor.org/packages/3.23/BiocParallel)*)
to combine the scores. It is recommended to use a parallel background if
you calculate more than 300 gene similarities. It also have an argument
in case you want to calculate the similarity scores of several sets.

### Gene similarities

To compare the function of two genes there is the `geneSim` function and
`mgeneSim` function for several comparisons. In this example we compare
the genes BRCA1 and BRCA2 and NAT2, which are the genes 672, 675 and 10
respectively in ENTREZID:

``` r

geneSim("672", "675", genesKegg)
## [1] 0.0824295
geneSim("672", "675", genesReact)
## [1] 1

mgeneSim(c("BRCA1" = "672", "BRCA2" = "675", "NAT2" = "10"), genesKegg)
##           BRCA1       BRCA2        NAT2
## BRCA1 1.0000000 0.082429501 0.000000000
## BRCA2 0.0824295 1.000000000 0.008241758
## NAT2  0.0000000 0.008241758 1.000000000
mgeneSim(c("BRCA1" = "672", "BRCA2" = "675", "NAT2" = "10"), genesReact)
##           BRCA1     BRCA2      NAT2
## BRCA1 1.0000000 1.0000000 0.2171848
## BRCA2 1.0000000 1.0000000 0.2171848
## NAT2  0.2171848 0.2171848 1.0000000
```

Note that for the same genes each database or list provided has
different annotations, which result on different similarity scores. In
this example BRCA1 has 3 and 27 pathways in KEGG and Reactome
respectively and BRCA2 has 1 and 59 pathways in KEGG and Reactome
respectively which results on different scores.

### Gene cluster similarities

There are two methods:

- Combining all the pathways for each cluster and compare between them.
- Calculate the similarity between genes of a cluster and the other
  cluster.

#### By pathways

As explained, in this method all the pathways of a cluster are compared
with all the pathways of the other cluster. If a method to combine
pathways similarities is not provided, all pathway similarities are
returned:

``` r

clusterSim(c("672", "675"), c("100", "10", "1"), genesKegg)
## [1] 0.04210526
clusterSim(c("672", "675"), c("100", "10", "1"), genesKegg, NULL)
##            00230       01100       05340 00232 00983
## 04120 0.00000000 0.000000000 0.011764706     0     0
## 03440 0.04210526 0.006908463 0.000000000     0     0
## 05200 0.00000000 0.008241758 0.005540166     0     0
## 05212 0.00000000 0.001666667 0.019047619     0     0

clusters <- list(cluster1 = c("672", "675"),
                 cluster2 = c("100", "10", "1"),
                 cluster3 = c("18", "10", "83"))
mclusterSim(clusters, genesKegg, "rcmax.avg")
##            cluster1   cluster2   cluster3
## cluster1 1.00000000 0.01587957 0.00256157
## cluster2 0.01587957 1.00000000 0.56412591
## cluster3 0.00256157 0.56412591 1.00000000
mclusterSim(clusters, genesKegg, "max")
##             cluster1   cluster2    cluster3
## cluster1 1.000000000 0.04210526 0.008241758
## cluster2 0.042105263 1.00000000 1.000000000
## cluster3 0.008241758 1.00000000 1.000000000
```

#### By genes

In this method first the similarities between each gene is calculated,
then the similarity between each group of genes is calculated. Requiring
two methods to combine values, the first one to combine pathways
similarities and the second one to combine genes similarities. If only
one is provided it returns the matrix of similarities of the genes of
each cluster:

``` r

clusterGeneSim(c("672", "675"), c("100", "10", "1"), genesKegg)
## [1] 0.02605425
clusterGeneSim(c("672", "675"), c("100", "10", "1"), genesKegg, "max")
##            100          10  1
## 672 0.01176471 0.000000000 NA
## 675 0.04210526 0.008241758 NA

mclusterGeneSim(clusters, genesKegg, c("max", "rcmax.avg"))
##             cluster1   cluster2    cluster3
## cluster1 1.000000000 0.02605425 0.006181319
## cluster2 0.026054248 1.00000000 1.000000000
## cluster3 0.006181319 1.00000000 1.000000000
mclusterGeneSim(clusters, genesKegg, c("max", "max"))
##             cluster1   cluster2    cluster3
## cluster1 1.000000000 0.04210526 0.008241758
## cluster2 0.042105263 1.00000000 1.000000000
## cluster3 0.008241758 1.00000000 1.000000000
```

Note the differences between
[`mclusterGeneSim()`](https://biocor.llrs.dev/reference/mclusterGeneSim.md)
and [`mclusterSim()`](https://biocor.llrs.dev/reference/mclusterSim.md)
in the similarity values of the clusters. If we set
`method = c("max", "max")` in
[`mclusterGeneSim()`](https://biocor.llrs.dev/reference/mclusterGeneSim.md)
then the similarity between the clusters is the same as
[`clusterSim()`](https://biocor.llrs.dev/reference/clusterSim.md).

### Converting similarities

If needed, Jaccard similarity can be calculated from Dice similarity
using [`D2J()`](https://biocor.llrs.dev/reference/conversions.md):

``` r

D2J(sim)
##               R-HSA-9732724 R-HSA-9856532 R-HSA-9031525 R-HSA-5339700
## R-HSA-9732724    1.00000000         0.000   0.000000000             0
## R-HSA-9856532    0.00000000         1.000   0.000000000             0
## R-HSA-9031525    0.00000000         0.000   1.000000000             0
## R-HSA-5339700    0.00000000         0.000   0.000000000             1
## R-HSA-9635644    0.00000000         0.000   0.000000000             0
## R-HSA-964739     0.00000000         0.000   0.000000000             0
## R-HSA-9702600    0.00000000         0.000   0.000000000             0
## R-HSA-9694614    0.00000000         0.025   0.000000000             0
## R-HSA-8877330    0.05882353         0.000   0.000000000             0
## R-HSA-211945     0.00000000         0.000   0.009090909             0
##               R-HSA-9635644 R-HSA-964739 R-HSA-9702600 R-HSA-9694614
## R-HSA-9732724             0            0             0         0.000
## R-HSA-9856532             0            0             0         0.025
## R-HSA-9031525             0            0             0         0.000
## R-HSA-5339700             0            0             0         0.000
## R-HSA-9635644             1            0             0         0.000
## R-HSA-964739              0            1             0         0.000
## R-HSA-9702600             0            0             1         0.000
## R-HSA-9694614             0            0             0         1.000
## R-HSA-8877330             0            0             0         0.000
## R-HSA-211945              0            0             0         0.000
##               R-HSA-8877330 R-HSA-211945
## R-HSA-9732724    0.05882353  0.000000000
## R-HSA-9856532    0.00000000  0.000000000
## R-HSA-9031525    0.00000000  0.009090909
## R-HSA-5339700    0.00000000  0.000000000
## R-HSA-9635644    0.00000000  0.000000000
## R-HSA-964739     0.00000000  0.000000000
## R-HSA-9702600    0.00000000  0.000000000
## R-HSA-9694614    0.00000000  0.000000000
## R-HSA-8877330    1.00000000  0.000000000
## R-HSA-211945     0.00000000  1.000000000
```

Also if one has a Jaccard similarity and wants a Dice similarity, can
use the [`J2D()`](https://biocor.llrs.dev/reference/conversions.md)
function.

## High volumes of gene similarities

We can compute the whole similarity of genes in KEGG or Reactome by
using :

``` r

## Omit those genes without a pathway
nas <- sapply(genesKegg, function(y){all(is.na(y)) | is.null(y)})
genesKegg2 <- genesKegg[!nas]
m <- mgeneSim(names(genesKegg2), genesKegg2, method  = "max")
```

It takes around 5 hours in one core but it requires high memory
available.

If one doesn’t have such a memory available can compute the similarities
by pieces, and then fit it in another matrix with:

``` r

sim <- AintoB(m, B)
```

Usually B is a matrix of size `length(genes)`, see `?AintoB()`.

## An example of usage

In this example I show how to use BioCor to analyse a list of genes by
functionality. With a list of genes we are going to see how similar are
those genes:

``` r

genes.id <- c("10", "15", "16", "18", "2", "9", "52", "3855", "3880", "644", 
              "81327", "9128", "2073", "2893", "5142", "60", "210", "81", 
              "1352", "88", "672", "675")
genes.id <- mapIds(org.Hs.eg.db, keys = genes.id, keytype = "ENTREZID", 
                   column = "SYMBOL")
## 'select()' returned 1:1 mapping between keys and columns
genes <- names(genes.id)
names(genes) <- genes.id
react <- mgeneSim(genes, genesReact)
## Warning in mgeneSim(genes, genesReact): Some genes are not in the list
## provided.
## We remove genes which are not in list (hence the warning):
nan <- genes %in% names(genesReact)
react <- react[nan, nan]
hc <- hclust(as.dist(1 - react))
plot(hc, main = "Similarities between genes")
```

![Dendrogram of the similarities of genes according to
Reactome.](BioCor_1_basics_files/figure-html/hclust1-1.png)

Gene clustering by similarities

Now we can see the relationship between the genes. We can group them for
a cluster analysis to visualize the relationship between the clusters:

``` r

mycl <- cutree(hc, h = 0.2)
clusters <- split(genes[nan], as.factor(mycl))
# Removing clusters of just one gene
(clusters <- clusters[lengths(clusters) >= 2])
## $`1`
##   NAT2  AANAT   NAT1  BLVRA   ALAD  COX10 
##   "10"   "15"    "9"  "644"  "210" "1352" 
## 
## $`2`
## AARS1  ACTB BRCA1 
##  "16"  "60" "672" 
## 
## $`3`
##   ABAT  GRIA4  ACTN2 
##   "18" "2893"   "88" 
## 
## $`4`
##   A2M ACTN4 
##   "2"  "81" 
## 
## $`5`
##   KRT7  KRT19 
## "3855" "3880" 
## 
## $`8`
##  ERCC5  BRCA2 
## "2073"  "675"
names(clusters) <- paste0("cluster", names(clusters))
## Remember we can use two methods to compare clusters
sim_clus1 <- mclusterSim(clusters, genesReact)
plot(hclust(as.dist(1 - sim_clus1)), 
     main = "Similarities between clusters by pathways")
```

![Dendrogram of clusters of genes according to
Reactome.](BioCor_1_basics_files/figure-html/hclust3-1.png)

Clustering using clusterSim

``` r

sim_clus2 <- mclusterGeneSim(clusters, genesReact)
plot(hclust(as.dist(1 - sim_clus2)), 
     main ="Similarities between clusters by genes")
```

![Dendrogram of clusters according to similarities between genes from
Reactome pathways.](BioCor_1_basics_files/figure-html/hclust3b-1.png)

Clustering using clusterGeneSim

Each method results in a different dendrogram as we can see on Figure
@ref(fig:hclust3) compared to Figure @ref(fig:hclust3b).

## Comparing with GOSemSim

In this section I will compare the functional similarity of BioCor with
the closely related package
*[GOSemSim](https://bioconductor.org/packages/3.23/GOSemSim)*. The genes
and gene clusters used were extracted from GOSemSim’s vignette, we only
change the ontology, instead of the molecular function, the biological
process will be used:

``` r

hsGO <- GOSemSim::godata('org.Hs.eg.db', ont = "BP", computeIC = FALSE)
## 
## Warning in GOSemSim::godata("org.Hs.eg.db", ont = "BP", computeIC = FALSE): use
## 'annoDb' instead of 'OrgDb'
## preparing gene to GO mapping data...
```

I will compare the functions geneSim from section [geneSim and mgeneSim
from
GOSemSim](https://bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html#genesim-and-mgenesim)
with both data sets from KEGG and Reactome:

``` r

goSemSim <- GOSemSim::geneSim("241", "251", semData = hsGO, 
                              measure = "Wang", combine="BMA")
# In case it is null
sim <- ifelse(is.na(goSemSim), 0, getElement(goSemSim, "geneSim"))
BioCor::geneSim("241", "251", genesReact, "BMA") - sim
## [1] 0.04819402

genes <- c("835", "5261","241", "994")
goSemSim <- GOSemSim::mgeneSim(genes, semData = hsGO, 
                   measure = "Wang", combine = "BMA",
                   verbose = FALSE, drop = NULL)
BioCor::mgeneSim(genes, genesReact, "BMA", round = TRUE) - goSemSim
##         835   5261    241    994
## 835   0.000 -0.341 -0.160 -0.182
## 5261 -0.341  0.000 -0.021 -0.345
## 241  -0.160 -0.021  0.000 -0.149
## 994  -0.182 -0.345 -0.149  0.000
```

We can observe there is more similarity according to the gene ontology
than according to the pathways. See FAQ [question 8](#conflict) about
the use of `BioCor::` and `GOSemSim::`.

If named characters are passed they are used to name the resulting
matrix:

``` r

genes <- c("CDC45", "MCM10", "CDC20", "NMU", "MMP1")
genese <- mapIds(org.Hs.eg.db, keys = genes, column = "ENTREZID", 
                 keytype = "SYMBOL")
## 'select()' returned 1:1 mapping between keys and columns
BioCor::mgeneSim(genese, genesReact, "BMA")
##            CDC45      MCM10     CDC20        NMU       MMP1
## CDC45 1.00000000 0.90721577 0.4870304 0.05000484 0.08838455
## MCM10 0.90721577 1.00000000 0.4691837 0.05921898 0.09515408
## CDC20 0.48703038 0.46918373 1.0000000 0.16108096 0.21760465
## NMU   0.05000484 0.05921898 0.1610810 1.00000000 0.11932506
## MMP1  0.08838455 0.09515408 0.2176046 0.11932506 1.00000000
```

We can further compare the cluster similarities from the [next section
of the
vignette](https://bioconductor.org/packages/release/bioc/vignettes/GOSemSim/inst/doc/GOSemSim.html#clustersim-and-mclustersim):

``` r

gs1 <- c("835", "5261","241", "994", "514", "533")
gs2 <- c("578","582", "400", "409", "411")
BioCor::clusterSim(gs1, gs2, genesReact, "BMA") - 
    GOSemSim::clusterSim(gs1, gs2, hsGO, measure = "Wang", combine = "BMA")
## [1] -0.2172804

x <- org.Hs.egGO
hsEG <- mappedkeys(x)
set.seed(123)
(clusters <- list(a=sample(hsEG, 20), b=sample(hsEG, 20), c=sample(hsEG, 20)))
## $a
##  [1] "417"    "5555"   "4540"   "2805"   "5108"   "57820"  "7005"   "10095" 
##  [9] "283951" "4100"   "80201"  "28795"  "54529"  "90459"  "4291"   "9356"  
## [17] "10099"  "51088"  "79879"  "84318" 
## 
## $b
##  [1] "2346"   "390245" "116987" "153562" "29761"  "51657"  "28585"  "159371"
##  [9] "442868" "4597"   "143379" "23220"  "5950"   "85417"  "25778"  "348235"
## [17] "55636"  "9415"   "693135" "157773"
## 
## $c
##  [1] "100861412" "57213"     "80021"     "374739"    "25953"     "91754"    
##  [7] "53"        "125981"    "10945"     "10794"     "256764"    "55901"    
## [13] "5276"      "10090"     "57515"     "26127"     "51751"     "54726"    
## [19] "79581"     "100302289"
BioCor::mclusterSim(clusters, genesReact, "BMA") - 
    GOSemSim::mclusterSim(clusters, hsGO, measure = "Wang", combine = "BMA")
## Warning in BioCor::mclusterSim(clusters, genesReact, "BMA"): Some genes are not
## in the list provided.
##            a          b          c
## a 0.00000000 0.09764821 0.13582351
## b 0.09764821 0.00000000 0.07851081
## c 0.13582351 0.07851081 0.00000000
```

## WGCNA and BioCor

*[WGCNA](https://CRAN.R-project.org/package=WGCNA)* uses the correlation
of the expression data of several samples to cluster genes. Sometimes,
from a biological point of view the interpretation of the resulting
modules is difficult, even more when some groups of genes end up not
having an enrichment in previously described functions. BioCor was
originally thought to be used to overcome this problem: to help
clustering genes, not only by correlation but also by functionality.

In order to have groups functionally related, functional similarities
can enhance the clustering of genes when combined with experimental
correlations. The resulting groups will reflect, not only the
correlation of the expression provided, but also the functionality known
of those genes.

We propose the following steps:

1.  Calculate the similarities for the expression data
2.  Calculate the similarities of the genes in the expression
3.  Combine the similarities
4.  Calculate the adjacency
5.  Identify modules with hierarchical clustering

Here we provide an example on how to use BioCor with WGCNA:

`sim` is a list where each element is a matrix of similarities between
genes Our normalized expression is in the `expr` variable, a matrix
where the samples are in the rows and genes in the columns.

``` r

expr.sim <- WGCNA::cor(expr) # or bicor

## Combine the similarities
similarity <- similarities(c(list(exp = expr.sim), sim), mean, na.rm = TRUE)

## Choose the softThreshold
pSFT <- pickSoftThreshold.fromSimilarity(similarity)

## Or any other function we want
adjacency <- adjacency.fromSimilarity(similarity, power = pSFT$powerEstimate)

## Once we have the similarities we can calculate the TOM with TOM
TOM <- TOMsimilarity(adjacency) ## Requires adjacencies despite its name 
dissTOM <- 1 - TOM
geneTree <- hclust(as.dist(dissTOM), method = "average")
## We can use a clustering tool to group the genes
dynamicMods <- cutreeHybrid(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 30)
moduleColors <- labels2colors(dynamicMods$labels)
```

Once the modules are identified using the functional similarities of
this package and the gene correlations, one can continue with the
workflow of WGCNA.

An important aspect in this process is deciding how to combine the
similarities and the expression data: - If the functional similarities
play a huge role, we will end up having only those genes closely related
to the same functions. - If the functional similarities play a low role,
it will be similarly to only use WGCNA, and the genes won’t be
functionally related.

For these reasons it is better to use weights between `0.5` and `1` for
expression if you use `weighted.sum` or similar functions.

There are several things to take into account when choosing a way to
combine: - The size of the gray or 0 modules (those who don’t show a
specific pattern) - The number and size of the modules created. - The
way the similarities are combined

Violin plots may help to view the differences in size and distribution
of the modules across different methods of combining the similarities.

## FAQ

### How is defined the pathway similarity?

BioCor uses the [Sørensen–Dice
index](https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient):
The dice similarity is the double of the genes shared by the pathways
divided by the number of genes in each pathway.

We can calculate the similarity between two pathways ($`x`$, $`w`$)
with:

``` math
Dice(x, w) = \frac{2 |x \cap w|}{|x| + |w|}
```

This is implemented in the `diceSim` function, which results is similar
to Jaccard index:

``` math
Jaccard(x, w) = \frac{|x \cap w|}{|x \cup w|}
```

Both Jaccard index and dice index are between 0 and 1 ($`[0, 1]`$). To
calculate the Jaccard index from the `diceSim` use the `D2J` function.

### Why does BioCor use the dice coefficient and not the Jaccard ?

We consider Dice coefficient better than Jaccard because it has higher
values for the same comparisons, which reflects that including a gene in
a pathway is not easily done.

### How does BioCor combine similarities between several pathways of two genes?

Although the recommend method is the “max” method, (set as default),
there are implemented other methods in `combineScores` of the
*[GOSemSim](https://bioconductor.org/packages/3.23/GOSemSim)* package
which I borrowed[^1].

### Why do you recommend using the max method to combine similarities scores for pathways?

The purpose of combining the scores is usually to find the relationships
between genes through their pathways. The higher the similarity is
between two pathway of two genes, the higher functionality do the genes
share, even if those genes have other non-related functions.

### How to detect which functional relationship is more important between two genes?

If two genes are involved in the same pathways usually they have (to
some extent, maybe indirect) interactions. To detect which relationship
is more important between two genes one could measure other similarities
scores and check the stoichiometry of the pathways and measure the
expression changes and correlation between them or use dynamic
simulations of the pathways.

### How to detect with which genes is my gene of interest related?

You can measure the [gene similarity](#geneSim) between those genes and
also measure the expression correlation of your gene of interest with
other genes.

### Why isn’t available a method for calculating GO similarities?

This is covered by the
*[GOSemSim](https://bioconductor.org/packages/3.23/GOSemSim)* package,
you can use it to produce a similarity matrix (i.e. use `mgeneSim`). You
can parallelize it with `foreach` package or `BiocParallel` if your list
of genes is big.

### I get an error! How do I solve this?

If the error is like this:

``` r
Error in FUN(X[[i]], ...) : 
  trying to get slot "geneAnno" from an object of a basic class ("list") with no slots
```

And you have loaded the `GOSemSim` library, R is calling the GOSemSim
function of the same name. Use `BioCor::` to call the function from
`BioCor` (f.ex:
[`BioCor::geneSim`](https://biocor.llrs.dev/reference/geneSim.md))

If the error is not previously described in the [support
forum](https:support.bioconductor.org), post a question there.

My apologies if you found a bug or an inconsistency between what
`BioCor` should do and what it actually does. Once you checked that it
is a bug, please let me know at the
*[issues](https://github.com/llrs/BioCor/issues)* page of Github.

## Session Info

``` r

sessionInfo()
## R version 4.6.1 (2026-06-24)
## Platform: x86_64-pc-linux-gnu
## Running under: Ubuntu 24.04.4 LTS
## 
## Matrix products: default
## BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
## LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
## 
## locale:
##  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
##  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
##  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
## [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
## 
## time zone: UTC
## tzcode source: system (glibc)
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] reactome.db_1.96.0   org.Hs.eg.db_3.23.1  AnnotationDbi_1.74.0
##  [4] IRanges_2.46.0       S4Vectors_0.50.2     Biobase_2.72.0      
##  [7] BiocGenerics_0.58.1  generics_0.1.4       BioCor_1.37.1       
## [10] BiocStyle_2.40.0    
## 
## loaded via a namespace (and not attached):
##  [1] KEGGREST_1.52.2     GOSemSim_2.38.3     xfun_0.60          
##  [4] bslib_0.12.0        htmlwidgets_1.6.4   lattice_0.22-9     
##  [7] vctrs_0.7.3         tools_4.6.1         yulab.utils_0.2.4  
## [10] parallel_4.6.1      RSQLite_3.53.3      blob_1.3.0         
## [13] pkgconfig_2.0.3     Matrix_1.7-5        desc_1.4.3         
## [16] graph_1.90.0        lifecycle_1.0.5     compiler_4.6.1     
## [19] textshaping_1.0.5   Biostrings_2.80.1   Seqinfo_1.2.0      
## [22] codetools_0.2-20    htmltools_0.5.9     sass_0.4.10        
## [25] yaml_2.3.12         pkgdown_2.2.1       crayon_1.5.3       
## [28] jquerylib_0.1.4     GO.db_3.23.1        BiocParallel_1.46.0
## [31] cachem_1.1.0        digest_0.6.39       bookdown_0.48      
## [34] fastmap_1.2.0       grid_4.6.1          cli_3.6.6          
## [37] XML_3.99-0.24       GSEABase_1.74.0     rappdirs_0.3.4     
## [40] bit64_4.8.4         rmarkdown_2.31      XVector_0.52.0     
## [43] httr_1.4.8          bit_4.6.0           otel_0.2.0         
## [46] ragg_1.5.2          png_0.1-9           memoise_2.0.1      
## [49] evaluate_1.0.5      knitr_1.51          rlang_1.3.0        
## [52] xtable_1.8-8        DBI_1.3.0           BiocManager_1.30.27
## [55] annotate_1.90.0     jsonlite_2.0.0      R6_2.6.1           
## [58] systemfonts_1.3.2   fs_2.1.0
```

[^1]: See the [Combining values section](#combining) and the help page
    of `combineScores`.
