# Advanced usage of BioCor

Abstract

Describes how to use the BioCor package to answer several biological
questions and how to use functional similarities with other measures.

## Introduction

In this vignette we assume that the reader is already familiarized with
the [introduction
vignette](https://bioconductor.org/packages/3.23/BioCor/vignettes/BioCor_1_basics.html),
but wants to know how can it help to answer other questions in other
situations

We follow the same convention about names of the objects used
`genesReact` and `genesKegg` are the list with information about the
pathways the genes are involved in.

## Merging similarities

If one calculates similarities with KEGG data and Reactome or other
input for the same genes or clusters BioCor provides a couple of
functions to merge them.

We can set a weight to each similarity input with `weighted.sum`,
multiply them also using a weight for each similarity (with
`weighted.prod`), doing the mean or just adding them up. Similarities
allow us to apply a function to combine the matrices of a list. Here we
use some of the genes used in the first vignette:

``` r

kegg <- mgeneSim(c("672", "675", "10"), genesKegg)
react <- mgeneSim(c("672", "675", "10"), genesReact)
## We can sum it adding a weight to each origin
weighted.sum(c(kegg["672", "675"], react["672", "675"]), w = c(0.3, 0.7))
## [1] 0.7247289

## Or if we want to perform for all the matrix
## A list of matrices to merge
sim <- list("kegg" = kegg, "react" = react)
similarities(sim, weighted.sum, w = c(0.3, 0.7))
##           672       675        10
## 672 1.0000000 0.7247289 0.1399444
## 675 0.7247289 1.0000000 0.1424169
## 10  0.1399444 0.1424169 1.0000000
similarities(sim, weighted.prod, w = c(0.3, 0.7))
##           672          675           10
## 672 0.2100000 0.0173101952 0.0000000000
## 675 0.0173102 0.2100000000 0.0003460163
## 10  0.0000000 0.0003460163 0.2100000000
similarities(sim, prod)
##           672         675          10
## 672 1.0000000 0.082429501 0.000000000
## 675 0.0824295 1.000000000 0.001647697
## 10  0.0000000 0.001647697 1.000000000
similarities(sim, mean)
##            672       675         10
## 672 1.00000000 0.5412148 0.09996025
## 675 0.54121475 1.0000000 0.10408113
## 10  0.09996025 0.1040811 1.00000000
```

This functions are similar to `weighted.mean`, except that first the
multiplication by the weights is done and then the `NA`s are removed:

``` r

weighted.mean(c(1, NA), w = c(0.5, 0.5), na.rm = TRUE)
## [1] 1
weighted.mean(c(1, 0.5, NA), w = c(0.5, 0.25, 0.25), na.rm = TRUE)
## [1] 0.8333333
weighted.sum(c(1, NA), w = c(0.5, 0.5))
## [1] 0.5
weighted.sum(c(1, 0.5, NA), w = c(0.5, 0.25, 0.25))
## [1] 0.625
weighted.prod(c(1, NA), w = c(0.5, 0.5))
## [1] 0.5
weighted.prod(c(1, 0.5, NA), w = c(0.5, 0.25, 0.25))
## [1] 0.0625
```

## Assessing a differential study

In this example we will use the functional similarities in a classical
differential study.

### Obtaining data

We start using data from the [RNAseq
workflow](https://bioconductor.org/help/workflows/rnaseqGene/#differential-expression-analysis)
and following the analysis comparing treated and untreated:

``` r

suppressPackageStartupMessages(library("airway"))
data("airway")
library("DESeq2")

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds$dex <- relevel(dds$dex, "untrt")
dds <- DESeq(dds)
## estimating size factors
## estimating dispersions
## gene-wise dispersion estimates
## mean-dispersion relationship
## final dispersion estimates
## fitting model and testing
res <- results(dds, alpha = 0.05)
summary(res)
## 
## out of 33469 with nonzero total read count
## adjusted p-value < 0.05
## LFC > 0 (up)       : 2211, 6.6%
## LFC < 0 (down)     : 1817, 5.4%
## outliers [1]       : 0, 0%
## low counts [2]     : 16687, 50%
## (mean count < 7)
## [1] see 'cooksCutoff' argument of ?results
## [2] see 'independentFiltering' argument of ?results
plot(res$log2FoldChange, -log10(res$padj),
    pch = 16, xlab = "log2FC",
    ylab = "-log10(p.ajd)", main = "Untreated vs treated"
)
logFC <- 2.5
abline(v = c(-logFC, logFC), h = -log10(0.05), col = "red")
```

![Volcano plot (log2FC on the horitzonal axis and log(p-value) on the
vertical axis) of the airway
dataset.](BioCor_2_advanced_files/figure-html/simulate-1.png)

Volcano plot. The airway data

As we can see here there are around 4000 genes differentially expressed
genes, some of them differentially expressed above 2^2.5.

### Selecting differentially expressed genes

Usually in such a study one selects genes above certain logFC or fold
change threshold, here we use the absolute value of 2.5:

``` r

fc <- res[abs(res$log2FoldChange) >= logFC & !is.na(res$padj), ]
fc <- fc[fc$padj < 0.05, ]
# Convert Ids (used later)
genes <- select(org.Hs.eg.db,
    keys = rownames(res), keytype = "ENSEMBL",
    column = c("ENTREZID", "SYMBOL")
)
## 'select()' returned 1:many mapping between keys and columns
genesFC <- genes[genes$ENSEMBL %in% rownames(fc), ]
genesFC <- genesFC[!is.na(genesFC$ENTREZID), ]
genesSim <- genesFC$ENTREZID
names(genesSim) <- genesFC$SYMBOL
genesSim <- genesSim[!duplicated(genesSim)]
# Calculate the functional similarity
sim <- mgeneSim(genes = genesSim, info = genesReact, method = "BMA")
## Warning in mgeneSim(genes = genesSim, info = genesReact, method = "BMA"): Some
## genes are not in the list provided.
```

Once the similarities for the selected genes are calculated we can now
visualize the effect of each method:

``` r

nas <- apply(sim, 1, function(x) {
    all(is.na(x))
})
sim <- sim[!nas, !nas]

MDSs <- cmdscale(1 - sim)
plot(MDSs, type = "n", main = "BMA similarity", xlab = "MDS1", ylab = "MDS2")
up <- genes[genes$ENSEMBL %in% rownames(fc)[fc$log2FoldChange >= logFC], "SYMBOL"]
text(MDSs, labels = rownames(MDSs), col = ifelse(rownames(MDSs) %in% up, "black", "red"))
abline(h = 0, v = 0)
legend("top", legend = c("Up-regulated", "Down-regulated"), fill = c("black", "red"))
```

![First two dimensions of a multi dimensional scaling method based on
the similarity of genes. Colored if the fold change is above 2.5
(red).](BioCor_2_advanced_files/figure-html/pval1-1.png)

Functional similarity of genes with logFC above 2,5. Similar genes
cluster together.

This plot illustrate that some differentially expressed genes are quite
similar according to their pathways. Suggesting that there might be a
relationship between them. Furthermore, some up-regulated genes seem
functionally related to down-regulated genes indicating a precise
regulation of the pathways where they are involved.

Note that here we are only using 74 genes from the original 127.

### Are differentially expressed genes selected by their functionality?

In the previous section we have seen that some differentially expressed
genes are functionally related and that they have a high logFC value.
Are genes differentially expressed more functional related than those
which aren’t differential expressed?

For simplicity we will use a subset of 400 genes represented again in a
volcano plot and we will look for the functional similarities between
those genes:

``` r

set.seed(250)
subRes <- res[!is.na(res$log2FoldChange), ]
subs <- sample.int(nrow(subRes), size = 400)
subRes <- subRes[subs, ]
g <- genes[genes$ENSEMBL %in% rownames(subRes), "ENTREZID"]
gS <- mgeneSim(g[g %in% names(genesReact)], genesReact, "BMA")
deg <- rownames(subRes[subRes$padj < 0.05 & !is.na(subRes$padj), ])
keep <- rownames(gS) %in% genes[genes$ENSEMBL %in% deg, "ENTREZID"]
plot(subRes$log2FoldChange, -log10(subRes$padj),
    pch = 16, xlab = "log2FC",
    ylab = "-log10(p.ajd)", main = "Untreated vs treated"
)
abline(v = c(-logFC, logFC), h = -log10(0.05), col = "red")
```

![Volcano plot of a subset of 400 genes.
](BioCor_2_advanced_files/figure-html/setting-1.png)

Volcano plot of the subset of 400 genes. This subset will be used in the
following sections

We can answer this by testing it empirically:

``` r

library("boot")
# The mean of genes differentially expressed
(scoreDEG <- mean(gS[!keep, keep], na.rm = TRUE))
## [1] 0.1228395
b <- boot(data = gS, R = 1000, statistic = function(x, i) {
    g <- !rownames(x) %in% rownames(x)[i]
    mean(x[g, i], na.rm = TRUE)
})
(p.val <- (1 + sum(b$t > scoreDEG)) / 1001)
## [1] 0.9170829
hist(b$t, main = "Distribution of scores", xlab = "Similarity score")
abline(v = scoreDEG, col = "red")
```

![Histogram of the scores of similarity between several genes. A
vertical red line indicates the score of between those differentially
expressed and those which
aren't.](BioCor_2_advanced_files/figure-html/cluster2-1.png)

Distribution of scores between differentially expressed genes and those
who aren’t. The line indicates the mean score of the similarity between
differentially expressed genes and those which aren’t differentially
expressed.

Comparing the genes differentially expressed and those who aren’t
doesn’t show that they are non-randomly selected (p-value 0.9170829).
The genes with a p-value below the threshold are not more closely
functionally related than all the other genes[^1].

### Are functionally related the selected differentially expressed genes?

We have seen that the genes differentially expressed aren’t selected by
their functionality. However they could be more functionally related
than the other genes. Are the differentially expressed genes more
functionally similar than it would be expected ?

``` r

(scoreW <- combineScores(gS[keep, keep], "avg"))
## [1] 0.1580637
b <- boot(data = gS, R = 1000, statistic = function(x, i) {
    mean(x[i, i], na.rm = TRUE)
})
(p.val <- (1 + sum(b$t > scoreW)) / 1001) # P-value
## [1] 0.01298701
hist(b$t, main = "Distribution of scores", xlab = "Similarity score")
abline(v = scoreW, col = "red")
```

![Distribution of similarity scores between expressed genes. A vertical
red line indicates the real similarity between those differentially
expressed genes.](BioCor_2_advanced_files/figure-html/pval2-1.png)

Distribution of the similarity within differentially expressed genes.
The line indicates the mean funtional similarity whitin them.

If we selected randomly the genes from our pool we would expect a score
around 0.1580637 with a probability of 0.012987. That means that the
differentially expressed genes is highly different compared to the other
genes if we use a significance threshold of 0.05 [^2].

### Influence of the fold change in the functionally similarity of the genes

We have seen that the genes differentially expressed are not selected by
functional similarity but they are functionally related. Now we would
like to know if selecting a fold change threshold affects the functional
similarity between them.

To know the relationship between the fold change and the similarity
between genes we have several methods:

``` r

s <- seq(0, max(abs(subRes$log2FoldChange)) + 0.05, by = 0.05)
l <- sapply(s, function(x) {
    deg <- rownames(subRes[abs(subRes$log2FoldChange) >= x, ])
    keep <- rownames(gS) %in% genes[genes$ENSEMBL %in% deg, "ENTREZID"]
    BetweenAbove <- mean(gS[keep, keep], na.rm = TRUE)
    AboveBelow <- mean(gS[keep, !keep], na.rm = TRUE)
    BetweenBelow <- mean(gS[!keep, !keep], na.rm = TRUE)
    c(
        "BetweenAbove" = BetweenAbove, "AboveBelow" = AboveBelow,
        "BetweenBelow" = BetweenBelow
    )
})
L <- as.data.frame(cbind(logfc = s, t(l)))
plot(L$logfc, L$BetweenAbove,
    type = "l", xlab = "abs(log2) fold change",
    ylab = "Similarity score",
    main = "Similarity scores along logFC", col = "darkred"
)
lines(L$logfc, L$AboveBelow, col = "darkgreen")
lines(L$logfc, L$BetweenBelow, col = "black")
legend("topleft",
    legend = c(
        "Between genes above and below threshold",
        "Whitin genes above threshold",
        "Whitin genes below threshold"
    ),
    fill = c("darkgreen", "darkred", "black")
)
```

![A line plot on the X axis the absolute log fold change of the
threshold used, on the vertical axis the similarity score. Three lines,
in red the similarity within genes above the threshold, in black those
below the threshold, in green between above the threshold and below the
threshold.](BioCor_2_advanced_files/figure-html/logfc1-1.png)

Similarity of genes along abs(logFC). Assessing the similarity of genes
according to their absolute log2 fold change.

The functional similarity of the genes above the threshold increases
with a more restrictive threshold, indicating that a logFC threshold
selects genes by their functionality, or in other words that genes
differentially expressed tend to be of related pathways. The similarity
between those genes above the threshold and below remains constant as
well as within genes below the threshold.

``` r

l <- sapply(s, function(x) {
    # Names of genes up and down regulated
    degUp <- rownames(subRes[subRes$log2FoldChange >= x, ])
    degDown <- rownames(subRes[subRes$log2FoldChange <= -x, ])

    # Translate to ids in gS
    keepUp <- rownames(gS) %in% genes[genes$ENSEMBL %in% degUp, "ENTREZID"]
    keepDown <- rownames(gS) %in% genes[genes$ENSEMBL %in% degDown, "ENTREZID"]

    # Calculate the mean similarity between each subgrup
    between <- mean(gS[keepUp, keepDown], na.rm = TRUE)

    c("UpVsDown" = between)
})
L <- as.data.frame(cbind("logfc" = s, "UpVsDown" = l))
plot(L$logfc, L$UpVsDown,
    type = "l",
    xlab = "abs(log2) fold change threshold",
    ylab = "Similarity score",
    main = "Similarity scores along logFC"
)
legend("topright",
    legend = "Up vs down regulated genes",
    fill = "black"
)
```

![Similarity score between genes up and down-regulated along the
threshold of log fold
change.](BioCor_2_advanced_files/figure-html/logfc2-1.png)

Functional similarity between the up-regulated and down-regulated genes.

The maximal functional similarity between genes up-regulated and
down-regulated are at 1.15 log2 fold change.

## Assessing a new pathway

Sometimes the top differentially expressed genes or some other key genes
are selected as a signature or a potential new group of related genes.
In those cases we can test how does the network of genes change if we
add them. Here we create a new pathway named `deg`, and we see the
effect on the functional similarity score for all the genes:

``` r

# Adding a new pathway "deg" to those genes
genesReact2 <- genesReact
diffGenes <- genes[genes$ENSEMBL %in% deg, "ENTREZID"]
genesReact2[diffGenes] <- sapply(genesReact[diffGenes], function(x) {
    c(x, "deg")
})
plot(ecdf(mgeneSim(names(genesReact), genesReact)))
curve(ecdf(mgeneSim(names(genesReact2), genesReact2)), color = "red")
```

This would take lot of time, for a illustration purpose we reduce to
some genes:

``` r

library("Hmisc")
genesReact2 <- genesReact
diffGenes <- genes[genes$ENSEMBL %in% deg, "ENTREZID"]
# Create the new pathway called deg
genesReact2[diffGenes] <- sapply(genesReact[diffGenes], function(x) {
    c(x, "deg")
})
ids <- unique(genes[genes$ENSEMBL %in% rownames(subRes), "ENTREZID"])
Ecdf(
    c(
        mgeneSim(ids, genesReact, method = "BMA"),
        mgeneSim(ids, genesReact2, method = "BMA")
    ),
    group = c(rep("Reactome", length(ids)^2), rep("Modified", length(ids)^2)),
    col = c("black", "red"), xlab = "Functional similarities",
    main = "Empirical cumulative distribution"
)
```

![Empirical cumulative distribution of the functional similarity with
the original data (colored in red) and with an added pathway (in
black).](BioCor_2_advanced_files/figure-html/newPathway2-1.png)

The effect of adding a new pathway to a functional similarity. In red
the same network as in black but with the added pathway.

## Merging sources of information

Sometimes we have several origin of information, either several
databases, or information from other programs… We can merge this in the
single object required by the function in `BioCor` using the function
`combineSources`[^3].

This functions helps to evaluate what happens when we add more pathway
information. For instance here we add the information in Kegg and the
information in Reactome and we visualize it using the same procedure as
previously:

``` r

genesKegg <- as.list(org.Hs.egPATH)
gSK <- mgeneSim(rownames(gS), genesKegg)
mix <- combineSources(genesKegg, genesReact)
## Warning in combineSources(genesKegg, genesReact): More than 50% of genes identifiers of a source are unique
## Check the identifiers of the genes
gSMix <- mgeneSim(rownames(gS), mix)
Ecdf(c(gS, gSK, gSMix),
    group = c(
        rep("Reactome", length(gS)), rep("Kegg", length(gSK)),
        rep("Mix", length(gSMix))
    ),
    col = c("black", "red", "blue"), xlab = "Functional similarities",
    main = "Empirical cumulative distribution."
)
## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
## collapsing to unique 'x' values
## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
## collapsing to unique 'x' values
## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
## collapsing to unique 'x' values
```

![Comparing the functional similarity by looking at the empirical
cumulative distribution. Kegg in black, Reactome in blue, and both mixed
in red.](BioCor_2_advanced_files/figure-html/combineSource-1.png)

Comparison of functional similarity in different databases.

When mixed, there is a huge increase in the genes that share a pathway
using the `max` method. Observe in next figure
(@ref(fig:combineSource2)) how does the method affects to the results:

``` r

gSK2 <- mgeneSim(rownames(gS), genesKegg, method = "BMA")
gS2 <- mgeneSim(rownames(gS), genesReact, method = "BMA")
gSMix2 <- mgeneSim(rownames(gS), mix, method = "BMA")
Ecdf(c(gS2, gSK2, gSMix2),
    group = c(
        rep("Reactome", length(gS)), rep("Kegg", length(gSK)),
        rep("Mix", length(gSMix))
    ),
    col = c("black", "red", "blue"), xlab = "Functional similarities (BMA)", main = "Empirical cumulative distribution."
)
## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
## collapsing to unique 'x' values
## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
## collapsing to unique 'x' values
## Warning in regularize.values(x, y, ties, missing(ties), na.rm = na.rm):
## collapsing to unique 'x' values
```

![Empirical cumulative distribution of pathways according to Kegg
(black), Reactome (blue) and a mix of both
(red).](BioCor_2_advanced_files/figure-html/combineSource2-1.png)

Comparison of functional similarity in different gene sets.

Now we can appreciate that most of the functional similarity is brought
by Reactome database

## Comparing with GO similarities

As suggested in the main vignette functional similarities can be
compared to semantic similarities such as those based on GO. Here a
comparison using the biological process from the gene ontologies is
done:

``` r

library("GOSemSim")
## 
## GOSemSim v2.38.3 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
## 
## Please cite:
## 
## Guangchuang Yu. Gene Ontology Semantic Similarity Analysis Using
## GOSemSim. In: Kidder B. (eds) Stem Cell Transcriptional Networks.
## Methods in Molecular Biology. 2020, 2117:207-215. Humana, New York, NY.
## 
## Attaching package: 'GOSemSim'
## The following objects are masked from 'package:BioCor':
## 
##     clusterSim, combineScores, geneSim, mclusterSim, mgeneSim
BP <- godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE)
## Warning in godata("org.Hs.eg.db", ont = "BP", computeIC = TRUE): use 'annoDb'
## instead of 'OrgDb'
## preparing gene to GO mapping data...
## preparing IC data...
gsGO <- GOSemSim::mgeneSim(rownames(gS), semData = BP, measure = "Resnik", verbose = FALSE)
keep <- rownames(gS) %in% rownames(gsGO)
hist(as.dist(gS[keep, keep] - gsGO),
    main = "Difference between functional similarity and biological process",
    xlab = "Functional similarity - biological process similarity"
)
```

![Functional similarities for the same genes compared between GO and
Reactome
annotation.](BioCor_2_advanced_files/figure-html/GOSemSim-1.png)

Comparison of similarities. Functional similarities compared to
biological process semantic similarity.

On this graphic we can observe that some genes have a large functional
similarity and few biological similarity. They are present together in
several pathways while they share few biological process, indicating
that they might be key elements of the pathways they are in. On the
other hand some other pairs of genes show higher biological process
similarity than functional similarity, indicating specialization or
compartmentalization of said genes.

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
##  [1] GOSemSim_2.38.3             Hmisc_5.2-6                
##  [3] boot_1.3-32                 DESeq2_1.52.0              
##  [5] airway_1.32.0               SummarizedExperiment_1.42.0
##  [7] GenomicRanges_1.64.0        Seqinfo_1.2.0              
##  [9] MatrixGenerics_1.24.0       matrixStats_1.5.0          
## [11] BioCor_1.37.1               org.Hs.eg.db_3.23.1        
## [13] AnnotationDbi_1.74.0        IRanges_2.46.0             
## [15] S4Vectors_0.50.2            Biobase_2.72.0             
## [17] BiocGenerics_0.58.1         generics_0.1.4             
## [19] BiocStyle_2.40.0           
## 
## loaded via a namespace (and not attached):
##  [1] DBI_1.3.0           gridExtra_2.3.1     GSEABase_1.74.0    
##  [4] rlang_1.3.0         magrittr_2.0.5      otel_0.2.0         
##  [7] compiler_4.6.1      RSQLite_3.53.3      png_0.1-9          
## [10] systemfonts_1.3.2   vctrs_0.7.3         reactome.db_1.96.0 
## [13] stringr_1.6.0       pkgconfig_2.0.3     crayon_1.5.3       
## [16] fastmap_1.2.0       backports_1.5.1     XVector_0.52.0     
## [19] rmarkdown_2.31      graph_1.90.0        ragg_1.5.2         
## [22] bit_4.6.0           xfun_0.60           cachem_1.1.0       
## [25] jsonlite_2.0.0      blob_1.3.0          DelayedArray_0.38.2
## [28] BiocParallel_1.46.0 parallel_4.6.1      cluster_2.1.8.2    
## [31] R6_2.6.1            bslib_0.12.0        stringi_1.8.9      
## [34] RColorBrewer_1.1-3  rpart_4.1.27        jquerylib_0.1.4    
## [37] Rcpp_1.1.2          bookdown_0.48       knitr_1.51         
## [40] base64enc_0.1-6     Matrix_1.7-5        nnet_7.3-20        
## [43] rstudioapi_0.19.0   abind_1.4-8         yaml_2.3.12        
## [46] codetools_0.2-20    lattice_0.22-9      KEGGREST_1.52.2    
## [49] S7_0.2.2            evaluate_1.0.5      foreign_0.8-91     
## [52] desc_1.4.3          Biostrings_2.80.1   BiocManager_1.30.27
## [55] checkmate_2.3.4     ggplot2_4.0.3       scales_1.4.0       
## [58] xtable_1.8-8        glue_1.8.1          tools_4.6.1        
## [61] data.table_1.18.6.1 annotate_1.90.0     locfit_1.5-9.12    
## [64] fs_2.1.0            XML_3.99-0.24       grid_4.6.1         
## [67] colorspace_2.1-3    htmlTable_2.5.0     Formula_1.2-6      
## [70] cli_3.6.6           rappdirs_0.3.4      textshaping_1.0.5  
## [73] S4Arrays_1.12.0     gtable_0.3.6        yulab.utils_0.2.5  
## [76] sass_0.4.10         digest_0.6.39       SparseArray_1.12.2 
## [79] htmlwidgets_1.6.4   farver_2.1.2        memoise_2.0.1      
## [82] htmltools_0.5.9     pkgdown_2.2.1       lifecycle_1.0.5    
## [85] httr_1.4.8          GO.db_3.23.1        bit64_4.8.4
```

[^1]: From 400 genes there are 127 with pathway information and only 25
    where significantly differentially expressed in this subset.

[^2]: Remember that this is a small subset.

[^3]: See the help page of `combineSources`
