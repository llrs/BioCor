Tests and Coverage
================
15 septiembre, 2018 12:21:00

-   [Coverage](#coverage)
-   [Unit Tests](#unit-tests)

This output is created by [covrpage](https://github.com/yonicd/covrpage).

Coverage
--------

Coverage summary is created using the [covr](https://github.com/r-lib/covr) package.

| Object                                        | Coverage (%) |
|:----------------------------------------------|:------------:|
| BioCor                                        |     87.25    |
| [R/geneSim.R](../R/geneSim.R)                 |     53.19    |
| [R/mgeneSim.R](../R/mgeneSim.R)               |     53.70    |
| [R/pathSim.R](../R/pathSim.R)                 |     69.23    |
| [R/mpathSim.R](../R/mpathSim.R)               |     86.49    |
| [R/combineScores.R](../R/combineScores.R)     |     86.55    |
| [R/mclusterSim.R](../R/mclusterSim.R)         |     91.38    |
| [R/clusterSim.R](../R/clusterSim.R)           |     95.00    |
| [R/clusterGeneSim.R](../R/clusterGeneSim.R)   |     98.55    |
| [R/AintoB.R](../R/AintoB.R)                   |    100.00    |
| [R/auxiliar.R](../R/auxiliar.R)               |    100.00    |
| [R/combineSources.R](../R/combineSources.R)   |    100.00    |
| [R/conversions.R](../R/conversions.R)         |    100.00    |
| [R/diceSim.R](../R/diceSim.R)                 |    100.00    |
| [R/mclusterGeneSim.R](../R/mclusterGeneSim.R) |    100.00    |
| [R/similarities.R](../R/similarities.R)       |    100.00    |
| [R/weighted.R](../R/weighted.R)               |    100.00    |
| [R/zzz.R](../R/zzz.R)                         |    100.00    |

<br>

Unit Tests
----------

Unit Test summary is created using the [testthat](https://github.com/r-lib/testthat) package.

|                        | file                                                     |    n|   time|  error|  failed|  skipped|  warning|
|------------------------|:---------------------------------------------------------|----:|------:|------:|-------:|--------:|--------:|
| test\_AintoB.R         | [test\_AintoB.R](testthat/test_AintoB.R)                 |    6|  0.314|      0|       0|        0|        0|
| test\_auxiliar.R       | [test\_auxiliar.R](testthat/test_auxiliar.R)             |   16|  0.040|      0|       0|        0|        0|
| test\_clusterGeneSim.R | [test\_clusterGeneSim.R](testthat/test_clusterGeneSim.R) |   49|  0.202|      0|       0|        0|        0|
| test\_clusterSim.R     | [test\_clusterSim.R](testthat/test_clusterSim.R)         |   57|  0.200|      0|       0|        0|        0|
| test\_combineScores.R  | [test\_combineScores.R](testthat/test_combineScores.R)   |   62|  0.104|      0|       0|        0|        0|
| test\_combineSources.R | [test\_combineSources.R](testthat/test_combineSources.R) |    5|  0.009|      0|       0|        0|        0|
| test\_conversions.R    | [test\_conversions.R](testthat/test_conversions.R)       |    7|  0.010|      0|       0|        0|        0|
| test\_diceSim.R        | [test\_diceSim.R](testthat/test_diceSim.R)               |    7|  0.010|      0|       0|        0|        0|
| test\_geneSim.R        | [test\_geneSim.R](testthat/test_geneSim.R)               |   23|  0.083|      0|       0|        0|        0|
| test\_mpathSim.R       | [test\_mpathSim.R](testthat/test_mpathSim.R)             |   25|  0.057|      0|       0|        0|        0|
| test\_pathSim.R        | [test\_pathSim.R](testthat/test_pathSim.R)               |    6|  0.011|      0|       0|        0|        0|
| test\_similarities.R   | [test\_similarities.R](testthat/test_similarities.R)     |   18|  0.028|      0|       0|        0|        0|
| test\_weighted.R       | [test\_weighted.R](testthat/test_weighted.R)             |   14|  0.017|      0|       0|        0|        0|

<details closed> <summary> Show Detailed Test Results </summary>

| file                                                             | context                       | test                                              | status |    n|   time|
|:-----------------------------------------------------------------|:------------------------------|:--------------------------------------------------|:-------|----:|------:|
| [test\_AintoB.R](testthat/test_AintoB.R#L10)                     | Testing AintoB functions      | AintoB                                            | PASS   |    6|  0.314|
| [test\_auxiliar.R](testthat/test_auxiliar.R#L8)                  | Testing auxiliar functions    | seq2mat puts a combination to the right place     | PASS   |    5|  0.013|
| [test\_auxiliar.R](testthat/test_auxiliar.R#L24)                 | Testing auxiliar functions    | combinadic                                        | PASS   |    3|  0.006|
| [test\_auxiliar.R](testthat/test_auxiliar.R#L35)                 | Testing auxiliar functions    | duplicateIndices                                  | PASS   |    3|  0.006|
| [test\_auxiliar.R](testthat/test_auxiliar.R#L50)                 | Testing auxiliar functions    | removeDup                                         | PASS   |    4|  0.012|
| [test\_auxiliar.R](testthat/test_auxiliar.R#L62_L63)             | Testing auxiliar functions    | inverseList                                       | PASS   |    1|  0.003|
| [test\_clusterGeneSim.R](testthat/test_clusterGeneSim.R#L6_L7)   | Testing clusterGeneSim        | clusterGeneSim                                    | PASS   |   13|  0.081|
| [test\_clusterGeneSim.R](testthat/test_clusterGeneSim.R#L31)     | Testing clusterGeneSim        | mclusterGeneSim                                   | PASS   |    9|  0.040|
| [test\_clusterGeneSim.R](testthat/test_clusterGeneSim.R#L53_L54) | Testing clusterGeneSim        | clusterGeneSim GeneSetCollection                  | PASS   |   16|  0.045|
| [test\_clusterGeneSim.R](testthat/test_clusterGeneSim.R#L75)     | Testing clusterGeneSim        | mclusterGeneSim GeneSetCollection                 | PASS   |   11|  0.036|
| [test\_clusterSim.R](testthat/test_clusterSim.R#L8)              | Testing clusterSim            | clusterSim                                        | PASS   |   13|  0.056|
| [test\_clusterSim.R](testthat/test_clusterSim.R#L30)             | Testing clusterSim            | mclusterSim                                       | PASS   |   16|  0.069|
| [test\_clusterSim.R](testthat/test_clusterSim.R#L72)             | Testing clusterSim            | clusterSim                                        | PASS   |   14|  0.036|
| [test\_clusterSim.R](testthat/test_clusterSim.R#L92)             | Testing clusterSim            | mclusterSim                                       | PASS   |   14|  0.039|
| [test\_combineScores.R](testthat/test_combineScores.R#L7)        | Testing combineScores         | combineScores                                     | PASS   |   42|  0.063|
| [test\_combineScores.R](testthat/test_combineScores.R#L77)       | Testing combineScores         | combineScoresPar                                  | PASS   |    1|  0.003|
| [test\_combineScores.R](testthat/test_combineScores.R#L84)       | Testing combineScores         | combineScoresPar equivalent to combineScores      | PASS   |   12|  0.025|
| [test\_combineScores.R](testthat/test_combineScores.R#L125)      | Testing combineScores         | reciprocal                                        | PASS   |    5|  0.009|
| [test\_combineScores.R](testthat/test_combineScores.R#L145)      | Testing combineScores         | BMA                                               | PASS   |    1|  0.002|
| [test\_combineScores.R](testthat/test_combineScores.R#L149)      | Testing combineScores         | rcmax                                             | PASS   |    1|  0.002|
| [test\_combineSources.R](testthat/test_combineSources.R#L9)      | Testing combineSources        | combineSources                                    | PASS   |    5|  0.009|
| [test\_conversions.R](testthat/test_conversions.R#L5)            | Testing conversions functions | Conversions                                       | PASS   |    7|  0.010|
| [test\_diceSim.R](testthat/test_diceSim.R#L9)                    | Testing diceSim               | diceSim                                           | PASS   |    7|  0.010|
| [test\_geneSim.R](testthat/test_geneSim.R#L9)                    | Testing geneSim               | geneSim                                           | PASS   |   11|  0.036|
| [test\_geneSim.R](testthat/test_geneSim.R#L26)                   | Testing geneSim               | mgeneSim                                          | PASS   |   12|  0.047|
| [test\_mpathSim.R](testthat/test_mpathSim.R#L6)                  | Testing mpathSim              | mpathSim                                          | PASS   |    7|  0.014|
| [test\_mpathSim.R](testthat/test_mpathSim.R#L23)                 | Testing mpathSim              | pathSims\_matrix                                  | PASS   |    3|  0.005|
| [test\_mpathSim.R](testthat/test_mpathSim.R#L37)                 | Testing mpathSim              | mpathSim for GeneSetCollections and list is equal | PASS   |   10|  0.032|
| [test\_mpathSim.R](testthat/test_mpathSim.R#L77)                 | Testing mpathSim              | mpathSim for GeneSetCollections                   | PASS   |    5|  0.006|
| [test\_pathSim.R](testthat/test_pathSim.R#L6)                    | Testing pathSim               | pathSim                                           | PASS   |    6|  0.011|
| [test\_similarities.R](testthat/test_similarities.R#L12)         | Testing Similarities function | addSimilarities                                   | PASS   |    9|  0.013|
| [test\_similarities.R](testthat/test_similarities.R#L32)         | Testing Similarities function | similarities                                      | PASS   |    9|  0.015|
| [test\_weighted.R](testthat/test_weighted.R#L9)                  | Testing weighted functions    | weighted.sum                                      | PASS   |    8|  0.010|
| [test\_weighted.R](testthat/test_weighted.R#L28)                 | Testing weighted functions    | weighted.prod                                     | PASS   |    6|  0.007|

</details>

<details> <summary> Session Info </summary>

| Field    | Value                         |
|:---------|:------------------------------|
| Version  | R version 3.5.1 (2018-07-02)  |
| Platform | x86\_64-pc-linux-gnu (64-bit) |
| Running  | Ubuntu 18.04.1 LTS            |
| Language | en\_US                        |
| Timezone | Europe/Madrid                 |

| Package  | Version |
|:---------|:--------|
| testthat | 2.0.0   |
| covr     | 3.1.0   |
| covrpage | 0.0.55  |

</details>

<!--- Final Status : pass --->
