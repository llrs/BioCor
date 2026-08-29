# BioCor: A package to calculate functional similarities

Calculates a functional similarity measure between gene identifiers
based on the pathways described on KEGG and REACTOME.

## Important functions

- **[`pathSim()`](https://biocor.llrs.dev/reference/pathSim.md)**:
  Calculates the similarity between two pathways.

- **[`geneSim()`](https://biocor.llrs.dev/reference/geneSim.md)**:
  Calculates the similarity (based on pathSim) between two genes.

- **[`clusterSim()`](https://biocor.llrs.dev/reference/clusterSim.md)**:
  Calculates the similarity between two clusters of genes by joining
  pathways of each gene.

- **[`clusterGeneSim()`](https://biocor.llrs.dev/reference/clusterGeneSim.md)**:
  Calculates the similarity between two clusters of genes by comparing
  the similarity between the genes of a cluster.

- **[`similarities()`](https://biocor.llrs.dev/reference/similarities.md)**:
  Allows to combine the value of matrices of similarities.

- **[`conversions()`](https://biocor.llrs.dev/reference/conversions.md)**:
  Two functions to convert similarity measures.

- **[`weighted()`](https://biocor.llrs.dev/reference/weighted.md)**:
  Functions provided to combine similarities.

## See also

Useful links:

- <https://bioconductor.org/packages/BioCor>

- <https://biocor.llrs.dev>

- Report bugs at <https://github.com/llrs/BioCor/issues>

## Author

**Maintainer**: Lluís Revilla Sancho <lluis.revilla@gmail.com>
([ORCID](https://orcid.org/0000-0001-9747-2570))

Authors:

- Lluís Revilla Sancho <lluis.revilla@gmail.com>
  ([ORCID](https://orcid.org/0000-0001-9747-2570))

Other contributors:

- Pau Sancho-Bru ([ORCID](https://orcid.org/0000-0001-5569-9259))
  \[thesis advisor\]

- Juan José Salvatella Lozano
  ([ORCID](https://orcid.org/0000-0001-7613-3908)) \[thesis advisor\]
