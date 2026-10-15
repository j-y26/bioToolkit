# Running enrichment score of a gene set along a ranked list

Running enrichment score of a gene set along a ranked list

## Usage

``` r
gsea_running_score(stats, set, p = 1)
```

## Arguments

- stats:

  Ranking metric, named by gene

- set:

  Genes of the set

- p:

  Weight exponent, default is 1 (GSEA, fgsea, clusterProfiler)

## Value

A data frame with one row per gene, largest metric first: \`rank\`,
\`gene\`, \`metric\`, \`running_es\` and \`in_set\`
