# Labels with gene names in italics

Gene symbols are set in italics by convention. This turns labels into
plotmath expressions, italic where they are genes and plain otherwise,
and can be passed directly as the \`labels\` of a discrete scale. Labels
are used as they are: the case of gene symbols (e.g. "CD3E" in human,
"Cd3e" in mouse and rat) is not changed.

## Usage

``` r
italic_labels(x, is.gene = TRUE)
```

## Arguments

- x:

  Character vector of labels

- is.gene:

  Whether each label is a gene: a single TRUE / FALSE for all, a logical
  vector along \`x\`, or a logical vector named by label. Labels missing
  from a named vector are plain

## Value

An expression vector, one element per label

## Examples

``` r
library(ggplot2)
df <- data.frame(gene = c("Cd3e", "Ms4a1", "nCount_RNA"), value = 1:3)
ggplot(df, aes(value, gene)) +
    geom_point() +
    scale_y_discrete(labels = function(x) italic_labels(x, x != "nCount_RNA"))

```
