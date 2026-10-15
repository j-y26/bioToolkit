# Give the nested patchworks of a patchwork a fixed size

patchwork aligns the panels of the plots it combines, but a nested
patchwork (e.g. \`a \| b\` in \`(a \| b) / c\`) has no panel of its own,
so its fixed size is lost: the outer layout is left with relative sizes,
and plots sharing a row or column with it are stretched. On every level
that holds a nested patchwork, each plot is therefore drawn into a
ggplot whose single panel has exactly that plot's size, centred in its
cell. Levels without nesting are left to patchwork's own alignment. Must
be called with the device the plot is measured on open.

## Usage

``` r
fix_nested_patchworks(plot)
```

## Arguments

- plot:

  A patchwork object

## Value

The patchwork, rebuilt if it holds nested patchworks
