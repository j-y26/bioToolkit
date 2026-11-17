# Palettes of the panels of an enrichment dot plot

Palettes of the panels of an enrichment dot plot

## Usage

``` r
enrichment_panel_palettes(cols, panels, single = FALSE)
```

## Arguments

- cols:

  NULL; one palette name or colour vector for all panels; palette names,
  one per panel; or a list of colour vectors

- panels:

  Panel names

- single:

  Whether the plot has one result; with \`cols = NULL\`, it is drawn in
  reversed magma
  ([`magma_r_trimmed()`](https://j-y26.github.io/bioToolkit/reference/magma_r_trimmed.md)).
  Default is FALSE

## Value

A list of colour vectors, named by panel
