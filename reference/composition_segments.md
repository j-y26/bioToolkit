# Segments of stacked composition bars

One row per segment (a group in a bar), placed and stacked. Bars are the
conditions of
[`summarize_group_meta()`](https://j-y26.github.io/bioToolkit/reference/summarize_group_meta.md);
each panel (\`facet\`) holds only the bars with cells in it. Positions
run on over the whole plot (the first panel's bars are 1..n, the next
panel's n + 1...), so that a bar name repeated in two panels gets its
own place in each.

## Usage

``` r
composition_segments(
  tab,
  value = c("percent", "fraction", "count"),
  order.by = NULL,
  first.at.start = FALSE
)
```

## Arguments

- tab:

  Result of
  [`summarize_group_meta()`](https://j-y26.github.io/bioToolkit/reference/summarize_group_meta.md)

- value:

  "percent" or "fraction" (the proportion) or "count"

- order.by:

  A group: the bars of each panel are sorted by decreasing proportion of
  that group (ties keep their order). NULL keeps the order of the levels

- first.at.start:

  Whether the first group is stacked at the start of the value axis
  (left of a horizontal bar). FALSE puts it at the end (top of a
  vertical bar), so that the stack reads in the order of the legend

## Value

A data frame: \`facet\` (if any), \`bar\`, \`pos\`, \`group\`,
\`count\`, \`total\`, \`prop\`, \`value\` (drawn), \`lower\` and
\`upper\` (ends of the segment on the value axis)
