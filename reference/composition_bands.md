# Bands joining the segments of neighbouring bars

For each pair of neighbouring bars of a panel and each group, a polygon
from the group's segment on the first bar (its edge facing the second)
to its segment on the second, in a smooth S ("curve", as in alluvial
plots) or a straight line. Groups empty in both bars get no band; bars
of different panels are never joined.

## Usage

``` r
composition_bands(seg, half, shape = c("curve", "straight"), n.points = 40)
```

## Arguments

- seg:

  Result of
  [`composition_segments()`](https://j-y26.github.io/bioToolkit/reference/composition_segments.md)

- half:

  Half the width of a bar, in units of the category axis

- shape:

  "curve" or "straight"

- n.points:

  Points along each curved edge

## Value

A data frame: \`facet\` (if any), \`group\`, \`id\` (one per polygon),
\`pos\` and \`val\` (the polygon's corners). No rows if there is nothing
to join
