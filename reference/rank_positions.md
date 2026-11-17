# Rank positions and tied ranks, within each split

Ties of \`rank_value\` are broken by \`value\`, then by the feature name
(radix sort, independent of the locale), so that the result does not
depend on the row order.

## Usage

``` r
rank_positions(df, decreasing)
```

## Arguments

- df:

  Table from
  [`rank_input()`](https://j-y26.github.io/bioToolkit/reference/rank_input.md)

- decreasing:

  Whether rank 1 is the largest value

## Value

\`df\` ordered by split and rank, with \`rank\` (position, 1 to n) and
\`rank_tied\` (competition rank: tied values share the smallest rank of
the run)
