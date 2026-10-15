# Mean expression and percentage of expressing cells per group

Mean expression and percentage of expressing cells per group

## Usage

``` r
group_expression_stats(expr, unit, average, threshold, log.data)
```

## Arguments

- expr:

  A matrix of cells x genes

- unit:

  A factor along the cells

- average:

  "linear" or "log"

- threshold:

  A cell expresses a gene above this value

- log.data:

  Whether \`expr\` holds log1p values (the "data" layer), for which
  \`average\` applies

## Value

A list of three matrices, units x genes: \`mean\`, \`pct\` and \`n\`
(number of cells)
