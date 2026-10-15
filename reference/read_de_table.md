# Read one DE table into standard columns

Read one DE table into standard columns

## Usage

``` r
read_de_table(data, method, cols, label = NULL, comparison.col = NULL)
```

## Arguments

- data:

  A DE table

- method:

  Method, or "auto"

- cols:

  Columns given explicitly, as a list named by standard column

- label:

  Name of the comparison, for messages, or NULL

- comparison.col:

  Column naming the comparison of each row, or NULL

## Value

A list with the standardized \`table\`, the \`method\` used and the
\`mean_scale\`
