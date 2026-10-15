# Plotting units of a dot plot: groups, or split and group combinations

Plotting units of a dot plot: groups, or split and group combinations

## Usage

``` r
dot_units(group, split)
```

## Arguments

- group:

  A factor along the cells (the lower level, compared)

- split:

  A factor along the cells (the higher level, blocks), or NULL

## Value

A list with \`unit\`, a factor along the cells, and \`table\`, one row
per unit with cells (\`unit\`, \`group\` and, when split, \`split\`), in
split then group level order
