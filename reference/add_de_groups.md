# Add a group to each comparison

Add a group to each comparison

## Usage

``` r
add_de_groups(res, groups)
```

## Arguments

- res:

  Standardized DE results with a comparison column

- groups:

  A data frame with comparison and group columns, or a character vector
  of groups named by comparison

## Value

\`res\` with a \`group\` factor, its levels in the order the groups
first appear along the comparisons
