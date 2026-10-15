# Colours for the levels of a grouping variable

Every level gets a colour, including levels with no cells, so that a
group gets the same colour in a plot of a subset of cells as in a plot
of all of them.

## Usage

``` r
group_colours(cols, group, arg = "cols")
```

## Arguments

- cols:

  NULL for the default colours; a qualitative palette name; colours
  named by group; or colours in the order of the group levels

- group:

  A factor

- arg:

  Name of the argument, for messages

## Value

A list with \`group\` and \`cols\`, named by group
