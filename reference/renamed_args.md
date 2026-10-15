# Accept deprecated argument names

Older versions of a few functions
([`format_sc_plots()`](https://j-y26.github.io/bioToolkit/reference/format_sc_plots.md),
[`manipulate_sc_plot()`](https://j-y26.github.io/bioToolkit/reference/manipulate_sc_plot.md),
[`add_dimension_arrows()`](https://j-y26.github.io/bioToolkit/reference/add_dimension_arrows.md)
and the DE table helpers) used snake_case argument names. Since 0.10.0
they use the package's dot.case names (see CONTRIBUTING.md). The old
names are still accepted: the caller collects them in \`...\` and passes
them here, which warns once (class "deprecatedWarning") and returns the
values under their new names, ready for `list2env(..., environment())`
in the caller.

## Usage

``` r
renamed_args(dots, renames, fn, given = character())
```

## Arguments

- dots:

  \`list(...)\` of the caller

- renames:

  Character vector of new argument names, named by the old names

- fn:

  Name of the calling function, for messages

- given:

  Names of the arguments the caller received (from
  `names(match.call())`), to refuse an old and a new name given together

## Value

A named list of the values of the old arguments, named by the new names
(empty when no old name was used)
