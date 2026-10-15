# A legend title in bold, also when it is a plotmath expression

The \`face\` of a theme element does not apply to plotmath expressions
(e.g. \`quote(-log\[10\] ~ "padj")\`), which are drawn plain unless
wrapped in \`bold()\`. Strings are returned as they are (the theme makes
them bold).

## Usage

``` r
bold_title(x)
```

## Arguments

- x:

  A string, a plotmath expression, or NULL

## Value

\`x\`, with expressions wrapped in \`bold()\` (italic ones in
\`bolditalic()\`)
