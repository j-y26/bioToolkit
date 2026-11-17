# Labels of adjusted p-values

Labels of adjusted p-values

## Usage

``` r
format_padj(x)
```

## Arguments

- x:

  Adjusted p-values (legend breaks)

## Value

Plain decimals ("0.05", "0.001") when every value is at least 0.001 (or
0); otherwise all in powers of ten ("10^-5", "2.5 x 10^-4"), as an
expression vector, so that one legend never mixes the two
