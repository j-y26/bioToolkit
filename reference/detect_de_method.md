# Detect the tool that produced a DE table from its column names

Detect the tool that produced a DE table from its column names

## Usage

``` r
detect_de_method(data, cols, where)
```

## Arguments

- data:

  A DE table

- cols:

  Columns given explicitly

- where:

  Location for messages

## Value

A method name, or "none" when the columns were all given explicitly
