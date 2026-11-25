# Stop if a column is not categorical

Factors, characters and logicals are categorical; so are whole numbers
(cluster or chromosome numbers). Other numbers are measurements.

## Usage

``` r
check_categorical(values, column)
```

## Arguments

- values:

  The column's values

- column:

  The column's name, for the message

## Value

\`values\`, invisibly
