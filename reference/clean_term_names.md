# Readable names for MSigDB-style gene set names

Readable names for MSigDB-style gene set names

## Usage

``` r
clean_term_names(x)
```

## Arguments

- x:

  Gene set names

## Value

Names without the collection prefix and underscores, in sentence case,
for names in upper case with underscores. Other names get a capital
first letter when their first word is all lower case (e.g. GO
descriptions, "extracellular matrix organization"); words with capitals
inside (e.g. "mRNA") are left as they are
