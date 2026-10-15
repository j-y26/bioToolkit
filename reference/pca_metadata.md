# Sample covariates in the order of the samples

Sample covariates in the order of the samples

## Usage

``` r
pca_metadata(metadata, samples, required = FALSE)
```

## Arguments

- metadata:

  Data frame with sample names as row names, or NULL

- samples:

  Sample names

- required:

  Whether \`metadata\` must be given

## Value

A data frame with one row per sample (no columns when \`metadata\` is
NULL)
