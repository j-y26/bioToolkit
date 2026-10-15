# Variance fractions of each feature, from linear models

For each feature, the share of its variance explained by each covariate
(the LMG / Shapley decomposition of R^2: the increase in R^2 when the
covariate is added, averaged over all orders in which the covariates can
enter the model) and the residual fraction (1 - R^2). All features are
fitted at once, one QR decomposition per subset of covariates.

## Usage

``` r
variance_fractions(m, meta, covariates)
```

## Arguments

- m:

  Features x samples matrix

- meta:

  Data frame of covariates, rows in the order of the columns of \`m\`

- covariates:

  Columns of \`meta\`

## Value

A data frame, features x (covariates, Residuals); each row sums to 1

## References

Lindeman, Merenda & Gold (1980) Introduction to Bivariate and
Multivariate Analysis. Groemping (2006) Relative importance for linear
regression in R: the package relaimpo. J Stat Softw 17(1)
