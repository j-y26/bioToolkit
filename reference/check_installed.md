# Check if a package is installed, with informative error message

Check if a package is installed, with informative error message

## Usage

``` r
check_installed(pkg, reason = NULL)
```

## Arguments

- pkg:

  Character string, the package name.

- reason:

  Optional character string explaining why the package is needed.

## Value

Invisible TRUE if package is available; otherwise stops with an error.
