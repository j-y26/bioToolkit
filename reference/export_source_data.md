# Export the data behind a plot

Writes the table each plot of the package is built from, so that it can
be submitted as source data with a figure. For a patchwork, each plot is
written separately: one CSV file per plot, or one sheet per plot of an
XLSX file.

The table is the "source_data" attribute of a plot when it has one, else
its data (\`plot\$data\`); row names (e.g. cell names) are kept as a
first column.

## Usage

``` r
export_source_data(plot, file = NULL)
```

## Arguments

- plot:

  A ggplot or patchwork object

- file:

  Output file, ".csv" or ".xlsx". For a patchwork written to CSV, the
  plot names (or numbers) are appended to the file name. If NULL
  (default), nothing is written and the tables are only returned

## Value

A named list of data frames, one per plot, invisibly

## Examples

``` r
library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
f <- tempfile(fileext = ".csv")
export_source_data(p, f)
head(read.csv(f))
#>                 row  mpg cyl disp  hp drat    wt  qsec vs am gear carb
#> 1         Mazda RX4 21.0   6  160 110 3.90 2.620 16.46  0  1    4    4
#> 2     Mazda RX4 Wag 21.0   6  160 110 3.90 2.875 17.02  0  1    4    4
#> 3        Datsun 710 22.8   4  108  93 3.85 2.320 18.61  1  1    4    1
#> 4    Hornet 4 Drive 21.4   6  258 110 3.08 3.215 19.44  1  0    3    1
#> 5 Hornet Sportabout 18.7   8  360 175 3.15 3.440 17.02  0  0    3    2
#> 6           Valiant 18.1   6  225 105 2.76 3.460 20.22  1  0    3    1
```
