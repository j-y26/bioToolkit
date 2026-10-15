# Check that colours stay distinct with colour vision deficiency

Simulates deuteranopia and protanopia (red-green colour blindness,
together about 8% of men of European descent) and, optionally,
tritanopia, with colorspace, and flags pairs of colours that become hard
to tell apart. Useful after picking colours for cell types.

Distances are Euclidean distances in CIELAB (Delta E 1976). About 2.3 is
just noticeable side by side; colours of small points or thin lines need
far more to be told apart, hence the default of 10. Colours already
close with normal vision are reported too.

## Usage

``` r
check_palette(cols, min.distance = 10, cvd = c("deutan", "protan"))
```

## Arguments

- cols:

  A palette name or a vector of colours

- min.distance:

  Pairs closer than this are flagged, default is 10

- cvd:

  Deficiencies to simulate, any of "deutan" and "protan" (default, both)
  and "tritan"

## Value

A data frame of the flagged pairs (\`colour1\`, \`colour2\`, \`vision\`,
\`distance\`), invisibly. A warning lists them when there are any

## Examples

``` r
if (FALSE) { # \dontrun{
check_palette("okabe-ito")
check_palette(bt_palette("scanpy-20")[1:8])
} # }
```
