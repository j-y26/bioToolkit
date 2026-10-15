# Cut a density (violin) at the top of the value axis

Cut a density (violin) at the top of the value axis

## Usage

``` r
cut_density(d, top)
```

## Arguments

- d:

  Rows of one violin after the stat (\`y\` sorted or not, with
  \`violinwidth\`), in the orientation of \`flip_data()\`

- top:

  Top of the axis, on the scale drawn

## Value

The rows up to \`top\`, closed by a row at \`top\` with the width
interpolated there
