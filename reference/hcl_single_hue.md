# A single-hue sequential palette, trimmed to a lightness range

Takes one of base R's HCL sequential palettes
([`grDevices::hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html),
Zeileis et al. 2020) and keeps the part between lightness \`light\` and
\`dark\` (CIE L\*, 0 = black, 100 = white). The untrimmed palettes run
from near white, which vanishes on a white or transparent background, to
near black, which hides the hue; L\* 80 is a clear tint (dots of the
lowest values stay visible on white) and L\* 30 still clearly coloured.
Every palette trimmed this way spans the same lightness, so that
different hues encode the same value with the same lightness.

## Usage

``` r
hcl_single_hue(name, n = 9, light = 80, dark = 30)
```

## Arguments

- name:

  Name of an [`hcl.colors()`](https://rdrr.io/r/grDevices/palettes.html)
  sequential palette

- n:

  Number of colours, default is 9

- light, dark:

  Lightness (L\*) of the light and dark ends, default 80 and 30

## Value

\`n\` colours from light to dark, taken from the palette (not
interpolated)
