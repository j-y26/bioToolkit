# List the palettes of the package

List the palettes of the package

## Usage

``` r
list_palettes(type = NULL)
```

## Arguments

- type:

  Types to list: any of "qualitative" (for groups), "sequential" and
  "diverging". If NULL (default), all

## Value

A data frame with the palette \`name\`, its \`type\`, its number of
colours \`n\` and its \`source\`, for citing

## Examples

``` r
list_palettes("qualitative")
#>                 name        type   n
#> 1          scanpy-20 qualitative  20
#> 2       scanpy-20-cb qualitative  20
#> 3         scanpy-102 qualitative 102
#> 4          okabe-ito qualitative   8
#> 5         tol-bright qualitative   6
#> 6        tol-vibrant qualitative   6
#> 7          tol-muted qualitative   9
#> 8          tol-light qualitative   8
#> 9  tol-high-contrast qualitative   3
#> 10               npg qualitative  10
#> 11              aaas qualitative  10
#> 12              nejm qualitative   8
#> 13            lancet qualitative   9
#> 14              jama qualitative   7
#> 15               jco qualitative  10
#> 16        tableau-10 qualitative  10
#> 17     polychrome-36 qualitative  36
#> 18          alphabet qualitative  26
#> 19            paired qualitative  12
#>                                                                                                                                                                             source
#> 1  matplotlib 'tab10' (the colours scanpy uses for up to 10 groups) followed by the 10 lighter 'tab20' tones; D3 category10 / category20 (M. Bostock). Default for up to 20 groups
#> 2                                       scanpy default_20 for 11-20 groups (Wolf et al. 2018, Genome Biol 19:15): tab20 without greys, three colours adjusted for colour-blindness
#> 3                                                                           scanpy default_102 ('godsnot_102', maximally distinct colours, black removed). Default above 20 groups
#> 4                                                                                     Okabe & Ito (2008), Color Universal Design; Wong (2011) Nat Methods 8:441. Colour-blind safe
#> 5                                                                         Paul Tol (2021), Colour schemes, SRON/EPS/TN/09-002, https://personal.sron.nl/~pault/. Colour-blind safe
#> 6                                                                         Paul Tol (2021), Colour schemes, SRON/EPS/TN/09-002, https://personal.sron.nl/~pault/. Colour-blind safe
#> 7                                                                         Paul Tol (2021), Colour schemes, SRON/EPS/TN/09-002, https://personal.sron.nl/~pault/. Colour-blind safe
#> 8                                           Paul Tol (2021), Colour schemes, SRON/EPS/TN/09-002, https://personal.sron.nl/~pault/. Colour-blind safe; light, for annotation strips
#> 9                                                  Paul Tol (2021), Colour schemes, SRON/EPS/TN/09-002, https://personal.sron.nl/~pault/. Colour-blind safe; distinct in greyscale
#> 10                                                                          ggsci (N. Xiao), after the colours of Nature Publishing Group figures; not an official journal palette
#> 11                                                                                   ggsci (N. Xiao), after the colours of Science (AAAS) figures; not an official journal palette
#> 12                                                                  ggsci (N. Xiao), after the colours of New England Journal of Medicine figures; not an official journal palette
#> 13                                                                                       ggsci (N. Xiao), after the colours of The Lancet figures; not an official journal palette
#> 14                                                                                             ggsci (N. Xiao), after the colours of JAMA figures; not an official journal palette
#> 15                                                                     ggsci (N. Xiao), after the colours of Journal of Clinical Oncology figures; not an official journal palette
#> 16                                                                                                                Tableau 10 (Tableau Software), as in grDevices::palette.colors()
#> 17                                                                              Coombes et al. (2019) J Stat Softw 90(1), as in grDevices::palette.colors(); Seurat's 'polychrome'
#> 18                                                                                Green-Armytage (2010), a colour alphabet, as in grDevices::palette.colors(); Seurat's 'alphabet'
#> 19                                                                                         ColorBrewer 'Paired' (C. Brewer), as in grDevices::palette.colors(); light / dark pairs
```
