# Look up a palette by name, ignoring case

A name ending in "\_r" (matplotlib's convention) is the palette of the
name without it, reversed: "viridis_r", "plasma_r", "npg_r"... Reversed
palettes are not stored, so every palette has one.

## Usage

``` r
find_palette(name)
```

## Arguments

- name:

  Palette name

## Value

The palette entry (colours reversed for an "\_r" name), or NULL when
there is none
