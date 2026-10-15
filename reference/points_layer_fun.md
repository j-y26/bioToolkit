# A function wrapping point layers so that only they are rasterized

A function wrapping point layers so that only they are rasterized

## Usage

``` r
points_layer_fun(raster, n_cells, threshold, dpi)
```

## Arguments

- raster:

  NULL to rasterize above \`threshold\` cells, or TRUE / FALSE

- n_cells:

  Number of cells plotted

- threshold:

  Number of cells above which \`raster = NULL\` rasterizes

- dpi:

  Resolution of the rasterized points

## Value

A function taking and returning a layer
