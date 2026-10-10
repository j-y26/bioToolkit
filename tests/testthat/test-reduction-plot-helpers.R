test_that("panel_layout fills by row and by column", {
    by_row <- panel_layout(5, ncol = 3, by.row = TRUE, default_ncol = 5)
    expect_equal(by_row[c("ncol", "nrow", "dir", "bottom_left")],
                 list(ncol = 3L, nrow = 2L, dir = "h", bottom_left = 4L))

    by_col <- panel_layout(5, ncol = 3, by.row = FALSE, default_ncol = 5)
    expect_equal(by_col[c("ncol", "nrow", "dir", "bottom_left")],
                 list(ncol = 3L, nrow = 2L, dir = "v", bottom_left = 2L))

    # Filling 4 panels by column in 2 rows leaves only 2 columns
    expect_equal(panel_layout(4, ncol = 3, by.row = FALSE, default_ncol = 4)$ncol, 2L)
    expect_equal(panel_layout(3, ncol = NULL, by.row = TRUE, default_ncol = 3)$nrow, 1L)
    expect_equal(panel_layout(2, ncol = 10, by.row = TRUE, default_ncol = 2)$ncol, 2L)
})

test_that("seurat_feature_ncol matches Seurat's defaults", {
    expect_equal(vapply(c(1, 2, 6, 7, 9, 10), seurat_feature_ncol, integer(1)),
                 c(1L, 2L, 2L, 3L, 3L, 4L))
})

test_that("resolve_cells accepts names and indices", {
    all <- paste0("c", 1:5)
    expect_equal(resolve_cells(NULL, all), all)
    expect_equal(resolve_cells(c("c4", "c2"), all), c("c2", "c4"))
    expect_equal(resolve_cells(c(1, 3), all), c("c1", "c3"))
    expect_equal(resolve_cells(list("c1", "c5"), all), c("c1", "c5"))
    expect_warning(out <- resolve_cells(c("c1", "x"), all), "not found")
    expect_equal(out, "c1")
})

test_that("resolve_dims keeps character dims as arrow labels", {
    expect_equal(resolve_dims(c(2, 3), "umap"), list(dims = 2:3, labels = c("UMAP-2", "UMAP-3")))
    expect_equal(resolve_dims(c("A", "B"), "umap"), list(dims = 1:2, labels = c("A", "B")))
    expect_error(resolve_dims(1, "umap"), "two dimension")
})

test_that("resolve_raster_arg maps the deprecated rasterize", {
    expect_null(resolve_raster_arg(NULL, list()))
    expect_warning(r <- resolve_raster_arg(NULL, list(rasterize = "auto")), "deprecated")
    expect_null(r)
    expect_warning(r <- resolve_raster_arg(NULL, list(rasterize = TRUE)), "deprecated")
    expect_true(r)
    expect_error(resolve_raster_arg("auto", list()), "must be NULL")
})

test_that("cutoff_value handles NA, numbers and quantiles", {
    v <- c(1, 2, 3, 4, 5)
    expect_equal(cutoff_value(NA, v, min), 1)
    expect_equal(cutoff_value(NA, v, max), 5)
    expect_equal(cutoff_value(2.5, v, min), 2.5)
    expect_equal(cutoff_value("q50", v, max), 3)
    expect_true(is.na(cutoff_value("q50", numeric(0), max)))
    expect_error(cutoff_value("p50", v, max), "quantiles")
    expect_error(cutoff_value("q150", v, max), "between")
})

test_that("resolve_na_cutoff separates zeros only for non-negative data", {
    expect_equal(resolve_na_cutoff("auto", c(0, 1, 2), "point"), 0)
    expect_null(resolve_na_cutoff("auto", c(-1, 1), "point"))
    expect_null(resolve_na_cutoff("auto", c(0, 1), "density"))
    expect_null(resolve_na_cutoff(NULL, c(0, 1), "point"))
    expect_equal(resolve_na_cutoff(0.5, c(0, 1), "point"), 0.5)
    expect_error(resolve_na_cutoff("zero", c(0, 1), "point"), "auto")
})

test_that("resolve_feature_palette resolves every kind of palette", {
    presets <- feature_palette_presets()
    expect_equal(resolve_feature_palette(NULL, FALSE), presets$bioToolkit)
    expect_equal(resolve_feature_palette(NULL, TRUE), presets$`blue-magenta`)
    expect_equal(resolve_feature_palette("blue-purple", TRUE), presets$`blue-purple`)
    expect_length(resolve_feature_palette("mako", FALSE), 256)
    expect_equal(resolve_feature_palette("G", FALSE), resolve_feature_palette("mako", FALSE))
    expect_equal(resolve_feature_palette(c("grey90", "red"), FALSE), c("grey90", "red"))
    expect_length(resolve_feature_palette(function(n) rep("red", n), FALSE), 256)
    expect_error(resolve_feature_palette("not-a-palette", FALSE), "Unknown palette")
    expect_error(resolve_feature_palette(c("red", "notacolour"), FALSE), "Invalid colour")
    skip_if_not_installed("RColorBrewer")
    expect_length(resolve_feature_palette("Reds", FALSE), 9)
})

test_that("default palette starts light and coloured, not grey", {
    low <- grDevices::col2rgb(feature_palette_presets()$bioToolkit[1])[, 1]
    high <- grDevices::col2rgb(feature_palette_presets()$bioToolkit[256])[, 1]
    expect_gt(max(low) - min(low), 100)
    expect_gt(sum(low), sum(high))
})

test_that("smooth_values keeps units and wkde follows the weights", {
    set.seed(1)
    x <- stats::runif(400)
    y <- stats::runif(400)
    lim <- c(-0.1, 1.1)
    smooth <- smooth_values(x, y, rep(2, 400), lim, lim, c(0.1, 0.1), "smooth", n_grid = 64)
    expect_equal(smooth, rep(2, 400), tolerance = 1e-8)

    w <- ifelse(x > 0.5, 1, 0)
    wkde <- smooth_values(x, y, w, lim, lim, c(0.05, 0.05), "wkde", n_grid = 64)
    expect_gt(mean(wkde[x > 0.7]), 10 * mean(wkde[x < 0.3]))
})

test_that("fitted_breaks keep labels apart on short bars", {
    breaks <- fitted_breaks(1.5, 11)(c(0.69, 1.95))
    expect_gte(length(breaks), 2)
    expect_lte(length(breaks), 3)
    expect_true(all(breaks >= 0.69 & breaks <= 1.95))
    expect_gte(length(fitted_breaks(1.5, 11)(c(1.05, 1.97))), 2)
})

test_that("legend_title_label builds automatic and custom titles", {
    is_gene <- c(Epcam = TRUE, S.Score = FALSE)
    features <- names(is_gene)
    expect_null(legend_title_label(NULL, "auto", "Epcam", features, is_gene, NULL))
    expect_equal(legend_title_label("auto", "auto", "Epcam", features, is_gene, NULL), quote(italic("Epcam")))
    expect_equal(legend_title_label("auto", "auto", "S.Score", features, is_gene, "(smoothed)"),
                 "S.Score (smoothed)")
    expect_equal(legend_title_label("auto", "bold", "Epcam", features, is_gene, NULL), "Epcam")
    expect_equal(legend_title_label("auto", "auto", features, features, is_gene, NULL), "Value")
    expect_equal(legend_title_label("auto", "auto", "Epcam", "Epcam", is_gene, NULL), quote(italic("Epcam")))
    expect_equal(legend_title_label(c(S.Score = "S", Epcam = "E"), "auto", "Epcam", features, is_gene, NULL), "E")
    expect_equal(legend_title_label(c("E", "S"), "auto", "S.Score", features, is_gene, NULL), "S")
})

test_that("strip_face sets gene names in italics", {
    expect_equal(strip_face("bold", TRUE, FALSE), "bold.italic")
    expect_equal(strip_face("plain", TRUE, FALSE), "italic")
    expect_equal(strip_face("bold", FALSE, FALSE), "bold")
    expect_equal(strip_face("bold", TRUE, TRUE), "bold")
})

test_that("overflow_fit enlarges the device only when a legend spills", {
    expect_equal(overflow_fit(legend = 2, centre = 5, total = 10, margin = 0.1), 10)
    # Spills 0.6 past the left edge with the margin, enlarged on both sides
    expect_equal(overflow_fit(legend = 3, centre = 1, total = 10, margin = 0.1), 11.2)
})
