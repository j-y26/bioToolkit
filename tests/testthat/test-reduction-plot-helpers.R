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

test_that("resolve_continuous_palette resolves every kind of palette", {
    expect_equal(resolve_continuous_palette(NULL, FALSE), bt_palette("plasma_r"))
    expect_equal(resolve_continuous_palette(NULL, FALSE, default = "viridis_r"), bt_palette("viridis_r"))
    expect_equal(resolve_continuous_palette(NULL, TRUE), bt_palette("blue-magenta"))
    expect_equal(resolve_continuous_palette("blue-purple", TRUE), bt_palette("blue-purple"))
    expect_length(resolve_continuous_palette("mako", FALSE), 256)
    # Low / high, and low / mid / high colours
    expect_equal(resolve_continuous_palette(c("grey90", "red"), FALSE), c("grey90", "red"))
    expect_equal(resolve_continuous_palette(c(low = "blue", mid = "white", high = "red"), TRUE),
                 c("blue", "white", "red"))
    expect_error(resolve_continuous_palette("not-a-palette", FALSE), "Unknown palette")
    expect_error(resolve_continuous_palette("npg", FALSE), "qualitative")
    expect_error(resolve_continuous_palette(c("red", "notacolour"), FALSE), "Invalid colour")
    skip_if_not_installed("RColorBrewer")
    expect_length(resolve_continuous_palette("Reds", FALSE), 9)
})

test_that("default palette starts light and coloured, not grey", {
    pal <- bt_palette("plasma_r")
    low <- grDevices::col2rgb(pal[1])[, 1]
    high <- grDevices::col2rgb(pal[256])[, 1]
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
    expect_equal(legend_title_label("auto", "auto", "Epcam", features, is_gene, NULL), quote(bolditalic("Epcam")))
    expect_equal(legend_title_label("auto", "auto", "S.Score", features, is_gene, "(smoothed)"),
                 "S.Score (smoothed)")
    expect_equal(legend_title_label("auto", "bold", "Epcam", features, is_gene, NULL), "Epcam")
    expect_equal(legend_title_label("auto", "auto", features, features, is_gene, NULL), "Value")
    expect_equal(legend_title_label("auto", "auto", "Epcam", "Epcam", is_gene, NULL), quote(bolditalic("Epcam")))
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

test_that("fill_holes fills enclosed gaps but not bays open to the edge", {
    ring <- matrix(FALSE, 7, 7)
    ring[2:6, 2:6] <- TRUE
    ring[4, 4] <- FALSE
    expect_true(all(fill_holes(ring)[2:6, 2:6]))

    bay <- ring
    bay[4, 4:6] <- FALSE
    expect_false(any(fill_holes(bay)[4, 4:6]))
    expect_false(any(fill_holes(ring)[1, ]))
})

test_that("cloud_outline circles a cloud at the requested distance", {
    set.seed(1)
    # A uniform disk of radius 1 inch
    r <- sqrt(stats::runif(3000))
    a <- stats::runif(3000, 0, 2 * pi)
    u <- r * cos(a)
    v <- r * sin(a)
    out <- cloud_outline(u, v, expand = 0.1, smooth = 0.05, coverage = 0.98)
    expect_equal(length(unique(out$ring)), 1)
    expect_equal(stats::median(sqrt(out$u^2 + out$v^2)), 1.1, tolerance = 0.05)

    # Stray cells far from the cloud are left outside
    with_strays <- cloud_outline(c(u, 3, 3.1), c(v, 3, -3), expand = 0.1, smooth = 0.05, coverage = 0.98)
    expect_equal(length(unique(with_strays$ring)), 1)
    expect_lt(max(with_strays$u), 1.3)

    # Two separate clouds of one group get one outline each
    two <- cloud_outline(c(u, u + 4), c(v, v), expand = 0.1, smooth = 0.05, coverage = 0.98)
    expect_equal(length(unique(two$ring)), 2)
})

test_that("mark_paths works in data units and per panel", {
    set.seed(2)
    df <- data.frame(
        x = c(stats::rnorm(500, 0, 0.5), stats::rnorm(500, 10, 0.5)),
        y = c(stats::rnorm(500, 0, 0.5), stats::rnorm(500, 0, 0.5)),
        mark = factor(rep(c("a", "b"), each = 500)),
        panel = factor(rep(c("p1", "p2"), 500))
    )
    paths <- mark_paths(df, by_panel = TRUE, xlim = c(-5, 15), ylim = c(-5, 5),
                        panel.width = 4, panel.height = 2,
                        expand = 0.05, smooth = 0.05, coverage = 0.98)
    expect_setequal(unique(paste(paths$mark, paths$panel)), c("a p1", "a p2", "b p1", "b p2"))
    expect_true(all(abs(paths$x[paths$mark == "a"]) < 4))
    expect_true(all(paths$x[paths$mark == "b"] > 6))
    # Every ring of every panel is its own path
    expect_equal(length(unique(paths$ring)), 4)
})

test_that("group_colours names colours by group, keeping all levels", {
    g <- factor(c("b", "a"), levels = c("a", "b", "c"))
    out <- group_colours(NULL, g)
    # Unused levels keep a colour, so that subsets are coloured like the whole
    expect_equal(levels(out$group), c("a", "b", "c"))
    expect_equal(out$cols, stats::setNames(bt_palette("scanpy-20", 3), c("a", "b", "c")))
    expect_equal(group_colours(c("red", "blue", "green"), g)$cols, c(a = "red", b = "blue", c = "green"))
    expect_equal(unname(group_colours("npg", g)$cols), bt_palette("npg", 3))
    # Colours picked by index and named
    picked <- stats::setNames(bt_palette("scanpy-20")[c(5, 1, 2)], c("c", "a", "b"))
    expect_equal(group_colours(picked, g)$cols, picked)
    expect_error(group_colours(c("red"), g, arg = "mark.cols"), "`mark.cols` has 1")
    expect_error(group_colours(c(a = "red"), g), "No colour given")
    expect_error(group_colours("mako", g), "sequential palette")
})


test_that("arrow_panel_index picks the panels of each edge of a grid", {
    # 5 panels in 3 columns, by row: row 1 = 1 2 3, row 2 = 4 5
    expect_equal(arrow_panel_index(5, 3, 2, TRUE, "bottom.left"), 4)
    expect_equal(arrow_panel_index(5, 3, 2, TRUE, "all"), 1:5)
    expect_equal(arrow_panel_index(5, 3, 2, TRUE, "none"), integer(0))
    expect_equal(arrow_panel_index(5, 3, 2, TRUE, "left"), c(1, 4))
    expect_equal(arrow_panel_index(5, 3, 2, TRUE, "right"), c(3, 5))
    expect_equal(arrow_panel_index(5, 3, 2, TRUE, "top"), 1:3)
    # Panel 3 has nothing below it: it is on the bottom edge
    expect_equal(arrow_panel_index(5, 3, 2, TRUE, "bottom"), 3:5)
    # By column: col 1 = 1 2, col 2 = 3 4, col 3 = 5
    expect_equal(arrow_panel_index(5, 3, 2, FALSE, "bottom.left"), 2)
    expect_equal(arrow_panel_index(5, 3, 2, FALSE, "right"), c(4, 5))
    expect_equal(arrow_panel_index(5, 3, 2, FALSE, "bottom"), c(2, 4, 5))
    expect_equal(panel_layout(5, 3, TRUE, 3)$by.row, TRUE)
})

test_that("format_sc_plots places arrows by arrow.panels and replaces existing ones", {
    mk <- function(i) {
        p <- ggplot2::ggplot(data.frame(x = 1:3, y = 1:3), ggplot2::aes(x, y)) + ggplot2::geom_point()
        add_dimension_arrows(p)
    }
    n_arrows <- function(p) {
        sum(vapply(p$layers, function(l) inherits(l$geom_params$grob, "bt_dimension_arrows"), logical(1)))
    }
    plots <- lapply(1:4, mk)
    out <- format_sc_plots(plots, ncol = 2)
    # bottom.left by default: one set, on plot 3; the arrows already there are
    # replaced, not added to
    expect_equal(vapply(patchwork_parts(out)$plots, n_arrows, numeric(1)), c(0, 0, 1, 0))
    out_l <- format_sc_plots(plots, ncol = 2, arrow.panels = "left")
    expect_equal(vapply(patchwork_parts(out_l)$plots, n_arrows, numeric(1)), c(1, 0, 1, 0))
    out_all <- format_sc_plots(plots, ncol = 2, arrow.panels = "all")
    expect_equal(vapply(patchwork_parts(out_all)$plots, n_arrows, numeric(1)), c(1, 1, 1, 1))
    # The old arguments are translated, with a warning
    expect_warning(old <- format_sc_plots(plots, ncol = 2, add.arrow = TRUE, bottom.left.only = FALSE),
                   "arrow.panels")
    expect_equal(vapply(patchwork_parts(old)$plots, n_arrows, numeric(1)), c(1, 1, 1, 1))
    expect_error(suppressWarnings(format_sc_plots(plots, add.arrow = FALSE, arrow.panels = "all")), "not both")
})
