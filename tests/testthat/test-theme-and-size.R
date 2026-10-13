# Shared theme (F3), fixed-size engine (F2) and source data (F6)

test_that("theme sizes are fixed and scale with a custom base size", {
    expect_equal(theme_sizes(), list(base = 11, small = 9, large = 12, line = 0.5, axis = 1))
    scaled <- theme_sizes(base.size = 22)
    expect_equal(scaled[c("base", "small", "large")], list(base = 22, small = 18, large = 24))
    # Line widths do not scale with the text
    expect_equal(scaled[c("line", "axis")], list(line = 0.5, axis = 1))
})

test_that("line widths in points convert to ggplot2 units", {
    # A linewidth of 1 is ggplot2::.pt device units of 0.75 pt each
    expect_equal(pt_to_linewidth(0.75 * ggplot2::.pt), 1)
    expect_equal(pt_to_size(ggplot2::.pt), 1)
})

test_that("theme_bioToolkit applies the theme sizes", {
    th <- theme_bioToolkit("axis")
    expect_equal(th$axis.text$size, 9)
    expect_equal(th$axis.title$size, 11)
    # 1 pt axes, while the thin lines (grid plots) stay at 0.5 pt
    expect_equal(th$axis.line$linewidth, pt_to_linewidth(1))
    expect_equal(th$axis.ticks$linewidth, pt_to_linewidth(1))
    expect_equal(theme_bioToolkit("axis", base.size = 8)$axis.text$size, 9 * 8 / 11)
    grid <- theme_bioToolkit("grid")
    # Grid plots have the same 1 pt axes
    expect_equal(grid$axis.line$linewidth, pt_to_linewidth(1))
    expect_equal(grid$axis.ticks$linewidth, pt_to_linewidth(1))
    expect_equal(grid$axis.text.x$angle, 45)
    emb <- theme_bioToolkit("embedding")
    expect_equal(emb$legend.text$size, 11)
    expect_equal(emb$strip.text$size, 12)
    # Transparent by default, opaque on request
    expect_equal(theme_bioToolkit()$plot.background$fill, "transparent")
    expect_equal(theme_bioToolkit(bg = "white")$plot.background$fill, "white")
})

test_that("italic_labels sets only genes in italics", {
    out <- italic_labels(c("Cd3e", "nCount_RNA"), c(TRUE, FALSE))
    expect_true(is.expression(out))
    expect_equal(out[[1]], quote(italic("Cd3e")))
    expect_equal(out[[2]], "nCount_RNA")
    named <- italic_labels(c("Cd3e", "x"), c(Cd3e = TRUE))
    expect_equal(named[[2]], "x")
    # Symbols that would break parse() are kept
    expect_equal(italic_labels("HLA-DRA")[[1]], quote(italic("HLA-DRA")))
})

test_that("panel_size_theme returns the expected units", {
    th <- panel_size_theme(2, 3, ncol = 2, nrow = 1, panel.spacing = 0.1)
    expect_equal(as.numeric(th$panel.widths), c(2, 2))
    expect_equal(grid::unitType(th$panel.widths), c("inches", "inches"))
    expect_equal(as.numeric(th$panel.heights), 3)
    expect_equal(as.numeric(th$panel.spacing), 0.1)
    # One width per facet column
    expect_equal(as.numeric(panel_size_theme(c(1, 2.5), 1)$panel.widths), c(1, 2.5))
    expect_error(panel_size_theme(-1), "positive")
})

test_that("a dot grid with 10 rows of 0.2 in measures 2 in tall", {
    library(ggplot2)
    df <- expand.grid(gene = paste0("g", 1:10), group = paste0("G", 1:4))
    df$size <- seq_len(nrow(df))
    p <- ggplot(df, aes(group, gene, size = size)) +
        geom_point() +
        scale_x_discrete(expand = grid_expansion()) +
        scale_y_discrete(expand = grid_expansion()) +
        theme_bioToolkit("grid") +
        panel_size_theme(grid_panel_size(4, 0.2), grid_panel_size(10, 0.2))
    # A null device, so that building the grob writes no Rplots.pdf
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    gt <- ggplotGrob(p)
    panel <- gt$layout[gt$layout$name == "panel", ]
    expect_equal(grid::convertHeight(gt$heights[panel$t], "in", valueOnly = TRUE), 2)
    expect_equal(grid::convertWidth(gt$widths[panel$l], "in", valueOnly = TRUE), 0.8)
    # Each category spans exactly one cell: the y range is 10 units
    expect_equal(diff(ggplot_build(p)$layout$panel_params[[1]]$y$continuous_range), 10)
    size <- fixed_gtable_size(gt)
    expect_gt(size[["height"]], 2)
})

test_that("save_plot sizes fixed plots and needs a size for others", {
    library(ggplot2)
    p <- ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_bioToolkit()
    f <- tempfile(fileext = ".pdf")
    expect_error(save_plot(p, f), "no fixed panel size")
    dims <- save_plot(p, f, width = 3, height = 2)
    expect_equal(unname(dims), c(3, 2))
    fixed <- p + panel_size_theme(2, 1.5)
    expect_message(dims <- save_plot(fixed, f, width = 9, height = 9), "ignored")
    expect_gt(dims[["width"]], 2)
    expect_lt(dims[["width"]], 3.5)
    expect_true(file.exists(f))
    png <- tempfile(fileext = ".png")
    save_plot(fixed, png, dpi = 100)
    skip_if_not_installed("png")
    skip_if_not_installed("jpeg")
    # The default transparent background is kept in PNG
    expect_equal(png::readPNG(png)[1, 1, 4], 0)
    # A patchwork of transparent plots stays transparent
    pw <- tempfile(fileext = ".png")
    save_plot(fixed | fixed, pw, dpi = 100)
    expect_equal(png::readPNG(pw)[1, 1, 4], 0)
    # JPEG has no transparency: the background is white, not black
    jpg <- tempfile(fileext = ".jpg")
    save_plot(fixed, jpg, dpi = 100)
    expect_gt(min(jpeg::readJPEG(jpg)[1, 1, ]), 0.95)
})

test_that("export_source_data writes one table per plot", {
    library(ggplot2)
    a <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
    b <- ggplot(mtcars, aes(hp, qsec)) + geom_point()
    attr(b, "source_data") <- data.frame(x = 1:3)
    f <- tempfile(fileext = ".csv")
    tabs <- export_source_data(a, f)
    expect_equal(nrow(utils::read.csv(f)), nrow(mtcars))
    # Row names are kept as a column
    expect_equal(names(tabs$plot)[1], "row")
    tabs <- export_source_data(a | b, f)
    expect_length(tabs, 2)
    expect_equal(tabs[[2]], data.frame(x = 1:3))
    files <- paste0(sub("\\.csv$", "", f), "_", names(tabs), ".csv")
    expect_true(all(file.exists(files)))
    expect_error(export_source_data(a, tempfile(fileext = ".txt")), "csv or .xlsx")
    skip_if_not_installed("writexl")
    x <- tempfile(fileext = ".xlsx")
    export_source_data(a | b, x)
    expect_true(file.exists(x))
})
