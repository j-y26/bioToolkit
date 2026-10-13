# extract_legend() and saving legends with save_plot()

legend_test_plot <- function(...) {
    set.seed(1)
    de <- data.frame(
        logFC = rnorm(300, sd = 1.5),
        PValue = runif(300)^4,
        row.names = paste0("Gene", 1:300)
    )
    de$FDR <- p.adjust(de$PValue, "BH")
    suppressMessages(plot_volcano(de, method = "edgeR", label = "none", ...))
}

# Width and height of a grob in inches, on a PDF device
grob_inches <- function(g) {
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off())
    c(grid::convertWidth(grid::grobWidth(g), "in", valueOnly = TRUE),
      grid::convertHeight(grid::grobHeight(g), "in", valueOnly = TRUE))
}

test_that("by default the legend is exactly as drawn in the plot", {
    p <- legend_test_plot()
    leg <- extract_legend(p)
    expect_s3_class(leg, "bt_legend")
    expect_s3_class(leg, "gtable")
    grDevices::pdf(NULL)
    gt <- ggplot2::ggplotGrob(p)
    grDevices::dev.off()
    in_plot <- gt$grobs[[grep("guide-box-right", gt$layout$name)]]
    expect_equal(grob_inches(leg), grob_inches(in_plot))
})

test_that("direction and theme settings change the legend", {
    p <- legend_test_plot()
    vertical <- grob_inches(extract_legend(p))
    horizontal <- grob_inches(extract_legend(p, direction = "horizontal"))
    # Three keys in a row: wider and less tall than in a column
    expect_gt(horizontal[1], vertical[1])
    expect_lt(horizontal[2], vertical[2])
    bigger <- grob_inches(extract_legend(p, legend.text = ggplot2::element_text(size = 20)))
    expect_gt(bigger[1], vertical[1])
    expect_error(extract_legend(p, direction = "diagonal"), "should be one of")
})

test_that("hidden legends are extracted and plots without legends error", {
    # annotate.counts = "arrow" hides the legend
    leg <- extract_legend(legend_test_plot(annotate.counts = "arrow"))
    expect_s3_class(leg, "bt_legend")
    no_legend <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg)) + ggplot2::geom_point()
    expect_error(extract_legend(no_legend), "no legend")
    expect_error(extract_legend(no_legend | no_legend), "not a patchwork")
})

test_that("save_plot saves a legend sized to fit, with a margin", {
    p <- legend_test_plot()
    leg <- extract_legend(p, direction = "horizontal")
    pad <- 2 * theme_sizes()$small / 72
    f <- tempfile(fileext = ".pdf")
    dims <- save_plot(leg, f)
    expect_true(file.exists(f))
    expect_equal(unname(dims), grob_inches(leg) + pad, tolerance = 1e-6)
    # PNG: re-built with the PNG device's fonts, so its size may differ from
    # the PDF's, but the file has exactly the measured size
    skip_if_not_installed("png")
    g <- tempfile(fileext = ".png")
    dims_png <- save_plot(leg, g, dpi = 100)
    img <- png::readPNG(g)
    expect_equal(dim(img)[2:1], round(unname(dims_png) * 100), tolerance = 1)
    # Transparent like the plot it came from
    expect_lt(min(img[, , 4]), 1)
})

test_that("the volcano panel is square by default", {
    p <- legend_test_plot(panel.width = 2.5)
    expect_equal(p$theme$panel.widths, grid::unit(2.5, "in"))
    expect_equal(p$theme$panel.heights, grid::unit(2.5, "in"))
})
