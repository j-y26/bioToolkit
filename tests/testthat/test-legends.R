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
    # Base pdf(), as grob_inches() measures on it (the default PDF device,
    # cairo, measures with the system font)
    dims <- save_plot(leg, f, device = grDevices::pdf)
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


test_that("legend blocks arrange legends in columns, built when drawn", {
    p <- ggplot2::ggplot(mtcars, ggplot2::aes(wt, mpg, colour = hp, size = cyl)) + ggplot2::geom_point()
    block <- legend_block_grob(list(list(p + ggplot2::guides(size = "none")),
                                    list(p + ggplot2::guides(colour = "none"))))
    expect_s3_class(block, "bt_legend_block")
    grDevices::pdf(NULL)
    lay <- legend_block_layout(block)
    grDevices::dev.off()
    expect_length(lay$grobs, 2)
    expect_gt(lay$width, 0)
    # Fixed width only: the figure takes the height of the rows it spans
    fig <- legend_block_figure(block)
    expect_gt(as.numeric(fig$theme$panel.widths), lay$width)
    expect_null(fig$theme$panel.heights)
    expect_gt(attr(fig, "block_height"), 0)
})
