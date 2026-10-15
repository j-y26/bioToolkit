# Shared helpers of the 2026-10-14 revision: group bars, diverging scales,
# bold legend titles

test_that("group_bar_layer draws a bar outside the panel, spanning its categories", {
    df <- data.frame(x = factor(c("a", "b", "c")), y = 1:3, g = factor(c("G1", "G1", "G2")))
    bars <- data.frame(g = factor(c("G1", "G2")), x = factor(c("a", "c"), levels = c("a", "b", "c")),
                       xend = factor(c("b", "c"), levels = c("a", "b", "c")), fill = c("black", "red"))
    p <- ggplot2::ggplot(df, ggplot2::aes(x, y)) +
        ggplot2::geom_point() +
        group_bar_layer(bars, side = "top", width = 3.5) +
        ggplot2::facet_grid(cols = ggplot2::vars(g), scales = "free_x", space = "free_x") +
        ggplot2::coord_cartesian(clip = "off")
    expect_no_error(ggplot2::ggplotGrob(p))
    expect_equal(p$layers[[2]]$geom_params$side, "top")
    expect_error(group_bar_layer(bars, side = "middle"))
    m <- as.numeric(group_bar_margin("top", 3.5))
    expect_equal(m[3], GROUP_BAR_GAP + 3.5 + GROUP_BAR_TEXT_GAP)
    expect_error(check_bar_width(-1), "positive")
})

test_that("group_bar_colours: black by default, else group colours", {
    expect_equal(unname(group_bar_colours(NULL, c("a", "b"), "x")), c("black", "black"))
    expect_equal(group_bar_colours(c(b = "red", a = "blue"), c("a", "b"), "x")[["a"]], "blue")
    expect_equal(unname(group_bar_colours("okabe-ito", c("a", "b"), "x")), unname(bt_palette("okabe-ito", 2)))
})

test_that("diverging_limits trims each end to the data on its side", {
    # The user's example: data -3 to 10, cutoff 5 -> -3 to 5
    expect_equal(diverging_limits(c(-3, 2, 10), -5, 5), c(-3, 5))
    expect_equal(diverging_limits(c(-3, 2, 10)), c(-3, 10))
    # The centre is always inside
    expect_equal(diverging_limits(c(1, 4)), c(0, 4))
    expect_equal(diverging_limits(c(0, 0)), c(-1, 1))
    expect_equal(diverging_limits(numeric(0)), c(-1, 1))
})

test_that("rescale_diverging maps the centre to the middle and each end to its own end", {
    r <- rescale_diverging(0)
    expect_equal(r(c(-3, 0, 5), from = c(-3, 5)), c(0, 0.5, 1))
    expect_equal(r(c(-1.5, 2.5), from = c(-3, 5)), c(0.25, 0.75))
    # A side without room maps to its end colour; NA stays NA
    expect_equal(r(c(NA, 0, 2), from = c(0, 4)), c(NA, 0.5, 0.75))
})

test_that("bold_title makes plotmath titles bold", {
    expect_equal(bold_title("Genes"), "Genes")
    expect_null(bold_title(NULL))
    expect_equal(bold_title(quote(-log[10] ~ "padj")), quote(bold(-log[10] ~ "padj")))
    expect_equal(bold_title(quote(italic("Cd3e"))), quote(bolditalic("Cd3e")))
    expect_equal(bold_title(quote(bold("x"))), quote(bold("x")))
})

test_that("legend titles are bold in the theme", {
    expect_equal(theme_bioToolkit()$legend.title$face, "bold")
    expect_equal(theme_bioToolkit("embedding")$legend.title$face, "bold")
})

test_that("ggplot_colors replaces ggplotColours", {
    expect_equal(ggplot_colors(3), scales::hue_pal()(3))
    expect_warning(old <- ggplotColours(3), "ggplot_colors")
    expect_equal(old, ggplot_colors(3))
})
