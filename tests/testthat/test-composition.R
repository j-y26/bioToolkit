# Composition scatter (S3)

comp_meta <- function() {
    # Two conditions, two samples each, two stages; known counts
    data.frame(
        ct = c(rep("A", 6), rep("B", 2), rep("A", 2), rep("B", 6), "C"),
        geno = c(rep("WT", 8), rep("Mut", 9)),
        sample = c("w1", "w1", "w1", "w2", "w2", "w2", "w1", "w2", "m1", "m2", "m1", "m1", "m1", "m2", "m2", "m2", "m2"),
        stage = "P0",
        row.names = paste0("c", 1:17)
    )
}

test_that("summarize_group_meta counts every group in every condition", {
    out <- summarize_group_meta(comp_meta(), group.by = "ct", split.by = "geno")
    wt <- out[out$condition == "WT", ]
    expect_equal(wt$count, c(6, 2, 0))
    expect_equal(wt$prop, c(6, 2, 0) / 8)
    mut <- out[out$condition == "Mut", ]
    expect_equal(mut$prop, c(2, 6, 1) / 9)
    by_sample <- summarize_group_meta(comp_meta(), "ct", "geno", facet.by = "stage", sample.by = "sample")
    w1 <- by_sample[by_sample$sample == "w1" & by_sample$condition == "WT", ]
    expect_equal(w1$prop, c(3, 1, 0) / 4)
    expect_true(all(c("facet", "sample") %in% names(by_sample)))
})

test_that("composition_wide: pooled proportions, or the mean and SEM of samples", {
    tab <- summarize_group_meta(comp_meta(), "ct", "geno")
    w <- composition_wide(tab, "WT", "Mut", "percent")
    expect_equal(w$x_value, c(6, 2, 0) / 8)
    expect_equal(w$y_count, c(2, 6, 1))
    ts <- summarize_group_meta(comp_meta(), "ct", "geno", sample.by = "sample")
    ws <- composition_wide(ts, "WT", "Mut", "percent")
    a_w <- c(3 / 4, 3 / 4)
    expect_equal(ws$x_value[ws$group == "A"], mean(a_w))
    m <- c(m1 = 1 / 4, m2 = 1 / 5)
    expect_equal(ws$y_value[ws$group == "A"], mean(m))
    expect_equal(ws$y_sem[ws$group == "A"], stats::sd(m) / sqrt(2))
    expect_equal(ws$y_n_samples[1], 2)
})

test_that("plot_composition_scatter picks the axes and handles zeros on log axes", {
    skip_if_not_installed("ggrepel")
    meta <- comp_meta()
    expect_message(p <- plot_composition_scatter(meta, "ct", "geno", x = "WT", y = "Mut"), "x: WT, y: Mut")
    src <- attr(p, "source_data")
    expect_equal(src$x_value, c(6, 2, 0) / 8)
    expect_equal(as.numeric(p$theme$panel.heights), as.numeric(p$theme$panel.widths))
    expect_error(plot_composition_scatter(meta, "ct", "geno", x = "KO"), "not a level")
    expect_message(plot_composition_scatter(meta, "ct", "geno", x = "WT", y = "Mut", log = TRUE),
                   "cannot be drawn on log axes: C")
    q <- suppressMessages(plot_composition_scatter(meta, "ct", "geno", groups = c("A", "B"), label = "A"))
    expect_equal(nrow(attr(q, "source_data")), 2)
    expect_no_error(ggplot2::ggplotGrob(q))
})
