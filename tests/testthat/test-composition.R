# Composition scatter

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

test_that("plot_composition_scatter leaves out groups absent from a panel", {
    skip_if_not_installed("ggrepel")
    meta <- data.frame(
        ct = c("A", "A", "B", "B", "C", "C"),
        geno = c("WT", "Mut", "WT", "Mut", "WT", "Mut"),
        lineage = c("L1", "L1", "L1", "L1", "L2", "L2")
    )
    p <- suppressMessages(plot_composition_scatter(meta, "ct", "geno", facet.by = "lineage"))
    src <- attr(p, "source_data")
    expect_equal(nrow(src), 3)
    expect_setequal(paste(src$facet, src$group), c("L1 A", "L1 B", "L2 C"))
})


test_that("plot_composition_scatter draws SEM, SD or 95% CI bars with caps", {
    meta <- comp_meta()
    p <- suppressMessages(plot_composition_scatter(meta, "ct", "geno", x = "WT", y = "Mut",
                                                   sample.by = "sample", label = FALSE))
    src <- attr(p, "source_data")
    expect_equal(src$y_error, src$y_sem)
    bars <- Filter(function(l) inherits(l$geom, "GeomErrorbar"), p$layers)
    expect_length(bars, 2)
    # Caps (a width) and 0.75 pt black lines
    expect_true(all(vapply(bars, function(l) l$geom_params$width > 0, logical(1))))
    expect_equal(bars[[1]]$aes_params$colour, "black")
    expect_equal(bars[[1]]$aes_params$linewidth, pt_to_linewidth(0.75))
    psd <- suppressMessages(plot_composition_scatter(meta, "ct", "geno", x = "WT", y = "Mut",
                                                     sample.by = "sample", error = "sd", label = FALSE))
    expect_equal(attr(psd, "source_data")$x_error, attr(psd, "source_data")$x_sd)
    pci <- suppressMessages(plot_composition_scatter(meta, "ct", "geno", x = "WT", y = "Mut",
                                                     sample.by = "sample", error = "ci95", label = FALSE))
    s2 <- attr(pci, "source_data")
    expect_equal(s2$y_error, stats::qt(0.975, 1) * s2$y_sem)
    expect_no_error(ggplot2::ggplotGrob(pci))
})


# Table input and an empty split.by (composition bars)

test_that("summarize_group_meta takes tables and an empty split.by", {
    m <- comp_meta()
    tb <- tibble::as_tibble(m)
    out <- summarize_group_meta(tb, group.by = "ct", split.by = "geno")
    expect_equal(out$count[out$condition == "WT"], c(6, 2, 0))
    all <- summarize_group_meta(m, group.by = "ct", unit = "genes")
    expect_equal(levels(all$condition), "All genes")
    expect_equal(all$count, c(8, 8, 1))
    dup <- m
    rownames(dup) <- NULL
    dup <- rbind(dup, dup)
    expect_equal(sum(summarize_group_meta(dup, "ct")$count), 34)
    expect_error(summarize_group_meta(dup, "ct", cells = c("c1", "c2")), "row indices")
})

test_that("summarize_group_meta checks categories and names the unit", {
    m <- comp_meta()
    m$score <- seq(0.5, 8.5, by = 0.5)
    expect_error(summarize_group_meta(m, "score"), "not categorical")
    m$cluster <- rep(c(10L, 2L, 1L), length.out = 17)
    expect_equal(levels(summarize_group_meta(m, "cluster")$group), c("1", "2", "10"))
    m$geno[1:2] <- NA
    expect_message(summarize_group_meta(m, "ct", "geno", unit = "genes"), "2 genes without a value")
})

test_that("summarize_group_meta takes anything as.data.frame() converts", {
    # A character matrix goes through the same conversion as a Bioconductor
    # DataFrame or a GRanges
    mat <- as.matrix(comp_meta()[c("ct", "geno")])
    expect_equal(nrow(summarize_group_meta(mat, "ct", "geno")), 6)
    expect_error(summarize_group_meta(function(x) x, "ct"), "as.data.frame")
})

test_that("composition_segments stacks bars and places them", {
    tab <- summarize_group_meta(comp_meta(), "ct", "geno")
    seg <- composition_segments(tab, "percent")
    wt <- seg[seg$bar == "WT", ]
    expect_equal(wt$upper[wt$group == "A"], 1)
    expect_equal(wt$lower[wt$group == "A"], 2 / 8)
    expect_equal(wt$upper[wt$group == "B"], 2 / 8)
    expect_equal(tapply(seg$value, seg$bar, sum)[["Mut"]], 1)
    left <- composition_segments(tab, "percent", first.at.start = TRUE)
    expect_equal(left$lower[left$bar == "WT" & left$group == "A"], 0)
    cnt <- composition_segments(tab, "count")
    expect_equal(max(cnt$upper[cnt$bar == "Mut"]), 9)
    frac <- composition_segments(tab, "fraction")
    expect_equal(frac$value, seg$value)
    ord <- composition_segments(tab, "percent", order.by = "B")
    expect_equal(unique(as.character(ord$bar[order(ord$pos)])), c("Mut", "WT"))
    expect_error(composition_segments(tab, "percent", order.by = "Z"), "level")
})

test_that("composition_segments gives repeated bar names their own place per panel", {
    m <- rbind(comp_meta(), transform(comp_meta(), stage = "P7"))
    rownames(m) <- paste0("c", seq_len(nrow(m)))
    tab <- summarize_group_meta(m, "ct", "sample", facet.by = "stage")
    seg <- composition_segments(tab, "percent")
    expect_equal(length(unique(seg$pos)), 8)
    expect_true(max(seg$pos[seg$facet == "P0"]) < min(seg$pos[seg$facet == "P7"]))
})

test_that("composition_bands joins neighbours within a panel only", {
    tab <- summarize_group_meta(comp_meta(), "ct", "sample", facet.by = "geno")
    seg <- composition_segments(tab, "percent")
    bands <- composition_bands(seg, half = 0.3)
    # WT panel: w1, w2 (C empty in both -> 2 bands); Mut: m1, m2 (3 bands)
    expect_equal(length(unique(bands$id)), 5)
    expect_equal(nrow(bands[bands$id == bands$id[1], ]), 80)
    straight <- composition_bands(seg, half = 0.3, shape = "straight")
    expect_equal(nrow(straight[straight$id == straight$id[1], ]), 4)
    one <- composition_segments(summarize_group_meta(comp_meta(), "ct"), "percent")
    expect_equal(nrow(composition_bands(one, half = 0.3)), 0)
})

# Layer of a plot drawn by a given geom
layer_of <- function(p, geom) {
    which(vapply(p$layers, function(l) inherits(l$geom, geom), TRUE))
}

test_that("plot_composition_bars draws stacked bars with bands", {
    p <- plot_composition_bars(comp_meta(), "ct", "sample", facet.by = "geno", connect = "curve")
    expect_s3_class(p, "ggplot")
    b <- ggplot2::ggplot_build(p)
    rects <- b$data[[layer_of(p, "GeomRect")]]
    expect_equal(nrow(rects), 12)
    src <- attr(p, "source_data")
    expect_equal(as.numeric(tapply(src$value, src$bar, sum)), rep(1, 4))
    expect_false("pos" %in% names(src))
    expect_length(layer_of(p, "GeomPolygon"), 1)
    none <- plot_composition_bars(comp_meta(), "ct", "sample", facet.by = "geno", connect = "none")
    expect_length(layer_of(none, "GeomPolygon"), 0)
})

test_that("plot_composition_bars shows percent, fraction or count on the value axis", {
    axis_labels <- function(p) {
        b <- ggplot2::ggplot_build(p)
        b$layout$panel_params[[1]]$y$get_labels()
    }
    pct <- plot_composition_bars(comp_meta(), "ct", "geno")
    expect_equal(axis_labels(pct), c("0%", "25%", "50%", "75%", "100%"))
    expect_equal(ggplot2::get_labs(pct)$y, "% of cells")
    frac <- plot_composition_bars(comp_meta(), "ct", "geno", value = "fraction")
    expect_equal(axis_labels(frac), c("0", "0.25", "0.5", "0.75", "1"))
    expect_equal(ggplot2::get_labs(frac)$y, "Fraction of cells")
    cnt <- plot_composition_bars(comp_meta(), "ct", "geno", value = "count")
    expect_equal(ggplot2::get_labs(cnt)$y, "Cells")
    expect_equal(max(attr(cnt, "source_data")$upper), 9)
})

test_that("plot_composition_bars flips, counts and labels", {
    p <- plot_composition_bars(
        comp_meta(),
        "ct",
        "geno",
        flip = TRUE,
        value = "count",
        unit = "genes",
        show.n = TRUE,
        label = TRUE
    )
    expect_equal(ggplot2::get_labs(p)$x, "Genes")
    b <- ggplot2::ggplot_build(p)
    txt <- do.call(rbind, lapply(b$data, function(d) if ("label" %in% names(d)) d[c("label", "colour")]))
    expect_true(all(c("8", "9") %in% txt$label))
    expect_false(any(txt$colour %in% c("grey50", "gray50")))
    small <- plot_composition_bars(comp_meta(), "ct", "geno", label = TRUE, label.min = 0.5)
    sb <- ggplot2::ggplot_build(small)
    seg_txt <- sb$data[[layer_of(small, "GeomText")[1]]]
    expect_equal(sort(seg_txt$label), c("67%", "75%"))
})

test_that("plot_composition_bars keeps colours of unused levels and warns on many groups", {
    m <- comp_meta()
    m$ct <- factor(m$ct, levels = c("A", "B", "C", "D"))
    cols <- c(
        A = "#1F77B4",
        B = "#FF7F0E",
        C = "#2CA02C",
        D = "#D62728"
    )
    p <- plot_composition_bars(m, "ct", "geno", cols = cols)
    expect_false("D" %in% attr(p, "source_data")$group)
    b <- ggplot2::ggplot_build(p)
    expect_true("#2CA02C" %in% b$data[[layer_of(p, "GeomRect")]]$fill)
    many <- data.frame(g = paste0("g", 1:40), s = "x")
    expect_warning(plot_composition_bars(many, "g", "s"), "40 levels")
    expect_error(plot_composition_bars(comp_meta(), "ct", "geno", bar.width = 0), "bar.width")
})

test_that("plot_composition_bars saves as PDF and PNG", {
    p <- plot_composition_bars(comp_meta(), "ct", "sample", facet.by = "geno", show.n = TRUE)
    f <- tempfile(fileext = ".pdf")
    on.exit(unlink(c(f, sub("pdf$", "png", f))))
    expect_no_warning(save_plot(p, f))
    expect_no_warning(save_plot(p, sub("pdf$", "png", f)))
})

test_that("plot_composition_bars labels only segments the text fits in", {
    p <- plot_composition_bars(
        comp_meta(),
        "ct",
        "geno",
        flip = TRUE,
        label = TRUE,
        label.min = 0,
        bar.length = 0.5
    )
    b <- ggplot2::ggplot_build(p)
    seg_txt <- b$data[[layer_of(p, "GeomText")[1]]]
    expect_equal(sort(seg_txt$label), c("67%", "75%"))
})

test_that("plot_composition_bars ends a count axis at a break above the longest bar", {
    p <- plot_composition_bars(comp_meta(), "ct", "geno", value = "count", show.n = TRUE)
    b <- ggplot2::ggplot_build(p)
    y <- b$layout$panel_params[[1]]$y
    expect_equal(tail(y$get_labels(), 1), "10")
    expect_gte(y$continuous_range[2], 10)
})

test_that("plot_composition_bars: no bands by default, black borders, legend options", {
    p <- plot_composition_bars(comp_meta(), "ct", "geno")
    expect_length(layer_of(p, "GeomPolygon"), 0)
    b <- ggplot2::ggplot_build(p)
    expect_true(all(b$data[[layer_of(p, "GeomRect")]]$colour == "black"))
    expect_equal(ggplot2::get_labs(p)$fill, "ct")
    q <- plot_composition_bars(
        comp_meta(),
        "ct",
        "geno",
        legend.title = "Cell type",
        legend.position = "bottom",
        legend.ncol = 2
    )
    expect_equal(ggplot2::get_labs(q)$fill, "Cell type")
    expect_equal(q$theme$legend.position, "bottom")
    expect_equal(q$theme$legend.title.position, "top")
    expect_equal(q$guides$guides$fill$params$ncol, 2)
    expect_error(plot_composition_bars(comp_meta(), "ct", "geno", legend.position = "middle"), "legend.position")
    none <- plot_composition_bars(comp_meta(), "ct", "geno", legend.title = "")
    expect_s3_class(none$theme$legend.title, "element_blank")
})

test_that("horizontal bars leave room for the last value label before a right legend", {
    p <- plot_composition_bars(comp_meta(), "ct", "geno", flip = TRUE)
    expect_equal(as.numeric(p$theme$legend.box.spacing), 0.25)
    v <- plot_composition_bars(comp_meta(), "ct", "geno")
    expect_false(identical(v$theme$legend.box.spacing, p$theme$legend.box.spacing))
})

# Final review fixes

test_that("plot functions keep the numbered-rows error for `cells` names", {
    tb <- tibble::as_tibble(comp_meta())
    expect_error(plot_composition_bars(tb, "ct", "geno", cells = c("c1", "c2")), "row indices")
    expect_error(plot_composition_scatter(tb, "ct", "geno", cells = c("c1", "c2")), "row indices")
})

test_that("a single-bar panel is wide enough for its name", {
    m <- comp_meta()
    m$panel <- ifelse(m$sample == "m2", "Stimulated_long_name", "Control")
    m$bar <- ifelse(m$sample == "m2", "m2", m$sample)
    p <- plot_composition_bars(m, "ct", "bar", facet.by = "panel")
    widths <- as.numeric(p$theme$panel.widths)
    # "Stimulated_long_name": 20 characters of bold 9 pt text
    expect_gte(widths[2], 20 * 0.6 * 9 / 72)
    b <- ggplot2::ggplot_build(p)
    rects <- b$data[[layer_of(p, "GeomRect")]]
    bar_w <- unique(round(rects$xmax - rects$xmin, 6))
    expect_length(bar_w, 1)
    # Only the panel's own bar is named under it
    x2 <- b$layout$panel_params[[2]]$x
    drawn <- !is.na(x2$get_breaks())
    expect_equal(x2$get_labels()[drawn], "m2")
})

test_that("a count axis has whole-number breaks", {
    tiny <- data.frame(g = c("a", "b", "a"), s = c("x", "x", "y"))
    p <- plot_composition_bars(tiny, "g", "s", value = "count")
    b <- ggplot2::ggplot_build(p)
    breaks <- b$layout$panel_params[[1]]$y$get_breaks()
    breaks <- breaks[!is.na(breaks)]
    expect_true(all(breaks == round(breaks)))
})

test_that("show.n adds room without shortening the bars", {
    plain <- plot_composition_bars(comp_meta(), "ct", "geno")
    with_n <- plot_composition_bars(comp_meta(), "ct", "geno", show.n = TRUE)
    expect_equal(as.numeric(plain$theme$panel.heights), 2.5)
    h <- as.numeric(with_n$theme$panel.heights)
    expect_gt(h, 2.5)
    b <- ggplot2::ggplot_build(with_n)
    range_y <- diff(b$layout$panel_params[[1]]$y$continuous_range)
    # One unit of the value axis (a whole bar) is still 2.5 in
    expect_equal(h / range_y, 2.5)
})

test_that("border.width sets the bar borders; the value axis line is off by default", {
    p <- plot_composition_bars(comp_meta(), "ct", "geno")
    b <- ggplot2::ggplot_build(p)
    expect_equal(unique(b$data[[layer_of(p, "GeomRect")]]$linewidth), pt_to_linewidth(1))
    thin <- plot_composition_bars(comp_meta(), "ct", "geno", border.width = 0.25)
    tb <- ggplot2::ggplot_build(thin)
    expect_equal(unique(tb$data[[layer_of(thin, "GeomRect")]]$linewidth), pt_to_linewidth(0.25))
    expect_error(plot_composition_bars(comp_meta(), "ct", "geno", border.width = -1), "border.width")
    # Vertical bars: no y-axis line, x-axis line kept; horizontal: the reverse
    expect_s3_class(p$theme$axis.line.y, "element_blank")
    expect_false(inherits(p$theme$axis.line.x, "element_blank"))
    h <- plot_composition_bars(comp_meta(), "ct", "geno", flip = TRUE)
    expect_s3_class(h$theme$axis.line.x, "element_blank")
    with_line <- plot_composition_bars(comp_meta(), "ct", "geno", value.axis.line = TRUE)
    expect_false(inherits(with_line$theme$axis.line.y, "element_blank"))
})

test_that("connect.lighten sets how light the bands are", {
    band_fill <- function(p) {
        b <- ggplot2::ggplot_build(p)
        unique(b$data[[layer_of(p, "GeomPolygon")]]$fill)
    }
    cols <- c(
        A = "#000000",
        B = "#000000",
        C = "#000000"
    )
    p <- plot_composition_bars(comp_meta(), "ct", "geno", connect = "curve", cols = cols)
    expect_equal(band_fill(p), lighten_colour("#000000", 0.8))
    q <- plot_composition_bars(
        comp_meta(),
        "ct",
        "geno",
        connect = "curve",
        cols = cols,
        connect.lighten = 0.5
    )
    expect_equal(band_fill(q), "#808080")
    expect_error(
        plot_composition_bars(comp_meta(), "ct", "geno", connect = "curve", connect.lighten = 1.5),
        "connect.lighten"
    )
})
