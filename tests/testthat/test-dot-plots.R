# Dot-grid engine (G0), DotPlot_bioToolkit() (G1) and assay handling

grid_fixture <- function() {
    data.frame(
        x = factor(c("a", "b", "a", "b"), levels = c("a", "b")),
        y = factor(c("g1", "g1", "g2", "g2"), levels = c("g1", "g2")),
        size = c(50, 80, NA, 20),
        fill = c(1, 2, 3, 4)
    )
}

draw_grid <- function(df, ...) {
    dot_grid_plot(df, x.labels = c(a = "A", b = "B"), y.labels = c(g1 = "Cd3e", g2 = "Ms4a1"),
                  size.limits = c(0, 100), size.breaks = c(50, 100), size.title = "s",
                  fill.colours = c("white", "red"), fill.limits = c(0, 4), fill.title = "f", ...)
}

test_that("dot_grid_plot lists rows top to bottom and keeps rows without dots", {
    df <- grid_fixture()
    df$size[df$y == "g2"] <- NA
    p <- draw_grid(df, y.italic = c(g1 = TRUE, g2 = TRUE))
    b <- ggplot2::ggplot_build(p)
    y <- b$layout$panel_params[[1]]$y
    # The first level is drawn at the top, i.e. last on the scale; g2 has no
    # dot but keeps its row
    expect_equal(y$get_limits(), c("g2", "g1"))
    # Named italic flags give italic gene names (regression: names were used
    # as a lookup by label)
    expect_true(all(vapply(y$get_labels(), function(l) is.call(l) && identical(l[[1]], as.name("italic")),
                           logical(1))))
    expect_equal(nrow(b$data[[2]]), 2)
})

test_that("dot_grid_plot panels are cells x categories, with strip and brackets", {
    df <- grid_fixture()
    df$x_group <- factor(c("G1", "G2", "G1", "G2"))
    df$strip <- c("red", "red", "blue", "blue")
    p <- draw_grid(df, strip.axis = "y", cell.width = 0.2, cell.height = 0.25)
    w <- as.numeric(p$theme$panel.widths)
    h <- as.numeric(p$theme$panel.heights)
    # Strip facet (0.1 in, plus the axis gap of 0.04 in scaled to it) then
    # one column per group (one cell each, plus the axis gap);
    # rows: 2 cells plus the axis gap (bottom); the group bars are outside
    expect_equal(w, c(0.1 * (1 + 0.04 / 0.2), 0.2 + 0.04, 0.2 + 0.04))
    expect_equal(h, 2 * 0.25 + 0.04)
    bars <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), p$layers)
    expect_length(bars, 1)
    expect_equal(bars[[1]]$geom_params$side, "top")
    expect_equal(as.character(bars[[1]]$data$x_group), c("G1", "G2"))
    # Dot markers: a column one cell wide, holding one dot per row
    pd <- draw_grid(df, strip.axis = "y", strip.style = "dot", cell.width = 0.2, cell.height = 0.25)
    expect_equal(as.numeric(pd$theme$panel.widths)[1], 0.2 * (1 + 0.04 / 0.2))
    marker <- Filter(function(l) inherits(l$geom, "GeomPoint") && identical(l$aes_params$shape, 16), pd$layers)
    expect_length(marker, 1)
    expect_setequal(marker[[1]]$aes_params$colour, c("red", "blue"))
})

test_that("dot_grid_plot draws its axes and grid lines on the data panels only", {
    df <- grid_fixture()
    df$x_group <- factor(c("G1", "G2", "G1", "G2"))
    df$strip <- c("red", "red", "blue", "blue")
    p <- draw_grid(df, strip.axis = "y")
    # The theme axes only place the labels
    expect_s3_class(p$theme$axis.line, "ggplot2::element_blank")
    expect_s3_class(p$theme$axis.ticks, "ggplot2::element_blank")
    axes <- Filter(function(l) inherits(l$geom, "GeomGridAxis"), p$layers)
    expect_equal(unname(vapply(axes, function(l) l$geom_params$side, character(1))), c("left", "bottom"))
    # y-axis: one tick per row, on the first data column (not the strip
    # facet); x-axis: one tick per column, each in its group's panel
    expect_equal(nrow(axes[[1]]$data), 2)
    expect_equal(as.character(unique(axes[[1]]$data$x_group)), "G1")
    expect_setequal(as.character(axes[[2]]$data$x_group), c("G1", "G2"))
    # Grid lines (grey segments) are never in the strip facet
    grid_layers <- Filter(function(l) inherits(l$geom, "GeomSegment") && identical(l$aes_params$colour, "grey92"),
                          p$layers)
    expect_length(grid_layers, 2)
    for (l in grid_layers) expect_false(".strip" %in% as.character(l$data$x_group))
    # It renders
    expect_no_error(ggplot2::ggplotGrob(p))
})

test_that("cluster_order keeps groups in order", {
    m <- rbind(a = c(0, 0), b = c(5, 5), c = c(0.1, 0), d = c(9, 9), e = c(9, 8.9), f = c(0, 0.2))
    out <- cluster_order(m, factor(c("x", "x", "x", "y", "y", "y"), levels = c("x", "y")))
    expect_setequal(out[1:3], c("a", "b", "c"))
    expect_setequal(out[4:6], c("d", "e", "f"))
    # Similar rows are neighbours
    expect_equal(abs(match("a", out) - match("c", out)), 1)
})

test_that("group_expression_stats matches a hand computation", {
    expr <- cbind(g = c(0, log1p(1), log1p(3), 0), h = c(1, 1, 1, 1))
    unit <- factor(c("A", "A", "B", "B"))
    lin <- group_expression_stats(expr, unit, "linear", 0, log.data = TRUE)
    expect_equal(lin$mean["A", "g"], log1p(mean(c(0, 1))))
    expect_equal(lin$mean["B", "g"], log1p(mean(c(3, 0))))
    expect_equal(lin$pct[, "g"], c(A = 50, B = 50))
    lg <- group_expression_stats(expr, unit, "log", 0, log.data = TRUE)
    expect_equal(lg$mean["A", "g"], mean(c(0, log1p(1))))
    # Other layers are averaged as they are
    cnt <- group_expression_stats(expr, unit, "linear", 0, log.data = FALSE)
    expect_equal(cnt$mean["B", "g"], mean(c(log1p(3), 0)))
    expect_equal(lin$n, c(A = 2L, B = 2L))
})

test_that("scale_columns z-scores and zeroes constant genes", {
    z <- scale_columns(cbind(a = c(1, 2, 3), b = c(5, 5, 5)))
    expect_equal(as.vector(z[, "a"]), c(-1, 0, 1))
    expect_equal(as.vector(z[, "b"]), c(0, 0, 0))
})

test_that("dot_fill_limits: z-scores centred on 0, each end at most 2.5 or its data", {
    # Each end trimmed to the data on its side (2026-10-14)
    expect_equal(dot_fill_limits(c(-1, 0.5), TRUE, NA, NA), c(-1, 0.5))
    expect_equal(dot_fill_limits(c(-4, 3), TRUE, NA, NA), c(-2.5, 2.5))
    expect_equal(dot_fill_limits(c(-0.3, 4), TRUE, -2, 2), c(-0.3, 2))
    expect_equal(dot_fill_limits(c(0.2, 3), FALSE, NA, NA), c(0, 3))
    expect_equal(dot_fill_limits(c(0.2, 3), FALSE, NA, 2), c(0, 2))
})

test_that("DotPlot_bioToolkit computes % expressing and means per group", {
    so <- sim_seurat()
    expect_message(p <- DotPlot_bioToolkit(so, features = c("Cd3e", "Ms4a1"), group.by = "ct"),
                   "assay 'RNA', layer 'data' \\(the default assay; others: ALT")
    src <- attr(p, "source_data")
    x <- SeuratObject::LayerData(so, assay = "RNA", layer = "data")["Cd3e", ]
    t_cells <- so$ct == "T"
    row <- src[src$group == "T" & src$feature == "Cd3e", ]
    expect_equal(row$pct_expressing, mean(x[t_cells] > 0) * 100)
    expect_equal(row$mean_expression, log1p(mean(expm1(x[t_cells]))))
    expect_equal(row$n_cells, 40)
    # z-scores across the three groups
    m <- src$mean_expression[src$feature == "Cd3e"]
    expect_equal(src$mean_expression_scaled[src$feature == "Cd3e"], as.vector(scale(m)))
    expect_equal(unique(src$assay), "RNA")
    expect_equal(unique(src$average), "linear")
})

test_that("DotPlot_bioToolkit reads the assay given, and reports missing genes", {
    so <- sim_seurat()
    p <- suppressMessages(DotPlot_bioToolkit(so, features = "Cd3e", group.by = "ct", assay = "ALT"))
    src <- attr(p, "source_data")
    x <- SeuratObject::LayerData(so, assay = "ALT", layer = "data")["Cd3e", ]
    expect_equal(src$mean_expression[src$group == "B"], log1p(mean(expm1(x[so$ct == "B"]))))
    expect_equal(unique(src$assay), "ALT")
    # Actb is only in RNA; CD3E only differs in case
    expect_warning(suppressMessages(DotPlot_bioToolkit(so, features = c("Cd3e", "Actb", "CD3E"),
                                                      group.by = "ct", assay = "ALT")),
                   "Actb \\(in assay RNA\\).*CD3E \\(did you mean Cd3e\\?\\)")
    expect_error(DotPlot_bioToolkit(so, features = "Cd3e", group.by = "ct", assay = "nope"), "not found")
})

test_that("DotPlot_bioToolkit groups genes, splits groups and keeps rows below pct.min", {
    so <- sim_seurat()
    p <- suppressMessages(DotPlot_bioToolkit(so, features = list(Tc = c("Cd3e", "Actb"), Bc = "Ms4a1"),
                                             group.by = "cond", split.by = "ct", pct.min = 99))
    src <- attr(p, "source_data")
    expect_equal(levels(src$feature_group), c("Tc", "Bc"))
    expect_equal(nrow(src), 3 * 2 * 3)
    expect_true(any(!src$drawn))
    b <- ggplot2::ggplot_build(p)
    # Six rows (2 conditions in each of 3 cell types) are kept, in 3 row
    # panels (blocks = split.by), named by condition (group.by) and marked
    # by the condition colours
    expect_equal(length(unique(b$layout$layout$ROW)), 3)
    y <- b$layout$panel_params[[1]]$y
    expect_setequal(unlist(lapply(y$get_labels(), as.character)), c("KO", "WT"))
    marker <- Filter(function(l) inherits(l$geom, "GeomPoint") && identical(l$aes_params$shape, 16), p$layers)
    expect_length(unique(marker[[1]]$aes_params$colour), 2)
    # Scaling needs two groups
    expect_message(DotPlot_bioToolkit(so, features = "Cd3e", group.by = "ct", cells = colnames(so)[1:40]),
                   "Only one group")
})

test_that("DotPlot_bioToolkit bars: gene groups on top, split.by on the right, or by the labels", {
    so <- sim_seurat()
    feats <- list(T = "Cd3e", B = c("Ms4a1", "Cd68"))
    p <- suppressMessages(DotPlot_bioToolkit(so, feats, group.by = "ct", split.by = "cond"))
    bars <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), p$layers)
    expect_equal(unname(vapply(bars, function(l) l$geom_params$side, "")), c("top", "right"))
    p2 <- suppressMessages(DotPlot_bioToolkit(so, feats, group.by = "ct", split.by = "cond", flip = TRUE,
                                              feature.group.side = "labels", split.by.cols = "okabe-ito"))
    bars2 <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), p2$layers)
    # Flipped: split.by blocks are columns (top), gene groups rows (labels: left)
    expect_equal(unname(vapply(bars2, function(l) l$geom_params$side, "")), c("top", "left"))
    expect_false(any(bars2[[1]]$data$fill == "black"))
    expect_no_error(ggplot2::ggplotGrob(p2))
    # The % expressing legend is black
    expect_equal(p$scales$get_scales("size")$guide$params$override.aes$fill, "black")
})

test_that("resolve_assay fails clearly on an assay without the layer", {
    so <- sim_seurat()
    so[["RAW"]] <- SeuratObject::CreateAssay5Object(
        counts = SeuratObject::LayerData(so, assay = "RNA", layer = "counts"))
    expect_error(resolve_assay(so, "RAW", "data"), "no 'data' layer.*NormalizeData")
    expect_silent(a <- resolve_assay(so, "ALT", "data", quiet = TRUE))
    expect_equal(a$key, "alt_")
    # A gene of another assay only is not read from it silently
    expect_error(fetch_features(so, "Actb", colnames(so), "data", a), "not in assay 'ALT'.*Actb \\(in RNA")
})

test_that("top_markers keeps each gene with its strongest group", {
    markers <- data.frame(
        gene = c("Cd3e", "Cd3d", "Ms4a1", "Cd79a", "Cd3e", "Lyz2"),
        cluster = factor(c("T", "T", "B", "B", "B", "B"), levels = c("T", "B")),
        avg_log2FC = c(3, 2.5, 4, 3.5, 0.5, -2),
        p_val = 1e-10, p_val_adj = c(1e-8, 1e-8, 1e-8, 1e-8, 1e-8, 1e-8)
    )
    out <- top_markers(markers, n = 2)
    expect_equal(out, list(T = c("Cd3e", "Cd3d"), B = c("Ms4a1", "Cd79a")))
    expect_equal(top_markers(markers, n = 3, unique = FALSE)$B, c("Ms4a1", "Cd79a", "Cd3e"))
})
