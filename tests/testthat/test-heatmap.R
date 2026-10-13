# Marker heatmap (S2)

test_that("heatmap cells: z-scores from all cells, then a reproducible downsample", {
    so <- sim_seurat()
    p <- suppressMessages(DoHeatmap_bioToolkit(so, c("Cd3e", "Ms4a1"), group.by = "ct", max.cells = 10, seed = 3))
    src <- attr(p, "source_data")
    expect_equal(length(unique(src$column)), 30)
    expect_true(all(table(src$group[src$feature == "Cd3e"]) == 10))
    x <- SeuratObject::LayerData(so, assay = "RNA", layer = "data")["Cd3e", ]
    z <- (x - mean(x)) / stats::sd(x)
    cd3e <- src[src$feature == "Cd3e", ]
    expect_equal(cd3e$value_scaled, unname(z[cd3e$column]))
    p2 <- suppressMessages(DoHeatmap_bioToolkit(so, c("Cd3e", "Ms4a1"), group.by = "ct", max.cells = 10, seed = 3))
    expect_identical(attr(p2, "source_data")$column, src$column)
})

test_that("heatmap groups: means in blocks of split.by, groups compared within", {
    so <- sim_seurat()
    p <- suppressMessages(DoHeatmap_bioToolkit(so, list(Tc = "Cd3e", Bc = "Ms4a1"), group.by = "cond",
                                               split.by = "ct", columns = "groups"))
    src <- attr(p, "source_data")
    # 3 blocks (cell types) x 2 groups (conditions) per gene
    expect_equal(nrow(src), 6 * 2)
    expect_equal(levels(droplevels(src$block)), c("T", "B", "Mac"))
    x <- SeuratObject::LayerData(so, assay = "RNA", layer = "data")["Cd3e", ]
    cells <- so$ct == "T" & so$cond == "WT"
    row <- src[src$feature == "Cd3e" & src$block == "T" & src$group == "WT", ]
    expect_equal(row$value, log1p(mean(expm1(x[cells]))))
    expect_equal(row$n_cells, sum(cells))
    m <- src$value[src$feature == "Cd3e"]
    expect_equal(src$value_scaled[src$feature == "Cd3e"], as.vector(scale(m)))
    expect_equal(unique(src$average), "linear")
    # Annotations must be constant within a column of group means
    expect_error(suppressMessages(DoHeatmap_bioToolkit(so, "Cd3e", group.by = "ct", columns = "groups",
                                                       annotations = "cond")), "several values")
})

test_that("heatmap layout: column widths, gaps, row heights and the label rule", {
    cols <- data.frame(block = factor(rep(c("A", "B"), c(6, 4))),
                       group = factor(c(rep(c("g1", "g2"), each = 3), rep("g1", 4))))
    cl <- heatmap_column_layout(cols, "cells", panel.width = 1, cell.width = 0.2)
    # Block A has a gap between its two groups
    expect_true(all(diff(cl$columns$x[cl$columns$block == "A"]) %in% c(1, 1 + (max(cl$columns$x[1:6]) - 6))))
    expect_equal(sum(cl$widths), 1, tolerance = 1e-8)
    gl <- heatmap_column_layout(cols[c(1, 4, 7), ], "groups", panel.width = 1, cell.width = 0.2)
    expect_equal(unname(gl$widths), c(0.4, 0.2))
    ft <- dot_feature_table(list(a = c("x", "y"), b = "z"))
    rl <- heatmap_row_layout(ft, cell.height = 0.13, panel.height = NULL)
    expect_equal(rl$row.height, 0.13)
    expect_equal(rl$rows$y, c(3, 2, 1))
    many <- dot_feature_table(paste0("g", 1:100))
    expect_equal(heatmap_row_layout(many, 0.13, NULL)$row.height, 0.06)
    so <- sim_seurat()
    expect_message(DoHeatmap_bioToolkit(so, c("Cd3e", "Ms4a1"), group.by = "ct", panel.height = 0.1),
                   "too thin to name")
})

test_that("heatmap colours: teal-yellow for z-scores, magma for positive values", {
    expect_equal(heatmap_fill_limits(c(-4, 3), TRUE, TRUE, NA, NA), c(-2.5, 2.5))
    expect_equal(heatmap_fill_limits(c(-0.4, 0.2), FALSE, TRUE, NA, NA), c(-0.4, 0.4))
    expect_equal(heatmap_fill_limits(c(0.5, 3), FALSE, FALSE, NA, NA), c(0, 3))
    pal <- bt_palette("teal-yellow")
    expect_equal(pal[c(1, 6, 11)], c("#00CED1", "#000000", "#FFE800"))
})

test_that("spread_positions keeps labels a line apart, in order", {
    out <- spread_positions(c(5, 5.2, 5.1, 1), gap = 1, lim = c(0.5, 10.5))
    expect_true(all(diff(sort(out)) >= 1 - 1e-9))
    expect_equal(order(out), order(c(5, 5.2, 5.1, 1)))
})

test_that("genomic regions are not gene names", {
    expect_equal(is_genomic_region(c("chr1-100-200", "chrX:5-10", "Nkx2-5", "Cd3e")), c(TRUE, TRUE, FALSE, FALSE))
})

test_that("heatmap palette: dark centre for cells, light centre for group means", {
    so <- sim_seurat()
    mid <- function(p) toupper(p$scales$get_scales("fill")$palette(0.5))
    pc <- suppressMessages(DoHeatmap_bioToolkit(so, c("Cd3e", "Ms4a1"), group.by = "ct"))
    pg <- suppressMessages(DoHeatmap_bioToolkit(so, c("Cd3e", "Ms4a1"), group.by = "ct", columns = "groups"))
    expect_equal(mid(pc), "#000000")
    expect_equal(mid(pg), toupper(bt_palette("blue-magenta")[6]))
})
