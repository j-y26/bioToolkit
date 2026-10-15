# Violin and ridge plots (S1)

test_that("feature_flags takes TRUE / FALSE or feature names", {
    f <- c("a", "b", "c")
    expect_equal(feature_flags(TRUE, f, "log"), c(a = TRUE, b = TRUE, c = TRUE))
    expect_equal(feature_flags("b", f, "log"), c(a = FALSE, b = TRUE, c = FALSE))
    expect_warning(feature_flags(c("b", "z"), f, "log"), "not plotted: z")
    expect_error(feature_flags(1, f, "log"), "TRUE, FALSE or feature names")
})

test_that("pct_label marks shares above 0 that round to 0", {
    expect_equal(pct_label(c(0, 0.002, 0.26)), c("0%", "<1%", "26%"))
})

test_that("sample_per_unit caps each unit and is reproducible", {
    cells <- data.frame(cell = paste0("c", 1:30), unit = factor(rep(c("a", "b"), c(25, 5))))
    s1 <- sample_per_unit(cells, 10, seed = 1)
    expect_equal(sum(s1 %in% cells$cell[cells$unit == "a"]), 10)
    expect_true(all(cells$cell[cells$unit == "b"] %in% s1))
    expect_identical(s1, sample_per_unit(cells, 10, seed = 1))
})

test_that("VlnPlot_bioToolkit mixes genes and metadata, reading genes from one assay", {
    so <- sim_seurat()
    so$score <- seq_len(ncol(so))
    expect_message(p <- VlnPlot_bioToolkit(so, c("Cd3e", "score"), group.by = "ct"), "assay 'RNA'")
    expect_s3_class(p, "patchwork")
    src <- export_source_data(p)
    expect_length(src, 2)
    cd3e <- src[[1]]
    x <- SeuratObject::LayerData(so, assay = "RNA", layer = "data")["Cd3e", ]
    expect_equal(cd3e$value, unname(x[cd3e$cell]))
    expect_equal(unique(cd3e$assay), "RNA")
    # Metadata only: no assay is read, nor reported
    expect_no_message(q <- VlnPlot_bioToolkit(so, "score", group.by = "ct"))
    expect_null(attr(q, "source_data")$assay)
})

test_that("VlnPlot_bioToolkit: zeros of genes, log axes and y.max", {
    so <- sim_seurat()
    p <- suppressMessages(VlnPlot_bioToolkit(so, "Cd3e", group.by = "ct", zeros = "annotate"))
    src <- attr(p, "source_data")
    expect_equal(src$drawn, src$value > 0)
    expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomText"), logical(1))))
    so$n <- c(0, seq_len(ncol(so) - 1))
    expect_message(p2 <- VlnPlot_bioToolkit(so, "n", group.by = "ct", log = TRUE), "1 cell\\(s\\) with 'n' <= 0")
    expect_false(attr(p2, "source_data")$drawn[1])
    expect_equal(value_limit("q50", "n", 1:9), 5)
    expect_equal(value_limit(c(n = 3), "n", 1:9), 3)
    expect_true(is.na(value_limit(c(m = 3), "n", 1:9)))
})

test_that("VlnPlot_bioToolkit sizes violins and draws split.by as blocks", {
    so <- sim_seurat()
    p <- suppressMessages(VlnPlot_bioToolkit(so, "Cd3e", group.by = "ct", cell.width = 0.3, panel.height = 2))
    expect_equal(as.numeric(p$theme$panel.widths), 3 * 0.3)
    expect_equal(as.numeric(p$theme$panel.heights), 2)
    # Horizontal: groups on y, first at the top
    pf <- suppressMessages(VlnPlot_bioToolkit(so, "Cd3e", group.by = "ct", flip = TRUE, cell.height = 0.2))
    expect_equal(as.numeric(pf$theme$panel.heights), 3 * 0.2)
    b <- ggplot2::ggplot_build(pf)
    expect_equal(b$layout$panel_params[[1]]$y$get_limits(), c("Mac", "B", "T"))
    # Blocks: one panel per condition (split.by), each with the cell types
    ps <- suppressMessages(VlnPlot_bioToolkit(so, "Cd3e", group.by = "ct", split.by = "cond", cell.width = 0.3))
    expect_equal(nrow(ggplot2::ggplot_build(ps)$layout$layout), 2)
    expect_equal(as.numeric(ps$theme$panel.widths), c(0.9, 0.9))
    # Stacked: one plot per feature
    pst <- suppressMessages(VlnPlot_bioToolkit(so, c("Cd3e", "Ms4a1", "Cd68"), group.by = "ct", stack = TRUE))
    expect_length(pst$patches$plots, 2)
})

test_that("VlnPlot_bioToolkit marks split.by blocks with bars and cuts violins at y.max", {
    so <- sim_seurat()
    ps <- suppressMessages(VlnPlot_bioToolkit(so, "Cd3e", group.by = "ct", split.by = "cond",
                                              split.by.cols = c(WT = "red", KO = "blue")))
    bars <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), ps$layers)
    expect_length(bars, 1)
    expect_equal(bars[[1]]$geom_params$side, "top")
    expect_setequal(bars[[1]]$data$fill, c("red", "blue"))
    pf <- suppressMessages(VlnPlot_bioToolkit(so, "Cd3e", group.by = "ct", split.by = "cond", flip = TRUE,
                                              split.by.side = "labels"))
    expect_equal(Filter(function(l) inherits(l$geom, "GeomGroupBar"), pf$layers)[[1]]$geom_params$side, "left")
    expect_no_error(ggplot2::ggplotGrob(pf))
    # Violins and boxes end at y.max; points above it are not drawn
    top <- 2
    py <- suppressMessages(VlnPlot_bioToolkit(so, "Cd3e", group.by = "ct", y.max = top))
    b <- ggplot2::ggplot_build(py)
    is_violin <- vapply(py$layers, function(l) inherits(l$geom, "GeomViolinTop"), logical(1))
    is_box <- vapply(py$layers, function(l) inherits(l$geom, "GeomBoxplotTop"), logical(1))
    expect_lte(max(b$data[[which(is_violin)]]$y), top)
    expect_lte(max(b$data[[which(is_box)]]$upper), top)
    expect_true(all(attr(py, "source_data")$value[0] <= top))
})

test_that("cut_density closes a violin at the top", {
    d <- data.frame(y = c(0, 1, 2, 3), violinwidth = c(0, 1, 0.5, 0))
    out <- cut_density(d, 1.5)
    expect_equal(out$y, c(0, 1, 1.5))
    expect_equal(out$violinwidth[3], 0.75)
    expect_equal(nrow(cut_density(d, -1)), 0)
    expect_equal(cut_density(d, 5), d)
})

test_that("VlnPlot_bioToolkit and RidgePlot_bioToolkit take a data frame", {
    meta <- data.frame(n = c(1:20, 101:120), mito = runif(40, 0, 0.1), s = rep(c("a", "b"), each = 20),
                       row.names = paste0("c", 1:40))
    p <- VlnPlot_bioToolkit(meta, c("n", "mito"), group.by = "s", log = "n", percent = "mito")
    expect_no_error(ggplot2::ggplotGrob(p[[1]]))
    expect_equal(p[[1]]$scales$get_scales("y")$get_transformation()$name, "log-10")
    expect_equal(sort(unique(as.character(attr(p[[1]], "source_data")$group))), c("a", "b"))
    skip_if_not_installed("ggridges")
    r <- RidgePlot_bioToolkit(meta, "n", group.by = "s", cell.height = 0.3)
    expect_no_error(ggplot2::ggplotGrob(r))
    # Rows numbered from the bottom: "a" (first) at the top
    expect_equal(r$scales$get_scales("y")$labels, c("b", "a"))
    expect_error(VlnPlot_bioToolkit(meta, "n", group.by = "zz"), "not a column")
})
