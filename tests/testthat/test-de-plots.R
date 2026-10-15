# plot_de_dots() (G2) and plot_multi_volcano() (G4)

multi_de_fixture <- function() {
    set.seed(2)
    make <- function(shift) {
        d <- data.frame(logFC = c(4, -3, 0.1, 2, stats::rnorm(96, shift)), logCPM = 5,
                        PValue = c(1e-12, 1e-9, 0.8, 1e-4, stats::runif(96)^3),
                        row.names = c("Cd3e", "Ms4a1", "Actb", "Cd68", paste0("G", 1:96)))
        d$FDR <- stats::p.adjust(d$PValue, "BH")
        d
    }
    list(A = make(0), B = make(0.5), C = make(-0.5))
}

test_that("plot_de_dots picks top genes and rings significant results", {
    de <- multi_de_fixture()
    p <- suppressMessages(plot_de_dots(de, n = 2))
    src <- attr(p, "source_data")
    expect_true(all(c("Cd3e", "Ms4a1") %in% src$feature))
    expect_equal(nrow(src), length(unique(src$feature)) * 3)
    expect_equal(src$significant, !is.na(src$padj) & src$padj < 0.05)
    # Comparison markers (dots, by default), faded and ringed layers; the
    # significance explained by one titled legend that spells out the cutoffs
    expect_equal(sum(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))), 3)
    p_none <- suppressMessages(plot_de_dots(de, n = 2, comparison.marker = "none"))
    expect_equal(sum(vapply(p_none$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))), 2)
    expect_equal(p$scales$get_scales("alpha")$name, "Significance")
    expect_equal(p$scales$get_scales("colour")$name, "Significance")
    expect_equal(significance_labels(0.05, 1), c("padj < 0.05\nand |log2FC| > 1", "Not significant"))
})

test_that("plot_de_dots: grey and hide modes, squished sizes and case warnings", {
    de <- multi_de_fixture()
    p <- suppressMessages(plot_de_dots(de, features = c("Cd3e", "Actb"), significance = "grey", size.max = 5))
    d <- p$data
    expect_true(all(is.na(d$fill[d$x == "Actb"])))
    lab <- ggplot2::ggplot_build(p)$plot$scales$get_scales("size")$get_labels()
    expect_match(utils::tail(lab, 1), "^≥5$")
    p2 <- suppressMessages(plot_de_dots(de, features = c("Cd3e", "Actb"), significance = "hide"))
    expect_true(all(is.na(p2$data$size[p2$data$x == "Actb"])))
    expect_warning(suppressMessages(plot_de_dots(de, features = c("CD3E", "Actb"))), "matching only in case: Cd3e")
})

test_that("plot_de_dots draws gene groups and comparison groups", {
    de <- multi_de_fixture()
    p <- suppressMessages(plot_de_dots(de, features = list(T = "Cd3e", B = c("Ms4a1", "Cd68")),
                                       groups = c(A = "g1", B = "g1", C = "g2")))
    # The comparison markers are a column facet of their own (".strip")
    expect_equal(setdiff(levels(p$data$x_group), ".strip"), c("T", "B"))
    expect_equal(levels(p$data$y_group), c("g1", "g2"))
    # Gene groups and comparison groups are marked by bars, black by default
    bars <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), p$layers)
    expect_length(bars, 2)
    expect_equal(bars[[1]]$geom_params$side, "top")
    expect_equal(bars[[2]]$geom_params$side, "right")
    expect_true(all(bars[[1]]$data$fill == "black"))
    p2 <- suppressMessages(plot_de_dots(de, features = list(T = "Cd3e", B = c("Ms4a1", "Cd68")),
                                        groups = c(A = "g1", B = "g1", C = "g2"),
                                        groups.cols = c(g1 = "red", g2 = "blue"), groups.side = "labels"))
    bars2 <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), p2$layers)
    expect_equal(bars2[[2]]$geom_params$side, "left")
    expect_equal(bars2[[2]]$data$fill, c("red", "blue"))
    expect_no_error(ggplot2::ggplotGrob(p2))
})

test_that("plot_de_dots colour limits reach each side's own extreme", {
    de <- multi_de_fixture()
    p <- suppressMessages(plot_de_dots(de, features = c("Cd3e", "Ms4a1")))
    lim <- p$scales$get_scales("fill")$limits
    sig <- attr(p, "source_data")
    sig <- sig[sig$significant, ]
    expect_equal(lim, range(sig$logfc))
})

test_that("plot_multi_volcano counts, strips, cutoff lines and offsets", {
    de <- multi_de_fixture()
    p <- suppressMessages(plot_multi_volcano(de, logfc.cutoff = 1, seed = 1))
    src <- attr(p, "source_data")
    std <- suppressMessages(standardize_de(de))
    expected_up <- tapply(std$padj < 0.05 & std$logfc > 1, std$comparison, sum)
    b <- ggplot2::ggplot_build(p)
    counts <- b$data[[which(vapply(p$layers, function(l) {
        inherits(l$geom, "GeomText") && identical(l$aes_params$colour, "firebrick4")
    }, logical(1)))[1]]]
    expect_equal(as.numeric(counts$label), as.vector(expected_up))
    # No band: dashed lines at +-logfc.cutoff, and a strip per comparison
    expect_false(any(vapply(p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1))))
    cut <- Filter(function(l) inherits(l$geom, "GeomHline"), p$layers)[[1]]
    expect_equal(sort(cut$data$yintercept), c(-1, 1))
    strips <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), p$layers)[[1]]
    expect_equal(strips$geom_params$side, "top")
    expect_equal(strips$data$fill, unname(bt_palette("scanpy-20", 3)))
    # Comparison names on top
    expect_equal(b$layout$panel_params[[1]]$x$get_labels(), c("A", "B", "C"))
    # Points stay within their column and are reproducible
    col <- as.integer(src$comparison)
    expect_true(all(abs(src$x - col) <= 0.4 + 1e-9))
    p2 <- suppressMessages(plot_multi_volcano(de, logfc.cutoff = 1, seed = 1))
    expect_equal(attr(p2, "source_data")$x, src$x)
    # Non-significant genes are drawn (in grey) by default, not with show.ns = FALSE
    expect_true(any(src$direction == "ns"))
    p3 <- suppressMessages(plot_multi_volcano(de, logfc.cutoff = 1, show.ns = FALSE))
    expect_true(all(attr(p3, "source_data")$direction != "ns"))
    # band.cols is the deprecated name of comparison.cols
    expect_warning(suppressMessages(plot_multi_volcano(de, band.cols = c(A = "red", B = "blue", C = "green"))),
                   "deprecated")
})

test_that("plot_multi_volcano colour limits are trimmed to each direction", {
    de <- multi_de_fixture()
    p <- suppressMessages(plot_multi_volcano(de, logfc.cutoff = 1))
    src <- attr(p, "source_data")
    lim <- p$scales$get_scales("colour")$limits
    up <- -log10(src$padj[src$direction == "up"])
    down <- -log10(src$padj[src$direction == "down"])
    expect_lte(lim[2], max(up) + 1e-9)
    expect_gte(lim[1], -max(down) - 1e-9)
})

test_that("plot_multi_volcano flips comparisons into rows", {
    de <- multi_de_fixture()
    p <- suppressMessages(plot_multi_volcano(de, logfc.cutoff = 1, flip = TRUE, label = "none"))
    b <- ggplot2::ggplot_build(p)
    # First comparison at the top; log2 fold change on x
    y <- b$layout$panel_params[[1]]$y
    expect_equal(y$get_labels()[order(-y$get_breaks())], c("A", "B", "C"))
    strips <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), p$layers)[[1]]
    expect_equal(strips$geom_params$side, "left")
    expect_true(inherits(Filter(function(l) inherits(l$geom, "GeomVline"), p$layers)[[1]]$geom, "GeomVline"))
    expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("plot_multi_volcano labels per column and orders comparisons", {
    de <- multi_de_fixture()
    p <- suppressMessages(plot_multi_volcano(de, label.n = 1, logfc.cutoff = 1))
    src <- attr(p, "source_data")
    per <- tapply(src$labelled, list(src$comparison, src$direction), sum)
    expect_true(all(per[, c("up", "down")] <= 1))
    p2 <- suppressMessages(plot_multi_volcano(de, order = "n_de", logfc.cutoff = 1))
    n_de <- table(attr(p2, "source_data")$comparison)
    expect_equal(levels(attr(p2, "source_data")$comparison), names(sort(n_de, decreasing = TRUE)))
    expect_error(suppressMessages(plot_multi_volcano(de$A)), "several comparisons")
})

test_that("plot_multi_volcano names the comparisons of grouped panels", {
    de <- multi_de_fixture()
    names(de) <- c("First long comparison", "Second long comparison", "Third long comparison")
    groups <- c("First long comparison" = "G2", "Second long comparison" = "G1",
                "Third long comparison" = "G2")
    p <- suppressMessages(plot_multi_volcano(de, groups = groups, logfc.cutoff = 1, label = "none"))
    b <- ggplot2::ggplot_build(p)
    # Each group's panel names its own comparisons (comparisons of a group
    # together; groups in their order of appearance)
    labs <- lapply(b$layout$panel_params, function(pp) pp$x$get_labels()[!is.na(pp$x$get_breaks())])
    expect_equal(labs[[1]], c("First long comparison", "Third long comparison"))
    expect_equal(labs[[2]], "Second long comparison")
    # Groups: black bars at the bottom (opposite the names) by default
    bars <- Filter(function(l) inherits(l$geom, "GeomGroupBar"), p$layers)
    expect_equal(bars[[2]]$geom_params$side, "bottom")
    expect_no_warning(ggplot2::ggplotGrob(p))
})

test_that("multi_volcano_offsets: sina width follows the density", {
    df <- data.frame(comparison = "A", direction = "up", logfc = c(rep(1, 50), 5),
                     neg_log10_padj = 1)
    df$logfc[1:50] <- df$logfc[1:50] + seq(0, 0.2, length.out = 50)
    off <- with_seed(1, multi_volcano_offsets(df, "sina", 0.4))
    expect_true(all(abs(off) <= 0.4))
    # The isolated gene at logfc 5 sits near the centre
    expect_lt(abs(off[51]), 0.05)
    pv <- multi_volcano_offsets(transform(df, neg_log10_padj = c(rep(1, 50), 4)), "pvalue", 0.4)
    expect_equal(pv[51], 0.4)
})

test_that("multi_volcano_offsets: one shape per side of 0, whatever the significance", {
    set.seed(2)
    lfc <- c(stats::rnorm(400, sd = 0.6), stats::rnorm(40, 2.5, 0.5), -stats::rexp(40) - 1)
    df <- data.frame(comparison = "A", logfc = lfc, neg_log10_padj = stats::runif(length(lfc), 0, 5))
    df$direction <- ifelse(abs(df$logfc) > 1 & df$neg_log10_padj > 1.3,
                           ifelse(df$logfc > 0, "up", "down"), "ns")
    # The direction (the cutoffs) does not change where points go
    off <- with_seed(1, multi_volcano_offsets(df, "sina", 0.4))
    off_all_ns <- with_seed(1, multi_volcano_offsets(transform(df, direction = "ns"), "sina", 0.4))
    expect_equal(off, off_all_ns)
    expect_true(all(abs(off) <= 0.4))
    # Each side fills the column at its widest, near 0 (no pinch at 0: the
    # density is reflected at the boundary)
    near0 <- abs(df$logfc) < 0.2
    expect_gt(max(abs(off[near0 & df$logfc > 0])), 0.3)
    expect_gt(max(abs(off[near0 & df$logfc < 0])), 0.3)
    # pvalue: grey genes are placed on the scale of their side
    pv <- multi_volcano_offsets(df, "pvalue", 0.4)
    up_side <- df$logfc >= 0
    expect_equal(pv[up_side], -0.4 + 0.8 * df$neg_log10_padj[up_side] / max(df$neg_log10_padj[up_side]))
})
