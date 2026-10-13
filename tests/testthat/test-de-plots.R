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
    # Faded and ringed layers, explained by one titled legend that spells out
    # the cutoffs
    expect_equal(sum(vapply(p$layers, function(l) inherits(l$geom, "GeomPoint"), logical(1))), 2)
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
    expect_equal(levels(p$data$x_group), c("T", "B"))
    expect_equal(levels(p$data$y_group), c("g1", "g2"))
})

test_that("plot_multi_volcano counts, band and offsets", {
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
    # Band spans exactly +-logfc.cutoff
    band <- b$data[[which(vapply(p$layers, function(l) inherits(l$geom, "GeomRect"), logical(1)))]]
    expect_equal(unique(band$ymax), 1)
    expect_equal(unique(band$ymin), -1)
    # Points stay within their column and are reproducible
    col <- as.integer(src$comparison)
    expect_true(all(abs(src$x - col) <= 0.4 + 1e-9))
    p2 <- suppressMessages(plot_multi_volcano(de, logfc.cutoff = 1, seed = 1))
    expect_equal(attr(p2, "source_data")$x, src$x)
    # Only DE genes unless show.ns
    expect_true(all(src$direction != "ns"))
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
