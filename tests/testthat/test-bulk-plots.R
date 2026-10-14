# Bulk plots (roadmap items B1 MA plot, B2 PCA, B4 variance explained,
# B5 quadrant plot)

# DE tables with known classes ------------------------------------------------

ma_fixture <- function() {
    data.frame(
        logFC = c(3, 2, -2.5, 0.5, -0.3, 1),
        logCPM = c(2, 8, 5, 3, 10, 1),
        PValue = c(1e-12, 1e-8, 1e-10, 0.2, 0.8, 0.01),
        FDR = c(1e-10, 1e-6, 1e-8, 0.4, 0.9, NA),
        row.names = c("Up1", "Up2", "Down1", "NS1", "NS2", "Filtered")
    )
}

quiet <- function(expr) suppressMessages(suppressWarnings(expr))

test_that("de_direction() classes features by both cutoffs, NA padj as not significant", {
    d <- de_direction(c(2, -2, 0.5, 2, NA), c(0.01, 0.01, 0.01, NA, 0.01), 0.05, 1)
    expect_equal(d, c("up", "down", "ns", "ns", "ns"))
})

test_that("plot_ma() draws log2 means as they are and linear means on a log10 axis", {
    p <- quiet(plot_ma(ma_fixture(), label = "none"))
    d <- p$data
    expect_false("Filtered" %in% d$feature)
    expect_equal(d$mean_plot, d$mean_expr)
    expect_equal(p$scales$get_scales("x")$get_transformation()$name, "identity")
    dir <- setNames(as.character(d$direction), d$feature)
    expect_equal(unname(dir[c("Up1", "Up2", "Down1", "NS1")]),
                 c("Up-regulated", "Up-regulated", "Down-regulated", "Not significant"))
    # Not significant first
    expect_true(max(which(d$direction == "Not significant")) < min(which(d$direction != "Not significant")))

    ds <- data.frame(log2FoldChange = c(2, -2, 0.1, 1), baseMean = c(100, 5000, 0, 20),
                     pvalue = c(1e-6, 1e-6, 0.5, 0.5), padj = c(1e-5, 1e-5, 0.6, 0.6),
                     row.names = c("A", "B", "Zero", "C"))
    p2 <- quiet(plot_ma(ds, label = "none"))
    expect_equal(p2$scales$get_scales("x")$get_transformation()$name, "log-10")
    # A mean of 0 cannot be drawn on a log axis
    expect_false("Zero" %in% p2$data$feature)
})

test_that("plot_ma() needs a mean expression and its scale", {
    seu <- data.frame(avg_log2FC = c(1, -1), p_val = c(0.01, 0.01), p_val_adj = c(0.02, 0.02),
                      row.names = c("A", "B"))
    expect_error(quiet(plot_ma(seu)), "no mean expression")
    custom <- data.frame(fc = c(1, -1, 0), pv = c(0.01, 0.01, 0.5), m = c(1, 2, 3), row.names = c("A", "B", "C"))
    expect_error(quiet(plot_ma(custom, logfc.col = "fc", padj.col = "pv", mean.col = "m")), "set `mean.scale`")
    p <- quiet(plot_ma(custom, logfc.col = "fc", padj.col = "pv", mean.col = "m", mean.scale = "log2",
                       label = "none"))
    expect_equal(nrow(p$data), 3)
})

test_that("plot_ma() labels the top features per direction and draws a trend", {
    d <- quiet(plot_ma(ma_fixture(), label.n = 1))$data
    expect_setequal(d$feature[d$labelled], c("Up1", "Down1"))
    p <- quiet(plot_ma(ma_fixture(), label = "none", trend = TRUE))
    expect_true(any(vapply(p$layers, function(l) inherits(l$geom, "GeomLine"), logical(1))))
})


# PCA ---------------------------------------------------------------------------

pca_fixture <- function() {
    set.seed(1)
    m <- matrix(rnorm(300 * 8, mean = 6), nrow = 300,
                dimnames = list(paste0("Gene", 1:300), paste0("S", 1:8)))
    m[1:60, 5:8] <- m[1:60, 5:8] + 3
    meta <- data.frame(treatment = rep(c("Control", "Treated"), each = 4),
                       batch = rep(c("A", "B"), 4),
                       depth = seq(1, 8),
                       row.names = colnames(m))
    list(m = m, meta = meta)
}

test_that("plot_pca() uses the most variable features and reports % variance in the titles", {
    f <- pca_fixture()
    p <- quiet(plot_pca(f$m, f$meta, color.by = "treatment", n.top = 100))
    top <- f$m[order(-apply(f$m, 1, var))[1:100], ]
    pr <- prcomp(t(top))
    pct <- pr$sdev^2 / sum(pr$sdev^2) * 100
    expect_equal(p$labels$x, sprintf("PC1 (%.1f%%)", pct[1]))
    expect_equal(p$labels$y, sprintf("PC2 (%.1f%%)", pct[2]))
    src <- attr(p, "source_data")
    expect_equal(abs(src$PC1), abs(unname(pr$x[, 1])))
    expect_true(all(c("sample", "treatment", "batch") %in% names(src)))
    # The treatment separates on PC1
    expect_true(all(sign(src$PC1[1:4]) != sign(src$PC1[5:8][1])))
})

test_that("plot_pca() matches metadata by sample name and checks inputs", {
    f <- pca_fixture()
    shuffled <- f$meta[rev(rownames(f$meta)), ]
    p <- quiet(plot_pca(f$m, shuffled, color.by = "treatment"))
    src <- attr(p, "source_data")
    expect_equal(src$treatment, f$meta[src$sample, "treatment"])
    expect_error(quiet(plot_pca(f$m, f$meta[1:6, ], color.by = "treatment")), "not in the row names")
    expect_error(quiet(plot_pca(f$m, f$meta, color.by = "nope")), "column of `metadata`")
    many <- f$meta
    many$id <- rownames(many)
    expect_error(quiet(plot_pca(f$m, many, shape.by = "id")), "at most 5 shapes")
    expect_error(quiet(plot_pca(f$m, dims = c(1, 1))), "two different")
})

test_that("plot_pca() takes prcomp results, flags raw counts, and adds a scree panel", {
    f <- pca_fixture()
    pr <- prcomp(t(f$m))
    p <- plot_pca(pr, f$meta, color.by = "depth")
    expect_true(is.numeric(attr(p, "source_data")$depth))
    counts <- round(2^f$m * 100)
    expect_warning(suppressMessages(plot_pca(counts)), "raw counts")
    s <- quiet(plot_pca(f$m, f$meta, color.by = "treatment", scree = TRUE, scree.n = 5))
    expect_s3_class(s, "patchwork")
    tabs <- export_source_data(s)
    expect_length(tabs, 2)
    expect_equal(nrow(tabs[[2]]), 5)
})

test_that("plot_pca() draws ellipses only for groups of 4 or more samples", {
    f <- pca_fixture()
    meta <- f$meta
    meta$group <- c(rep("big", 5), rep("small", 3))
    expect_message(plot_pca(f$m, meta, color.by = "group", ellipse = TRUE), "fewer than 4 samples: small")
})


# Variance explained ------------------------------------------------------------

test_that("variance_fractions() recovers known fractions in a balanced design", {
    # Two orthogonal factors: their sums of squares add up exactly
    meta <- data.frame(a = rep(c("x", "y"), each = 4), b = rep(c("u", "v"), 4))
    a <- ifelse(meta$a == "x", -1, 1)
    b <- ifelse(meta$b == "u", -1, 1)
    e <- c(0.3, -0.3, -0.3, 0.3, -0.3, 0.3, 0.3, -0.3)   # orthogonal to a and b
    m <- rbind(g1 = 2 * a + 1 * b + e, g2 = 0.5 * a + e)
    fr <- variance_fractions(m, meta, c("a", "b"))
    tss <- rowSums((m - rowMeans(m))^2)
    expect_equal(fr["g1", "a"], 4 * 8 / tss[["g1"]])
    expect_equal(fr["g1", "b"], 1 * 8 / tss[["g1"]])
    expect_equal(fr["g2", "b"], 0)
    # Orthogonal covariates: nothing shared, fractions sum to 1
    expect_equal(unname(rowSums(fr)), c(1, 1))
})

test_that("variance_fractions() splits R^2 by LMG: sums to 1, order-free, shared variance split", {
    set.seed(3)
    # Unbalanced design: depth partly follows condition
    meta <- data.frame(cond = rep(c("a", "b"), c(5, 7)), depth = c(rnorm(5, 1), rnorm(7, 2)))
    m <- matrix(rnorm(20 * 12), 20, 12, dimnames = list(paste0("g", 1:20), NULL))
    m <- m + outer(rnorm(20), meta$depth) + outer(rnorm(20), meta$cond == "b")
    fr <- variance_fractions(m, meta, c("cond", "depth"))
    expect_equal(unname(rowSums(fr)), rep(1, 20))
    expect_true(all(as.matrix(fr) >= -1e-12))
    fr2 <- variance_fractions(m, meta, c("depth", "cond"))
    expect_equal(fr2$cond, fr$cond)
    # Two covariates: share = mean of the R^2 gained first and second
    r2 <- function(f) summary(lm(f))$r.squared
    y <- m[1, ]
    a <- 0.5 * (r2(y ~ meta$cond) + r2(y ~ meta$cond + meta$depth) - r2(y ~ meta$depth))
    expect_equal(fr$cond[1], a)
})

test_that("variance_fractions() refuses confounded or saturating covariates", {
    meta <- data.frame(a = rep(c("x", "y"), each = 3), b = rep(c("p", "q"), each = 3), id = letters[1:6])
    m <- matrix(rnorm(12), 2, 6)
    expect_error(variance_fractions(m, meta, c("a", "b")), "confounded")
    expect_error(variance_fractions(m, meta, "id"), "different level for every sample")
})

test_that("plot_variance_explained() orders covariates by median, Residuals last", {
    f <- pca_fixture()
    p <- quiet(plot_variance_explained(f$m, f$meta, covariates = c("batch", "treatment")))
    lv <- levels(p$data$group)
    expect_equal(lv[length(lv)], "Residuals")
    expect_equal(lv[1], "treatment")
    src <- attr(p, "source_data")
    expect_equal(names(src), c("feature", "batch", "treatment", "Residuals"))
})

test_that("plot_variance_explained() reads variancePartition-like results and draws bars", {
    vp <- data.frame(Individual = c(0.5, 0.1), Tissue = c(0.2, 0.6), Residuals = c(0.3, 0.3),
                     row.names = c("Gene1", "Gene2"))
    p <- quiet(plot_variance_explained(vp))
    expect_equal(levels(p$data$group), c("Tissue", "Individual", "Residuals"))
    b <- quiet(plot_variance_explained(vp, features = c("Gene2", "Gene1")))
    expect_equal(levels(b$data$feature), c("Gene1", "Gene2"))
    expect_warning(suppressMessages(plot_variance_explained(vp, features = c("gene1", "Gene2"))),
                   "matching only in case: Gene1")
})

test_that("the pcs view gives the R^2 of each PC on each covariate", {
    f <- pca_fixture()
    p <- quiet(plot_variance_explained(f$m, f$meta, view = "pcs", n.pcs = 3, n.top = 300))
    expect_s3_class(p, "patchwork")
    tiles <- export_source_data(p)[[2]]
    pr <- prcomp(t(f$m))
    r2 <- summary(lm(pr$x[, 1] ~ f$meta$treatment))$r.squared
    expect_equal(tiles$r2[tiles$pc == "PC1" & tiles$covariate == "treatment"], r2)
    expect_equal(nrow(tiles), 3 * 3)
})


# Quadrant ----------------------------------------------------------------------

quadrant_fixture <- function() {
    genes <- c("CC1", "CC2", "DC1", "XO", "YO", "NS", "OnlyX")
    x <- data.frame(logFC = c(3, -2, 2.5, 2, 0.2, 0.1, 1), PValue = 1e-3,
                    FDR = c(1e-6, 1e-4, 1e-5, 1e-3, 0.5, 0.9, 1e-3), logCPM = 5, row.names = genes)
    y <- data.frame(logFC = c(2, -1.5, -2, 0.3, 1.8, -0.1), PValue = 1e-3,
                    FDR = c(1e-5, 1e-3, 1e-4, 0.4, 1e-4, 0.8), logCPM = 5, row.names = genes[1:6])
    list(x = x, y = y)
}

test_that("plot_fc_quadrant() classes shared features and counts quadrants", {
    f <- quadrant_fixture()
    p <- quiet(plot_fc_quadrant(f$x, f$y, name.x = "A", name.y = "B", label = "none"))
    d <- attr(p, "source_data")
    expect_false("OnlyX" %in% d$feature)
    cls <- setNames(as.character(d$category), d$feature)
    expect_equal(unname(cls[c("CC1", "CC2", "DC1", "XO", "YO", "NS")]),
                 c("Both, same direction", "Both, same direction", "Both, opposite", "A only", "B only", "Neither"))
    # Drawn from neither to concordant
    expect_equal(as.character(d$category[1]), "Neither")
    expect_equal(quadrant_counts(c(1, -1, -1, 1, 0), c(1, 1, -1, -1, 1)), c(1, 1, 1, 1))
    # Classes not asked for are grey
    p2 <- quiet(plot_fc_quadrant(f$x, f$y, categories = "both", label = "none"))
    d2 <- attr(p2, "source_data")
    expect_equal(as.character(d2$category[d2$feature == "XO"]), "Neither")
})

test_that("plot_fc_quadrant() shows the correlation and labels the top per quadrant", {
    f <- quadrant_fixture()
    p <- quiet(plot_fc_quadrant(f$x, f$y, label.n = 1))
    d <- attr(p, "source_data")
    expect_setequal(d$feature[d$labelled], c("CC1", "CC2", "DC1"))
    r <- cor(d$logfc_x, d$logfc_y, method = "spearman")
    expect_match(paste(deparse(p$labels$subtitle), collapse = ""), formatC(r, format = "f", digits = 2), fixed = TRUE)
    p2 <- quiet(plot_fc_quadrant(f$x, f$y, label = "selected", label.features = "NS", correlation = "none"))
    expect_null(p2$labels$subtitle)
    expect_equal(attr(p2, "source_data")$feature[attr(p2, "source_data")$labelled], "NS")
})

test_that("map_features() keeps one feature per gene", {
    d <- data.frame(feature = c("p1", "p2", "p3", "p4"), logfc = c(1, -3, 2, 1), padj = c(1e-5, 0.2, 1e-3, NA))
    map <- c(p1 = "G1", p2 = "G1", p3 = "G2", p4 = NA)
    m <- quiet(map_features(d, map, "min_padj"))
    expect_equal(m$feature_y[m$feature == "G1"], "p1")
    expect_equal(m$n_mapped[m$feature == "G1"], 2L)
    expect_false("p4" %in% m$feature_y)
    m2 <- quiet(map_features(d, data.frame(peak = names(map), gene = map), "max_abs"))
    expect_equal(m2$feature_y[m2$feature == "G1"], "p2")
})

test_that("plot_fc_quadrant() reads peaks through map.y and rejects duplicates", {
    f <- quadrant_fixture()
    peaks <- f$y
    rownames(peaks) <- paste0("peak", seq_len(nrow(peaks)))
    map <- setNames(rownames(f$y), rownames(peaks))
    p <- quiet(plot_fc_quadrant(f$x, peaks, map.y = map, label = "none"))
    expect_equal(nrow(attr(p, "source_data")), 6)
    dup <- standardize_de(f$y)
    dup$feature[2] <- dup$feature[1]
    expect_error(quiet(plot_fc_quadrant(f$x, dup)), "more than once")
})

test_that("quadrant_colours() accepts partial names and five colours", {
    expect_equal(quadrant_colours(c(ns = "grey90"))[["ns"]], "grey90")
    expect_equal(unname(quadrant_colours(c("red", "blue", "green", "orange", "grey"))),
                 c("red", "blue", "green", "orange", "grey"))
    expect_error(quadrant_colours(c("red", "blue")), "five colours")
    expect_error(quadrant_colours(c(both = "red")), "named by")
})
