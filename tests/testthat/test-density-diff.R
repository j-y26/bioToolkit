# Density difference plot (S4)

# Cells on a 2D embedding: condition A half at (0, 0) and half at (5, 5);
# condition B mostly at (5, 5); samples s1-s4 (A: s1, s2; B: s3, s4)
sim_density_object <- function(n_a = c(300, 300), n_b = c(100, 500), dup_a = 1) {
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if_not_installed("SeuratObject")
    set.seed(2)
    pts <- function(n, cx) cbind(stats::rnorm(n, cx, 0.6), stats::rnorm(n, cx, 0.6))
    a <- rbind(pts(n_a[1], 0), pts(n_a[2], 5))
    if (dup_a > 1) a <- a[rep(seq_len(nrow(a)), dup_a), ]
    b <- rbind(pts(n_b[1], 0), pts(n_b[2], 5))
    emb <- rbind(a, b)
    n <- nrow(emb)
    rownames(emb) <- paste0("c", seq_len(n))
    colnames(emb) <- c("umap_1", "umap_2")
    counts <- matrix(1, 2, n, dimnames = list(c("g1", "g2"), rownames(emb)))
    so <- suppressWarnings(SeuratObject::CreateSeuratObject(methods::as(counts, "dgCMatrix")))
    so[["umap"]] <- SeuratObject::CreateDimReducObject(embeddings = emb, key = "umap_", assay = "RNA")
    so$cond <- rep(c("A", "B"), c(nrow(a), nrow(b)))
    so$sample <- c(rep(c("s1", "s2"), length.out = nrow(a)), rep(c("s3", "s4"), length.out = nrow(b)))
    so
}

test_that("kde_grid smooths binned weights; condition_density sums to 1", {
    rd <- list(xlim = c(-3, 3), ylim = c(-3, 3))
    g <- kde_grid(c(0, 0), c(0, 0), c(1, 1), rd$xlim, rd$ylim, c(0.5, 0.5), n_grid = 64)
    expect_equal(dim(g), c(64, 64))
    # Peak at the cells
    peak <- which(g == max(g), arr.ind = TRUE)[1, ]
    expect_true(all(abs(peak - 32.5) <= 1))
    d <- data.frame(x = stats::rnorm(50), y = stats::rnorm(50))
    dens <- condition_density(d, rd, c(0.5, 0.5), 64, per_sample = FALSE)
    expect_equal(sum(dens), 1)
    expect_equal(attr(dens, "n_cells"), 50)
    # Per sample: the mean of each sample's normalized density
    d$sample <- rep(c("a", "b"), c(10, 40))
    ds <- condition_density(d, rd, c(0.5, 0.5), 64, per_sample = TRUE)
    one <- function(i) {
        g <- kde_grid(d$x[i], d$y[i], rep(1, length(i)), rd$xlim, rd$ylim, c(0.5, 0.5), 64)
        g / sum(g)
    }
    expect_equal(as.vector(ds), as.vector((one(1:10) + one(11:50)) / 2))
    expect_equal(attr(ds, "n_samples"), 2)
})

test_that("DensityDiffPlot_bioToolkit finds the region enriched in the compared condition", {
    so <- sim_density_object()
    expect_message(p <- DensityDiffPlot_bioToolkit(so, split.by = "cond", reference = "A"),
                   "A \\(600 cells\\), B \\(600 cells\\)")
    src <- attr(p, "source_data")
    near <- function(cx) abs(src$x - cx) < 0.5 & abs(src$y - cx) < 0.5
    expect_gt(stats::median(src$value[near(5)], na.rm = TRUE), 0)
    expect_lt(stats::median(src$value[near(0)], na.rm = TRUE), 0)
    expect_equal(levels(src$panel), "B vs A")
    expect_equal(as.numeric(p$theme$panel.heights), as.numeric(p$theme$panel.widths))
    # Raster: one row per grid point, masked points without a value
    r <- suppressMessages(DensityDiffPlot_bioToolkit(so, split.by = "cond", draw = "raster", n.grid = 64))
    rs <- attr(r, "source_data")
    expect_equal(nrow(rs), 64^2)
    expect_true(all(is.na(rs$value[rs$masked])))
    expect_error(DensityDiffPlot_bioToolkit(so, split.by = "cond", reference = "Z"), "Not level")
})

test_that("cell numbers do not change the map; samples weigh the same", {
    base <- sim_density_object()
    doubled <- sim_density_object(dup_a = 2)
    p1 <- suppressMessages(DensityDiffPlot_bioToolkit(base, split.by = "cond", reference = "A",
                                                       mask.quantile = 0))
    p2 <- suppressMessages(DensityDiffPlot_bioToolkit(doubled, split.by = "cond", reference = "A",
                                                       mask.quantile = 0))
    # Every A cell duplicated: the same densities, so the same values at the
    # B cells (whose positions are the same in both objects)
    v1 <- attr(p1, "source_data")
    v2 <- attr(p2, "source_data")
    b1 <- v1[v1$condition == "B", ]
    b2 <- v2[v2$condition == "B", ]
    # The bandwidth is set from all cells, which differ; compare loosely
    expect_gt(stats::cor(b1$value, b2$value), 0.99)
    s <- suppressMessages(DensityDiffPlot_bioToolkit(base, split.by = "cond", sample.by = "sample"))
    expect_true("density_compare" %in% names(attr(s, "source_data")))
})
