# Shared DE-table reader (roadmap item F5)

# One small fixture per method, with the columns each tool writes
de_fixtures <- function() {
    genes <- c("Cd3e", "Ms4a1", "Lyz2", "Actb")
    list(
        edgeR = data.frame(logFC = c(2, -1, 0.1, 0), logCPM = c(5, 6, 7, 8), F = 1:4,
                           PValue = c(1e-6, 1e-3, 0.5, 0.9), FDR = c(4e-6, 2e-3, 0.6, 0.9),
                           row.names = genes),
        DESeq2 = data.frame(baseMean = c(100, 50, 10, 1000), log2FoldChange = c(2, -1, 0.1, 0),
                            lfcSE = 0.2, stat = 1, pvalue = c(0, 1e-3, 0.5, NA),
                            padj = c(0, 2e-3, 0.6, NA), row.names = genes),
        limma = data.frame(logFC = c(2, -1, 0.1, 0), AveExpr = c(5, 6, 7, 8), t = 1,
                           P.Value = c(1e-6, 1e-3, 0.5, 0.9), adj.P.Val = c(4e-6, 2e-3, 0.6, 0.9),
                           B = 0, row.names = genes),
        seurat = data.frame(p_val = c(1e-6, 1e-3, 0.5, 0.9), avg_log2FC = c(2, -1, 0.1, 0),
                            pct.1 = 0.5, pct.2 = 0.2, p_val_adj = c(4e-6, 2e-3, 0.6, 0.9),
                            row.names = genes)
    )
}

test_that("auto picks the right method for each tool", {
    fx <- de_fixtures()
    for (m in names(fx)) {
        res <- suppressMessages(standardize_de(fx[[m]]))
        expect_equal(attr(res, "method"), m, info = m)
        expect_equal(names(res)[1:4], c("feature", "logfc", "pvalue", "padj"), info = m)
        expect_equal(res$feature, rownames(fx[[m]]), info = m)
        expect_equal(res$logfc, c(2, -1, 0.1, 0), info = m)
    }
    expect_message(standardize_de(fx$edgeR), "Detected edgeR")
})

test_that("mean expression and its scale follow the method", {
    fx <- de_fixtures()
    edger <- suppressMessages(standardize_de(fx$edgeR))
    expect_equal(edger$mean_expr, c(5, 6, 7, 8))
    expect_equal(attr(edger, "mean.scale"), "log2")
    deseq <- suppressMessages(standardize_de(fx$DESeq2))
    expect_equal(attr(deseq, "mean.scale"), "linear")
    seurat <- suppressMessages(standardize_de(fx$seurat))
    expect_false("mean_expr" %in% names(seurat))
    # Other columns are kept after the standard ones
    expect_true(all(c("pct.1", "pct.2") %in% names(seurat)))
})

test_that("explicit methods and columns override detection", {
    fx <- de_fixtures()
    d <- fx$DESeq2
    d$shrunk <- d$log2FoldChange / 2
    res <- suppressMessages(standardize_de(d, method = "DESeq2", logfc.col = "shrunk"))
    expect_equal(res$logfc, c(1, -0.5, 0.05, 0))
    custom <- data.frame(gene = c("a", "b"), lfc = c(1, -1), p = c(0.01, 0.2))
    res <- suppressMessages(standardize_de(custom, feature.col = "gene", logfc.col = "lfc", p.col = "p"))
    expect_equal(res$feature, c("a", "b"))
    # No adjusted p-values: computed with BH
    expect_equal(res$padj, stats::p.adjust(c(0.01, 0.2), "BH"))
    expect_error(standardize_de(custom), "Could not detect")
    expect_error(standardize_de(fx$edgeR, logfc.col = "nope"), "not found")
})

test_that("zero p-values are floored and NA p-values kept", {
    fx <- de_fixtures()
    expect_message(res <- standardize_de(fx$DESeq2), "replaced by")
    # A tenth of the smallest non-zero value
    expect_equal(res$pvalue[1], 1e-4)
    expect_equal(res$padj[1], 2e-4)
    expect_true(is.na(res$padj[4]))
    expect_message(standardize_de(fx$DESeq2), "no adjusted p-value")
})

test_that("a list of tables becomes comparisons, with optional groups", {
    fx <- de_fixtures()
    res <- suppressMessages(standardize_de(list("B: KO vs WT" = fx$edgeR, "T: KO vs WT" = fx$seurat),
                                           groups = c("B: KO vs WT" = "B cells", "T: KO vs WT" = "T cells")))
    expect_equal(levels(res$comparison), c("B: KO vs WT", "T: KO vs WT"))
    expect_equal(as.vector(table(res$comparison)), c(4, 4))
    expect_equal(levels(res$group), c("B cells", "T cells"))
    expect_equal(unname(attr(res, "method")), c("edgeR", "seurat"))
    # Seurat has no mean expression, so the scale is edgeR's
    expect_equal(attr(res, "mean.scale"), "log2")
    expect_error(standardize_de(list(fx$edgeR, fx$seurat)), "unique names")
    expect_error(suppressMessages(standardize_de(list(a = fx$edgeR), groups = c(b = "x"))), "No group")
})

test_that("a long table with a comparison column is read, using Seurat's gene column", {
    fx <- de_fixtures()
    long <- rbind(fx$seurat, fx$seurat)
    long$gene <- rep(rownames(fx$seurat), 2)
    long$cluster <- rep(c("0", "1"), each = 4)
    rownames(long) <- paste0(long$gene, rep(c("", "1"), each = 4))
    res <- suppressMessages(standardize_de(long, comparison.col = "cluster"))
    expect_equal(res$feature, long$gene)
    expect_equal(levels(res$comparison), c("0", "1"))
    expect_false("cluster" %in% names(res))
})

test_that("old Seurat natural-log fold changes are refused", {
    old <- data.frame(p_val = 0.01, avg_logFC = 1, p_val_adj = 0.05, row.names = "a")
    expect_error(standardize_de(old), "natural-log")
})
