# DE table helpers (rank_de_genes, generate_de_rank_file, filter_de_genes) and
# the deprecated-argument mapping used by them and the Seurat plot formatters
# (renamed to dot.case in 0.10.0)

de_edger <- data.frame(
    logFC = c(2.1, -1.4, 0.2, 1.5, NA), PValue = c(1e-8, 1e-4, 0.6, 0.01, NA),
    FDR = c(3e-8, 2e-4, 0.6, 0.04, NA), row.names = c("ISG15", "CD3E", "ACTB", "MX1", "XX")
)

test_that("filter_de_genes splits by direction and drops NA p-values", {
    res <- suppressMessages(filter_de_genes(de_edger, logfc.cutoff = 1))
    expect_equal(rownames(res$up_reg), c("ISG15", "MX1"))
    expect_equal(rownames(res$down_reg), "CD3E")
})

test_that("old argument names still work, with a deprecation warning", {
    new <- suppressMessages(filter_de_genes(de_edger, logfc.cutoff = 1))
    expect_warning(
        old <- suppressMessages(filter_de_genes(de_results = de_edger, log2fc_min = 1, de_method = "edgeR")),
        class = "deprecatedWarning"
    )
    expect_identical(old, new)
    # positional calls written for the old signature give the same result
    expect_identical(filter_de_genes(de_edger, 1, 0.05, "edgeR", "FDR", "logFC"), new)
    expect_error(filter_de_genes(de_edger, padj.col = "FDR", pvalue_col = "FDR"), "not both")
    expect_error(filter_de_genes(de_edger, foo = 1), "unused argument")
})

test_that("columns are found from any supported tool, case-insensitively", {
    dds <- data.frame(log2FoldChange = c(1, -2), pvalue = c(.01, .001), padj = c(.02, .002),
                      row.names = c("A", "B"))
    expect_equal(rank_de_genes(dds, method = "deseq2")$rank, c(1, 2))
    # a trimmed edgeR table (no PValue) is still recognised
    expect_message(rank_de_genes(de_edger[1:4, c("logFC", "FDR")]), "edgeR")
    expect_error(filter_de_genes(data.frame(x = 1)), "Could not detect")
})

test_that("generate_de_rank_file writes a two-column rank file", {
    f <- tempfile(fileext = ".rnk")
    rnk <- suppressMessages(generate_de_rank_file(de_edger[1:4, ], file = f))
    expect_equal(names(rnk), c("GeneName", "rank"))
    expect_equal(rnk$GeneName[1], "ISG15")
    expect_equal(nrow(utils::read.delim(f)), 4)
})

test_that("format_sc_plots accepts old and new names and keeps patchwork widths", {
    df <- data.frame(UMAP_1 = 1:10, UMAP_2 = 10:1)
    p <- ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2)) + ggplot2::geom_point()
    pw <- patchwork::wrap_plots(p, p, ncol = 2, widths = c(1, 2))
    expect_identical(format_sc_plots(pw)$patches$layout$widths, c(1, 2))
    expect_warning(format_sc_plots(list(p, p), ncol = 2, add_arrow = FALSE), class = "deprecatedWarning")
    expect_warning(add_dimension_arrows(p, arrow_color = "red"), class = "deprecatedWarning")
})

test_that("embedding plots take the dot.case arrow arguments in `...`", {
    new <- resolve_arrow_args(list(arrow.length = 1, text.font.size = 8), "DimPlot_bioToolkit")
    expect_identical(new, list(arrow_length = 1, text_font_size = 8))
    expect_warning(old <- resolve_arrow_args(list(arrow_length = 1), "DimPlot_bioToolkit"),
                   class = "deprecatedWarning")
    expect_identical(old, list(arrow_length = 1))
    expect_error(resolve_arrow_args(list(foo = 1), "DimPlot_bioToolkit"), "unused argument")
})
