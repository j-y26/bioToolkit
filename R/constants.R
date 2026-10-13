# Define global constants used in the package

# Column names in the DE result data frame from edgeR and DESeq2, used by
# rank_de_genes() and filter_de_genes(). Plots read DE tables with
# standardize_de() and DE_COLUMN_PRESETS below instead
PVAL_COL_EDGER <- "FDR"
FC_COL_EDGER <- "logFC"
PVAL_COL_DESEQ2 <- "padj"
FC_COL_DESEQ2 <- "log2FoldChange"
PVAL_COL_LIMMA <- "adj.P.Val"
FC_COL_LIMMA <- "logFC"
PVAL_COL_SEURAT <- "p_val_adj"
FC_COL_SEURAT <- "avg_log2FC"

# Tabix indexing params
MAX_TABIX_END <- 536870912L

# Columns of DE result tables per method, read by standardize_de() (roadmap
# item F5). `mean_scale` says how the mean expression column is expressed:
# "log2" values are plotted as they are, "linear" ones (DESeq2 baseMean) on a
# log10 axis. Seurat's FindMarkers() reports no mean expression.
DE_COLUMN_PRESETS <- list(
    edgeR = c(logfc = "logFC", pvalue = "PValue", padj = "FDR", mean_expr = "logCPM", mean_scale = "log2"),
    DESeq2 = c(logfc = "log2FoldChange", pvalue = "pvalue", padj = "padj", mean_expr = "baseMean", mean_scale = "linear"),
    limma = c(logfc = "logFC", pvalue = "P.Value", padj = "adj.P.Val", mean_expr = "AveExpr", mean_scale = "log2"),
    seurat = c(logfc = "avg_log2FC", pvalue = "p_val", padj = "p_val_adj", mean_expr = NA, mean_scale = NA)
)
