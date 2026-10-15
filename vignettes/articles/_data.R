# Public example data for the bioToolkit website guides
#
# Each loader downloads a public dataset once, processes it and caches the
# result as an RDS file, so that the guides load small, ready objects instead of
# repeating the processing. The cache lives in BT_DOCS_CACHE (environment
# variable; the GitHub Actions workflow keeps it between runs) or in the user
# cache directory. Delete the files to rebuild them; bump CACHE_VERSION when the
# processing changes so that cached files from an older version are not used.
#
# Datasets
#   Kang et al. 2018, Nat Biotechnol 36:89-94 (doi:10.1038/nbt.4042): PBMCs of
#     8 lupus patients, control or stimulated with interferon-beta for 6 h,
#     10x Genomics; from Bioconductor muscData::Kang18_8vs8(), with the
#     authors' cell type labels.
#   Alasoo et al. 2018, Nat Genet 50:424-431 (doi:10.1038/s41588-018-0046-7):
#     iPSC-derived macrophages of 6 donors, naive, IFN-gamma, Salmonella
#     (SL1344) or both; bulk RNA-seq, gene-level counts from Bioconductor
#     macrophage (tximeta).

CACHE_VERSION <- "v4"

docs_cache_dir <- function() {
    dir <- Sys.getenv("BT_DOCS_CACHE", tools::R_user_dir("bioToolkitDocs", "cache"))
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    dir
}

# Read a cached object, or build and cache it
cached <- function(name, build) {
    file <- file.path(docs_cache_dir(), paste0(name, "_", CACHE_VERSION, ".rds"))
    if (file.exists(file)) return(readRDS(file))
    obj <- build()
    saveRDS(obj, file)
    obj
}

# Cell types in a biological order (lymphoid, then myeloid, then the rest),
# with the authors' labels shortened as in Seurat's tutorials, so that labels
# on embeddings stay short
KANG_CELL_TYPES <- c(
    "CD4 T cells" = "CD4 T",
    "CD8 T cells" = "CD8 T",
    "NK cells" = "NK",
    "B cells" = "B",
    "CD14+ Monocytes" = "CD14 Mono",
    "FCGR3A+ Monocytes" = "FCGR3A Mono",
    "Dendritic cells" = "DC",
    "Megakaryocytes" = "Megakaryocyte"
)

# Singlets of the Kang data as a SingleCellExperiment, with tidy metadata
kang_singlets <- function() {
    sce <- suppressMessages(muscData::Kang18_8vs8())
    # Doublets and ambiguous droplets (demuxlet calls) are removed, and the
    # few cells without a cell type label in the published annotation
    sce <- sce[, sce$multiplets == "singlet" & !is.na(sce$cell)]
    md <- SummarizedExperiment::colData(sce)
    # Donors in numeric order (D101, D107, D1015, ...)
    md$donor <- factor(paste0("D", md$ind), levels = paste0("D", sort(unique(as.integer(md$ind)))))
    md$condition <- factor(ifelse(md$stim == "ctrl", "Control", "IFN-beta"),
                           levels = c("Control", "IFN-beta"))
    md$sample <- factor(paste(md$donor, md$condition, sep = "_"),
                        levels = paste(rep(levels(md$donor), each = 2), levels(md$condition), sep = "_"))
    md$cell_type <- factor(unname(KANG_CELL_TYPES[as.character(md$cell)]), levels = KANG_CELL_TYPES)
    SummarizedExperiment::colData(sce) <- md[c("donor", "condition", "sample", "cell_type")]
    sce
}

#' Kang 2018 PBMCs as a Seurat object (8,000 cells)
#'
#' A random 8,000 of the labelled singlets (seed 1; random, so the cell type proportions
#' of each sample are kept), genes detected in at least 10 cells. Log-normalized
#' (NormalizeData), 2,000 variable genes, 30 PCs. (The published matrix has
#' no mitochondrial genes, so there is no percent.mt QC metric.) Stimulated and control cells
#' were captured in separate 10x lanes, and interferon-beta changes the
#' expression of every cell type, so the two conditions separate on a UMAP of
#' the PCs ("umap.unintegrated"). The default "umap" is computed on Harmony
#' (integrated over condition), as in Seurat's own integration tutorial on this
#' dataset, so that each cell type forms one cluster.
docs_kang <- function() {
    cached("kang_seurat", function() {
        sce <- kang_singlets()
        set.seed(1)
        sce <- sce[, sort(sample(ncol(sce), 8000))]
        counts <- SingleCellExperiment::counts(sce)
        counts <- counts[Matrix::rowSums(counts > 0) >= 10, ]
        md <- as.data.frame(SummarizedExperiment::colData(sce))
        so <- SeuratObject::CreateSeuratObject(counts, meta.data = md, project = "Kang2018")
        so <- Seurat::NormalizeData(so, verbose = FALSE)
        so <- Seurat::FindVariableFeatures(so, nfeatures = 2000, verbose = FALSE)
        so <- Seurat::ScaleData(so, verbose = FALSE)
        # Exact PCA (approx = FALSE): the default approximate (irlba) PCA
        # occasionally failed numerically on this matrix (a first PC with a
        # standard deviation of 1e16), which scrambled the UMAP
        so <- Seurat::RunPCA(so, npcs = 30, approx = FALSE, verbose = FALSE, seed.use = 1)
        sdev <- Seurat::Stdev(so, "pca")
        if (!all(is.finite(sdev)) || sdev[1] > 100) stop("PCA failed (PC1 sd = ", sdev[1], ")")
        so <- Seurat::RunUMAP(so, dims = 1:20, reduction = "pca", seed.use = 1,
                              reduction.name = "umap.unintegrated",
                              reduction.key = "rawUMAP_", verbose = FALSE)
        so <- harmony::RunHarmony(so, group.by.vars = "condition", verbose = FALSE)
        so <- Seurat::RunUMAP(so, dims = 1:20, reduction = "harmony", seed.use = 1,
                              reduction.name = "umap", reduction.key = "UMAP_",
                              verbose = FALSE)
        Seurat::Idents(so) <- "cell_type"
        # scale.data is only needed for the PCA; dropping it keeps the cache small
        so <- Seurat::DietSeurat(so, layers = c("counts", "data"),
                                       dimreducs = c("pca", "harmony", "umap", "umap.unintegrated"))
        so
    })
}

#' Pseudobulk DE of IFN-beta vs control per Kang cell type (edgeR)
#'
#' Counts of all singlets summed per donor, condition and cell type; one test
#' per cell type with at least 10 cells in every sample, paired by donor
#' (design ~ donor + condition), filterByExpr, TMM, quasi-likelihood F test.
#' Samples, not cells, are the replicates, so the p-values are not inflated by
#' pseudo-replication (Squair et al. 2021, Nat Commun 12:5692).
#'
#' @return A named list of edgeR topTags tables (all genes), one per cell type
docs_kang_de <- function() {
    cached("kang_pseudobulk_de", function() {
        sce <- kang_singlets()
        counts <- SingleCellExperiment::counts(sce)
        md <- as.data.frame(SummarizedExperiment::colData(sce))
        out <- list()
        for (ct in levels(md$cell_type)) {
            cells <- md$cell_type == ct
            n <- table(droplevels(md$sample[cells]))
            if (length(n) < nlevels(md$sample) || any(n < 10)) next
            sample <- droplevels(md$sample[cells])
            # Sum counts per sample: counts x (cells -> samples) indicator matrix
            ind <- Matrix::sparse.model.matrix(~ 0 + sample)
            pb <- as.matrix(counts[, cells] %*% ind)
            colnames(pb) <- levels(sample)
            info <- unique(md[cells, c("sample", "donor", "condition")])
            info <- info[match(colnames(pb), info$sample), ]
            y <- edgeR::DGEList(pb)
            design <- stats::model.matrix(~ donor + condition, data = info)
            keep <- edgeR::filterByExpr(y, group = info$condition)
            y <- edgeR::normLibSizes(y[keep, , keep.lib.sizes = FALSE])
            fit <- edgeR::glmQLFit(y, design)
            res <- edgeR::glmQLFTest(fit, coef = "conditionIFN-beta")
            out[[ct]] <- as.data.frame(edgeR::topTags(res, n = Inf))
        }
        out
    })
}

#' MSigDB Hallmark gene sets (human symbols)
#'
#' msigdbr (>= 10) downloads the whole MSigDB archive from zenodo.org on first
#' use; on the GitHub runners this connection sometimes times out, so the
#' download is retried a few times and only the small Hallmark table is cached
#' (later builds do not contact zenodo at all).
#'
#' @return A data frame with `gs_name` and `gene_symbol` (one row per gene and
#'         set, as msigdbr::msigdbr(species = "Homo sapiens", collection = "H"))
docs_hallmark <- function() {
    cached("hallmark", function() {
        tries <- 3
        for (i in seq_len(tries)) {
            h <- tryCatch(
                msigdbr::msigdbr(species = "Homo sapiens", collection = "H"),
                error = function(e) e
            )
            if (!inherits(h, "error")) break
            message("msigdbr download failed (try ", i, " of ", tries, "): ", conditionMessage(h))
            # Give up after the last try; otherwise wait before the next one
            if (i == tries) stop(h)
            Sys.sleep(30)
        }
        as.data.frame(h[c("gs_name", "gene_symbol")])
    })
}

#' Alasoo 2018 macrophages: counts and sample information
#'
#' Gene-level counts (tximeta estimates, rounded) with gene symbols as row
#' names (Ensembl genes without a symbol dropped; for duplicated symbols, the
#' gene with the most counts kept), and the sample table.
#'
#' @return A list with `counts` (genes x 24 samples) and `samples` (data
#'         frame: donor, condition, harvest and RNA extraction dates)
docs_macrophage <- function() {
    cached("macrophage", function() {
        e <- new.env()
        utils::data("gse", package = "macrophage", envir = e)
        gse <- e$gse
        counts <- round(SummarizedExperiment::assay(gse, "counts"))
        symbol <- SummarizedExperiment::rowData(gse)$SYMBOL
        keep <- !is.na(symbol) & nzchar(symbol)
        counts <- counts[keep, ]
        symbol <- symbol[keep]
        # One row per symbol: the gene with the highest total count
        o <- order(-rowSums(counts))
        counts <- counts[o, ]
        symbol <- symbol[o]
        counts <- counts[!duplicated(symbol), ]
        rownames(counts) <- symbol[!duplicated(symbol)]
        storage.mode(counts) <- "integer"

        cd <- as.data.frame(SummarizedExperiment::colData(gse))
        samples <- data.frame(
            donor = factor(sub("_.*", "", cd$line_id)),
            condition = factor(cd$condition_name, levels = c("naive", "IFNg", "SL1344", "IFNg_SL1344")),
            harvest_date = factor(cd$macrophage_harvest),
            extraction_date = factor(cd$rna_extraction),
            row.names = paste(sub("_.*", "", cd$line_id), cd$condition_name, sep = "_")
        )
        colnames(counts) <- rownames(samples)
        list(counts = counts, samples = samples)
    })
}
