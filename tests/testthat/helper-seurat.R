# A small simulated Seurat object for the dot plot and assay tests: 120 cells
# of three cell types (40 each), two conditions, five genes with known
# expression, a normalized "RNA" assay and a second assay "ALT" (three genes,
# doubled counts) so that assay selection can be checked.
sim_seurat <- function() {
    testthat::skip_if_not_installed("Seurat")
    testthat::skip_if_not_installed("SeuratObject")
    set.seed(1)
    genes <- c("Cd3e", "Ms4a1", "Cd68", "Actb", "Neg1")
    ct <- rep(c("T", "B", "Mac"), each = 40)
    counts <- matrix(stats::rpois(5 * 120, 1), 5, 120, dimnames = list(genes, paste0("c", 1:120)))
    counts["Cd3e", ct == "T"] <- counts["Cd3e", ct == "T"] + stats::rpois(40, 6)
    counts["Ms4a1", ct == "B"] <- counts["Ms4a1", ct == "B"] + stats::rpois(40, 6)
    counts["Cd68", ct == "Mac"] <- counts["Cd68", ct == "Mac"] + stats::rpois(40, 6)
    counts["Neg1", ] <- 0
    so <- suppressWarnings(SeuratObject::CreateSeuratObject(methods::as(counts, "dgCMatrix")))
    so <- Seurat::NormalizeData(so, verbose = FALSE)
    so$ct <- factor(ct, levels = c("T", "B", "Mac"))
    so$cond <- rep(c("WT", "KO"), 60)
    so[["ALT"]] <- SeuratObject::CreateAssay5Object(counts = methods::as(counts[1:3, ] * 2, "dgCMatrix"))
    so <- Seurat::NormalizeData(so, assay = "ALT", verbose = FALSE)
    so
}
