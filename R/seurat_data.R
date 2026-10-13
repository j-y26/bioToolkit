# Fetching expression from Seurat objects
#
# Design
# ------
# A Seurat object often holds several assays with the same genes: "RNA"
# (LogNormalize), "SCT" (SCTransform), "integrated" (corrected values for a
# few thousand genes), ATAC gene activities... Seurat functions read the
# *default* assay unless told otherwise, so the same call can plot different
# numbers depending on what `DefaultAssay()` was last set to, and FetchData()
# silently falls back to other assays for genes missing from the default one.
# Every plot of expression in this package therefore:
#   * takes an `assay` argument (NULL = the default assay, as in Seurat), and
#     says which assay and layer it read, every time (`resolve_assay()`);
#   * fails with a clear message when that assay lacks the layer (e.g. a
#     Seurat v5 "RNA" assay with only counts, never normalized);
#   * fetches genes with the assay key prefix (e.g. "sct_Cdh5"), so that the
#     values come from that assay and nowhere else (`fetch_features()`). This
#     also reads Seurat v5 split layers ("data.sample1", "data.sample2", ...)
#     as one, as FetchData() joins them;
#   * records the assay and layer in its source data.
#
# This mirrors SCpubr's check_and_set_assay() (assay = NULL -> default
# assay), but without changing the default assay of the object and with the
# choice reported.


#' Resolve the assay and layer expression is read from
#'
#' @param object A Seurat object
#'
#' @param assay Assay name, or NULL for the default assay
#'
#' @param slot Layer name (e.g. "data", "counts")
#'
#' @param quiet Whether to skip the message naming the assay and layer
#'
#' @return A list with the `assay` name, its `key` (e.g. "sct_"), the
#'         `layer`, its `features` (row names) and whether the assay was the
#'         `default`
#'
#' @keywords internal
resolve_assay <- function(object,
                          assay,
                          slot,
                          quiet = FALSE) {
    assays <- SeuratObject::Assays(object)
    is_default <- is.null(assay)
    if (is_default) assay <- SeuratObject::DefaultAssay(object)
    if (!(is.character(assay) && length(assay) == 1 && assay %in% assays)) {
        stop("Assay '", paste(assay, collapse = " "), "' not found; the object has: ",
             paste(assays, collapse = ", "), call. = FALSE)
    }
    if (!(is.character(slot) && length(slot) == 1 && nzchar(slot))) {
        stop("`slot` must be a single layer name, e.g. \"data\"", call. = FALSE)
    }
    # Layers() warns when nothing matches; the absence is reported below. A
    # Seurat v5 assay may hold one layer per sample ("data.sample1", ...),
    # which match the search too
    available <- SeuratObject::Layers(object[[assay]])
    found <- suppressWarnings(SeuratObject::Layers(object[[assay]], search = slot))
    if (length(found) == 0) {
        hint <- if (slot == "data") " (normalize it with Seurat::NormalizeData(), or choose another assay)" else ""
        stop("Assay '", assay, "' has no '", slot, "' layer", hint, ". Its layers: ",
             paste(utils::head(available, 6), collapse = ", "),
             if (length(available) > 6) ", ...", call. = FALSE)
    }
    if (!quiet) {
        others <- setdiff(assays, assay)
        message("Expression from assay '", assay, "', layer '", slot, "'",
                if (is_default) " (the default assay",
                if (is_default && length(others)) paste0("; others: ", paste(others, collapse = ", "),
                                                         "; set `assay` to choose"),
                if (is_default) ")")
    }
    list(
        assay = assay,
        key = SeuratObject::Key(object[[assay]]),
        layer = slot,
        features = rownames(object[[assay]]),
        default = is_default
    )
}


#' Fetch features for a fixed set of cells
#'
#' Genes of the assay are fetched with its key prefix (e.g. "sct_Cdh5"), so
#' that they come from that assay only. Metadata columns and anything else
#' \code{Seurat::FetchData()} accepts (reduction dimensions, features already
#' prefixed with a key) are fetched as they are; a metadata column takes
#' precedence over a gene of the same name, as in FetchData().
#'
#' @param object A Seurat object
#'
#' @param features Features to fetch
#'
#' @param cells Cell names, in the order they are wanted
#'
#' @param slot Layer to fetch expression from
#'
#' @param assay The result of \code{resolve_assay()}, or NULL to let
#'              FetchData() search the default assay
#'
#' @return A data frame with one column per feature and one row per cell
#'
#' @keywords internal
fetch_features <- function(object,
                           features,
                           cells,
                           slot,
                           assay = NULL) {
    vars <- features
    if (!is.null(assay)) {
        meta <- colnames(object[[]])
        in_assay <- features %in% assay$features & !(features %in% meta)
        # A gene of another assay only: say where it is, rather than letting
        # FetchData() read it from that assay without notice
        elsewhere <- features[!in_assay & !(features %in% meta) & !has_assay_key(object, features)]
        elsewhere <- elsewhere[vapply(elsewhere, function(f) {
            any(vapply(SeuratObject::Assays(object), function(a) f %in% rownames(object[[a]]), logical(1)))
        }, logical(1))]
        if (length(elsewhere)) {
            where <- vapply(elsewhere, function(f) {
                paste(Filter(function(a) f %in% rownames(object[[a]]), SeuratObject::Assays(object)), collapse = ", ")
            }, character(1))
            stop("Feature(s) not in assay '", assay$assay, "': ",
                 paste0(elsewhere, " (in ", where, ")", collapse = ", "),
                 ". Set `assay`, or prefix them with the assay key (e.g. \"",
                 SeuratObject::Key(object[[strsplit(where[1], ", ")[[1]][1]]]), elsewhere[1], "\")",
                 call. = FALSE)
        }
        vars[in_assay] <- paste0(assay$key, features[in_assay])
    }
    data <- withCallingHandlers(
        Seurat::FetchData(object, vars = vars, cells = cells, layer = slot),
        warning = function(w) {
            msg <- conditionMessage(w)
            if (grepl("missing data|were not found", msg)) invokeRestart("muffleWarning")
        }
    )
    missing <- features[!(vars %in% colnames(data))]
    if (length(missing)) {
        stop("Feature(s) not found: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    data <- data[match(cells, rownames(data)), vars, drop = FALSE]
    colnames(data) <- features
    for (f in features) {
        if (!is.numeric(data[[f]])) stop("Feature '", f, "' is not numeric", call. = FALSE)
    }
    data
}


#' Whether feature names start with the key of an assay of the object
#'
#' @param object A Seurat object
#'
#' @param features Feature names
#'
#' @return A logical vector along `features`
#'
#' @keywords internal
has_assay_key <- function(object,
                          features) {
    keys <- vapply(SeuratObject::Assays(object), function(a) SeuratObject::Key(object[[a]]), character(1))
    vapply(features, function(f) any(startsWith(f, keys)), logical(1), USE.NAMES = FALSE)
}


#' Split genes into those of an assay and those missing from it
#'
#' Missing genes are reported in a warning, naming those that match a gene of
#' the assay only in case (gene symbols differ in case between species: human
#' CD3E, mouse and rat Cd3e) and the assays that do hold the others.
#'
#' @param object A Seurat object
#'
#' @param genes Gene names
#'
#' @param assay The result of \code{resolve_assay()}
#'
#' @return The genes found in the assay, in the order given
#'
#' @keywords internal
genes_in_assay <- function(object,
                           genes,
                           assay) {
    found <- genes %in% assay$features
    if (!all(found)) {
        missing <- unique(genes[!found])
        notes <- vapply(missing, function(g) {
            near <- assay$features[tolower(assay$features) == tolower(g)]
            if (length(near)) return(paste0(g, " (did you mean ", near[1], "?)"))
            other <- Filter(function(a) g %in% rownames(object[[a]]), setdiff(SeuratObject::Assays(object), assay$assay))
            if (length(other)) paste0(g, " (in assay ", paste(other, collapse = ", "), ")") else g
        }, character(1))
        warning(length(missing), " feature(s) not in assay '", assay$assay, "' are left out: ",
                paste(notes, collapse = ", "), call. = FALSE)
    }
    genes[found]
}


#' Whether features are genes, to set their names in italics
#'
#' @param object A Seurat object
#'
#' @param features Feature names, possibly prefixed by an assay key
#'
#' @return A logical vector named by feature
#'
#' @keywords internal
feature_is_gene <- function(object,
                            features) {
    assays <- Seurat::Assays(object)
    genes <- unique(unlist(lapply(assays, function(a) rownames(object[[a]]))))
    keys <- vapply(assays, function(a) SeuratObject::Key(object[[a]]), character(1))
    meta <- colnames(object[[]])
    is_gene <- vapply(features, function(f) {
        if (f %in% meta) return(FALSE)
        if (f %in% genes) return(TRUE)
        key <- keys[startsWith(f, keys)]
        length(key) > 0 && any(substring(f, nchar(key) + 1) %in% genes)
    }, logical(1))
    stats::setNames(is_gene, features)
}
