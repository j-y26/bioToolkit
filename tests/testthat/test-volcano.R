# Volcano plot (roadmap item R2); de_dot_plot() (R3) was replaced by plot_de_dots() (G2)

# A small edgeR-like table with known significance:
#   Up1..Up3     significant, logFC > 0 (Up1 most significant)
#   Down1, Down2 significant, logFC < 0
#   Weak         significant but |logFC| = 0.5
#   NS1, NS2     not significant
#   Filtered     NA adjusted p-value (DESeq2-style independent filtering)
volcano_fixture <- function() {
    data.frame(
        logFC = c(3, 2, 1.5, -2.5, -1.2, 0.5, 2.5, -0.3, 1),
        logCPM = 5,
        PValue = c(1e-12, 1e-8, 1e-5, 1e-10, 1e-4, 1e-3, 0.2, 0.8, 0.01),
        FDR = c(1e-10, 1e-6, 1e-3, 1e-8, 1e-2, 0.02, 0.4, 0.9, NA),
        row.names = c("Up1", "Up2", "Up3", "Down1", "Down2", "Weak", "NS1", "NS2", "Filtered")
    )
}

quiet_volcano <- function(...) suppressMessages(plot_volcano(...))

test_that("directions follow the adjusted p-value and logFC cutoffs", {
    p <- quiet_volcano(volcano_fixture(), label = "none")
    d <- p$data
    dir <- setNames(as.character(d$direction), d$feature)
    expect_equal(unname(dir[c("Up1", "Up2", "Up3", "Weak")]), rep("Up-regulated", 4))
    expect_equal(unname(dir[c("Down1", "Down2")]), rep("Down-regulated", 2))
    expect_equal(unname(dir[c("NS1", "NS2")]), rep("Not significant", 2))
    # NA adjusted p-value: not drawn
    expect_false("Filtered" %in% d$feature)
    # y is the -log10 adjusted p-value
    expect_equal(d$neg_log10_padj, -log10(d$padj))

    p1 <- quiet_volcano(volcano_fixture(), label = "none", logfc.cutoff = 1)
    expect_equal(as.character(p1$data$direction[p1$data$feature == "Weak"]), "Not significant")
})

test_that("non-significant points are drawn first", {
    d <- quiet_volcano(volcano_fixture(), label = "none")$data
    ns <- which(d$direction == "Not significant")
    expect_true(max(ns) < min(which(d$direction != "Not significant")))
})

test_that("label = 'top' labels the top features per direction, plus selected ones", {
    d <- quiet_volcano(volcano_fixture(), label.n = 1)$data
    expect_setequal(d$feature[d$labelled], c("Up1", "Down1"))
    d <- quiet_volcano(volcano_fixture(), label.n = 2, label.features = "NS2")$data
    expect_setequal(d$feature[d$labelled], c("Up1", "Up2", "Down1", "Down2", "NS2"))
    # Ranked by absolute logFC instead
    d <- quiet_volcano(volcano_fixture(), label.n = 1, label.rank = "logfc")$data
    expect_setequal(d$feature[d$labelled], c("Up1", "Down1"))
    d <- quiet_volcano(volcano_fixture(), label.n = 0)$data
    expect_false(any(d$labelled))
})

test_that("label = 'selected' labels only the selected features, significant or not", {
    d <- quiet_volcano(volcano_fixture(), label = "selected", label.features = c("NS1", "Up3"))$data
    expect_setequal(d$feature[d$labelled], c("NS1", "Up3"))
    expect_error(quiet_volcano(volcano_fixture(), label = "selected"), "needs `label.features`")
})

test_that("label = 'none' draws no labels", {
    p <- quiet_volcano(volcano_fixture(), label = "none", label.features = "Up1")
    expect_false(any(p$data$labelled))
    has_repel <- vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1))
    expect_false(any(has_repel))
})

test_that("missing label features warn, pointing to case mismatches", {
    expect_warning(
        quiet_volcano(volcano_fixture(), label = "selected", label.features = c("up1", "Up2")),
        "matching only in case: Up1"
    )
})

test_that("colours accept partial names, three unnamed colours, and reject bad input", {
    expect_equal(volcano_colours(NULL), c(up = "#C91F28", down = "#174EA6", ns = "grey75"))
    expect_equal(volcano_colours(c(ns = "grey90"))[["ns"]], "grey90")
    expect_equal(unname(volcano_colours(c("red", "blue", "grey"))), c("red", "blue", "grey"))
    expect_error(volcano_colours(c(sig = "red")), "named 'up'")
    expect_error(volcano_colours(c("red", "blue")), "three colours")
})

test_that("y limits leave the requested headroom in inches", {
    lim <- volcano_y_limits(y.max = 10, panel.height = 3, headroom = 0.5)
    per_inch <- diff(lim$limits) / 3
    expect_equal((lim$limits[2] - lim$points.top) / per_inch, 0.5)
    expect_true(lim$arrow.y > lim$points.top && lim$arrow.y < lim$limits[2])
    expect_error(volcano_y_limits(10, panel.height = 0.5, headroom = 0.5), "too small")
})

test_that("the volcano has a fixed panel and keeps vector labels", {
    p <- quiet_volcano(volcano_fixture(), annotate.counts = "arrow")
    # A null device, so that building the grob writes no Rplots.pdf
    grDevices::pdf(NULL)
    on.exit(grDevices::dev.off(), add = TRUE)
    gt <- ggplot2::ggplotGrob(p)
    size <- fixed_gtable_size(gt)
    expect_true(all(is.finite(unlist(size))))
    expect_equal(p$theme$panel.widths, grid::unit(3, "in"))
    expect_error(quiet_volcano(list(a = volcano_fixture(), b = volcano_fixture())), "one comparison")
})

test_that("labels are bold, and bold italic for genes", {
    face <- function(p) {
        l <- p$layers[vapply(p$layers, function(l) inherits(l$geom, "GeomTextRepel"), logical(1))][[1]]
        l$aes_params$fontface
    }
    expect_equal(face(quiet_volcano(volcano_fixture())), "bold.italic")
    expect_equal(face(quiet_volcano(volcano_fixture(), italic.genes = FALSE)), "bold")
})
