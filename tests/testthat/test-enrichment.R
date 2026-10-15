# standardize_enrichment(), plot_enrichment_grid(), plot_enrichment_dots(),
# plot_enrichment_bars() and plot_gsea_rank() (G3)

fgsea_fixture <- function(seed = 1) {
    set.seed(seed)
    d <- data.frame(pathway = c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA",
                                paste0("HALLMARK_SET_", 1:10)),
                    pval = c(1e-6, 1e-4, stats::runif(10)), NES = c(2.5, -2, stats::rnorm(10)), size = 100)
    d$padj <- stats::p.adjust(d$pval, "BH")
    d$leadingEdge <- I(rep(list(c("Nfkbia", "Jun")), 12))
    d
}

test_that("standardize_enrichment reads fgsea, clusterProfiler and GSEA desktop tables", {
    fg <- standardize_enrichment(fgsea_fixture())
    expect_equal(fg$term[1:2], c("TNFA signaling via NFKB", "Hypoxia"))
    expect_equal(fg$term_id[1], "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
    expect_equal(fg$leading_edge[1], "Nfkbia/Jun")
    expect_equal(unique(fg$type), "gsea")

    cp <- data.frame(ID = "GOBP_DNA_REPAIR", Description = "GOBP_DNA_REPAIR", setSize = 50, NES = 1.8,
                     pvalue = 1e-3, p.adjust = 0.01, core_enrichment = "A/B")
    out <- standardize_enrichment(cp)
    expect_equal(out[, c("term", "nes", "padj", "size")],
                 data.frame(term = "DNA repair", nes = 1.8, padj = 0.01, size = 50))

    ora <- data.frame(ID = "x", Description = "cell migration", GeneRatio = "31/72", BgRatio = "458/2803",
                      pvalue = 1e-6, p.adjust = 1e-4, geneID = "A/B", Count = 31)
    out <- standardize_enrichment(ora)
    expect_equal(out$gene_ratio, 31 / 72)
    expect_equal(out$size, 458)
    expect_equal(out$type, "ora")
    # GO descriptions get a capital first letter
    expect_equal(out$term, "Cell migration")

    desk <- data.frame(NAME = "HALLMARK_HYPOXIA", SIZE = 190, ES = 0.5, NES = 1.6,
                       `NOM p-val` = 0.001, `FDR q-val` = "0.02", check.names = FALSE)
    out <- standardize_enrichment(desk)
    expect_equal(out$padj, 0.02)
    expect_equal(out$pvalue, 0.001)

    two <- standardize_enrichment(list(KO = fgsea_fixture(1), Het = fgsea_fixture(2)))
    expect_equal(levels(two$comparison), c("KO", "Het"))
    expect_error(standardize_enrichment(data.frame(a = 1)), "Could not detect")
})

test_that("clean_term_names keeps non-MSigDB names and restores acronyms", {
    expect_equal(clean_term_names(c("HALLMARK_PI3K_AKT_MTOR_SIGNALING", "GOBP_DNA_REPAIR", "My own set")),
                 c("PI3K AKT MTOR signaling", "DNA repair", "My own set"))
    expect_equal(clean_term_names(c("extracellular matrix organization", "mRNA processing")),
                 c("Extracellular matrix organization", "mRNA processing"))
    expect_equal(replace_terms("Extracellular matrix organization", c("Extracellular matrix" = "ECM")),
                 "ECM organization")
    expect_error(replace_terms("a", "b"), "named character vector")
})

test_that("plot_enrichment_grid picks the top sets of each comparison", {
    res <- list(KO = fgsea_fixture(1), Het = fgsea_fixture(2))
    p <- plot_enrichment_grid(res, n = 1)
    src <- attr(p, "source_data")
    expect_setequal(unique(src$term), "TNFA signaling via NFKB")
    expect_equal(nrow(src), 2)
    p2 <- plot_enrichment_grid(res, terms = c("HALLMARK_HYPOXIA", "TNFA signaling via NFKB"))
    expect_setequal(unique(attr(p2, "source_data")$term), c("Hypoxia", "TNFA signaling via NFKB"))
    expect_error(plot_enrichment_grid(res, fill = "nes", size = "count"), "no `count`")
})

ora_fixture <- function(seed = 1) {
    set.seed(seed)
    k <- sample(5:40, 12)
    d <- data.frame(ID = paste0("GO", 1:12), Description = paste("gene set", 1:12), Count = k,
                    GeneRatio = paste0(k, "/200"), BgRatio = "100/20000",
                    pvalue = c(1e-8, 1e-6, stats::runif(10)^4))
    d$p.adjust <- stats::p.adjust(d$pvalue, "BH")
    d$geneID <- "A/B"
    d
}

test_that("plot_enrichment_dots draws one panel per result, each in its own palette", {
    # One result: a single ggplot, no strip
    p1 <- plot_enrichment_dots(ora_fixture(), n = 5)
    expect_s3_class(p1, "ggplot")
    src1 <- attr(p1, "source_data")
    expect_true(all(src1$significant))
    expect_lte(nrow(src1), 5)
    # Gene sets from the top: largest gene ratio first
    expect_equal(levels(p1$data$id), rev(src1$term_id))
    expect_true(all(diff(src1$x_gene_ratio) <= 0))

    two <- list(Up = ora_fixture(1), Down = ora_fixture(2))
    p2 <- plot_enrichment_dots(two, n = 4)
    expect_s3_class(p2, "patchwork")
    expect_equal(unique(attr(p2, "source_data")$panel), c("Up", "Down"))
    pals <- enrichment_panel_palettes(NULL, c("Up", "Down"))
    expect_equal(pals$Up, bt_palette("seq-purple"))
    expect_equal(pals$Down, bt_palette("seq-blue"))
    expect_equal(enrichment_panel_palettes(c(Down = "seq-red", Up = "seq-green"), c("Up", "Down"))$Down,
                 bt_palette("seq-red"))
    expect_equal(enrichment_panel_palettes(c("white", "red"), c("Up", "Down"))$Down, c("white", "red"))
    expect_error(enrichment_panel_palettes(list(Up = "seq-red"), c("Up", "Down")), "Down")

    # The panels, then the legends in two columns, in one flat layout (it
    # can be printed)
    parts <- patchwork_parts(p2)$plots
    p_up <- parts[[1]]
    p_down <- parts[[2]]
    grDevices::pdf(NULL)
    expect_no_error(print(p2))
    grDevices::dev.off()
    # Legend order is fixed: colour bar first, size second, in every panel
    for (pp in list(p_up, p_down)) {
        g <- pp$scales$get_scales("colour")$guide
        if (!identical(g, "none")) expect_equal(g$params$order, 1)
    }
    expect_equal(p_up$scales$get_scales("size")$guide$params$order, 2)
    # With shared limits, the size legend is drawn once
    expect_identical(p_down$scales$get_scales("size")$guide, "none")
    # Panels share colour limits
    expect_equal(p_up$scales$get_scales("colour")$limits, p_down$scales$get_scales("colour")$limits)
})

test_that("plot_enrichment_dots takes GSEA results, with counts from the leading edge", {
    p <- plot_enrichment_dots(fgsea_fixture(), n = 2)
    src <- attr(p, "source_data")
    expect_equal(src$count, c(2, 2))
    expect_equal(src$gene_ratio, c(2, 2) / 100)
    p2 <- plot_enrichment_dots(fgsea_fixture(), n = 2, x = "nes", order.by = "padj")
    expect_equal(attr(p2, "source_data")$x_nes, c(2.5, -2))
    expect_error(plot_enrichment_dots(ora_fixture(), x = "nes"), "no NES")
})

test_that("plot_enrichment_bars writes names opposite their bars", {
    p <- plot_enrichment_bars(fgsea_fixture(), n.up = 5, n.down = 5)
    d <- p$data
    expect_equal(d$term, c("Hypoxia", "TNFA signaling via NFKB"))
    expect_equal(sign(d$label_x), -sign(d$nes))
    expect_equal(d$hjust, c(0, 1))
    p2 <- plot_enrichment_bars(fgsea_fixture(), n.up = 3, n.down = 3, show.ns = TRUE)
    expect_lte(nrow(p2$data), 6)
    expect_true(any(p2$data$significant == "Not significant"))
    # Faded bars are explained in a legend
    expect_equal(p2$scales$get_scales("alpha")$name, "Significance")
})

test_that("label_row_height grows rows for wrapped names", {
    expect_equal(label_row_height(c("a", "b"), 0.2, 9), 0.2)
    expect_gt(label_row_height(c("a\nb\nc"), 0.2, 9), 0.2)
})

test_that("the running enrichment score follows GSEA's definition", {
    stats <- stats::setNames(c(3, 2, 1, -1, -2, -3), paste0("G", 1:6))
    run <- gsea_running_score(stats, c("G1", "G2"))
    # Hits weighted by |metric| (3 and 2 of 5), misses 1/4 each
    expect_equal(run$running_es, c(0.6, 1, 0.75, 0.5, 0.25, 0))
    expect_equal(run$in_set, c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE))
    expect_error(gsea_running_score(stats, "nope"), "None of")
})

test_that("plot_gsea_rank draws a gene set from an fgsea table with its ranked list", {
    set.seed(1)
    genes <- paste0("Gene", 1:500)
    stats <- sort(stats::setNames(stats::rnorm(500), genes), decreasing = TRUE)
    sets <- list(SET_UP = genes[c(1:20, 300)], SET_DOWN = genes[481:500])
    es_up <- max(gsea_running_score(stats, sets$SET_UP)$running_es)
    res <- data.frame(pathway = names(sets), pval = c(1e-4, 1e-3), padj = c(2e-4, 1e-3),
                      ES = c(es_up, -0.9), NES = c(2.1, -1.9), size = c(21, 20))
    p <- plot_gsea_rank(res, "SET_UP", stats = stats, gene.sets = sets, comparison = "KO vs WT")
    expect_s3_class(p, "patchwork")
    src <- attr(p, "source_data")
    expect_equal(sum(src$in_set), 21)
    # Title: the name, then the ID; subtitle: the comparison
    expect_equal(p$patches$annotation$title, "Set up\nSET_UP")
    expect_equal(p$patches$annotation$subtitle, "KO vs WT")
    # The first five genes of the leading edge (set genes up to the peak)
    # are named, by rank
    peak <- which.max(src$running_es)
    expect_equal(src$leading_edge, src$in_set & src$rank <= peak)
    expect_equal(src$gene[src$labelled], utils::head(src$gene[src$leading_edge], 5))
    # Named at the small text size (9 pt), in the metric panel (the last plot
    # of the patchwork, i.e. `p` itself)
    texts <- Filter(function(l) inherits(l$geom, "GeomText") && any(grepl("Gene", l$aes_params$label)),
                    p$layers)
    expect_equal(texts[[1]]$aes_params$size, pt_to_size(9))
    p_sel <- plot_gsea_rank(res, "SET_UP", stats = stats, gene.sets = sets, label = "selected",
                            label.features = "Gene300")
    expect_equal(attr(p_sel, "source_data")$gene[attr(p_sel, "source_data")$labelled], "Gene300")
    p_none <- plot_gsea_rank(res, "SET_UP", stats = stats, gene.sets = sets, label = "none", band = "none")
    expect_false(any(attr(p_none, "source_data")$labelled))
    expect_length(Filter(function(l) inherits(l$geom, "GeomRect"), p_none$patches$plots[[2]]$layers), 0)
    # Running score, gene positions over the metric band, ranked metric
    expect_length(p$patches$plots, 2)
    band <- Filter(function(l) inherits(l$geom, "GeomRect"), p$patches$plots[[2]]$layers)[[1]]$data
    expect_lt(nrow(band), 120)
    expect_equal(min(band$xmin), 0.5)
    expect_equal(max(band$xmax), 500.5)
    # Without the metric (the same for every gene set of a comparison)
    p_short <- plot_gsea_rank(res, "SET_UP", stats = stats, gene.sets = sets, show.metric = FALSE)
    expect_length(p_short$patches$plots, 1)
    # A reported ES that the ranked list does not give is flagged
    expect_warning(plot_gsea_rank(res, "SET_DOWN", stats = stats, gene.sets = sets), "differs")
    expect_error(plot_gsea_rank(res, "SET_UP"), "ranked list")
    expect_error(plot_gsea_rank(res, "nope", stats = stats, gene.sets = sets), "not found")
})

test_that("plot_enrichment_dots lines up panels side by side and collects legends", {
    ora <- function(n) {
        d <- data.frame(Description = paste("gene set", seq_len(n)), Count = 5 + seq_len(n),
                        GeneRatio = paste0(5 + seq_len(n), "/200"), BgRatio = "100/20000",
                        pvalue = 10^-(seq_len(n) + 2))
        d$p.adjust <- d$pvalue
        d
    }
    p <- plot_enrichment_dots(list(A = ora(6), B = ora(2)), n = 6, nrow = 1)
    # Both panels as tall as the longest: B gets 4 empty rows, without names
    parts <- patchwork_parts(p)$plots
    panels <- Filter(function(q) inherits(q$data, "data.frame") && "id" %in% names(q$data), parts)
    heights <- vapply(panels, function(q) as.numeric(q$theme$panel.heights), numeric(1))
    expect_equal(length(unique(round(heights, 6))), 1)
    b <- ggplot2::ggplot_build(panels[[2]])
    expect_equal(length(b$layout$panel_params[[1]]$y$get_breaks()), 2)
    # Legends in two columns beside the panels: colour bars, then sizes
    block <- parts[[length(parts)]]$layers[[2]]$geom_params$grob$children[[1]]
    expect_s3_class(block, "bt_legend_block")
    expect_length(block$columns, 2)
    expect_no_warning(save_plot(p, tempfile(fileext = ".pdf")))
    # Stacked: collected by patchwork, as before 1.0.0
    p_stacked <- plot_enrichment_dots(list(A = ora(6), B = ora(2)), n = 6, nrow = 1, legend.layout = "stacked")
    expect_identical(p_stacked$patches$layout$guides, "collect")
})
