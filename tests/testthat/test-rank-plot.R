# plot_rank(): ranking, ties, labels and source data

# Seven features with a tie of three (B, C, D at 2) and a signed metric
rank_fixture <- function() {
    data.frame(
        feature = c("A", "B", "C", "D", "E", "F", "G"),
        score = c(5, 2, 2, 2, 0.5, -1, -3),
        other = c(1, 3, 1, 2, 1, 1, 1),
        group = c("x", "x", "x", "y", "y", "y", "y")
    )
}

quiet_rank <- function(...) suppressMessages(plot_rank(...))

test_that("features are ranked from the largest value by default", {
    p <- quiet_rank(rank_fixture(), "score", feature.col = "feature", label = "none")
    src <- attr(p, "source_data")
    expect_equal(src$feature[order(src$rank)][1], "A")
    expect_equal(src$feature[order(src$rank)][7], "G")
    p2 <- quiet_rank(rank_fixture(), "score", feature.col = "feature", decreasing = FALSE, label = "none")
    src2 <- attr(p2, "source_data")
    expect_equal(src2$feature[src2$rank == 1], "G")
})

test_that("ties get distinct positions, a shared tied rank, and do not depend on row order", {
    d <- rank_fixture()
    p <- quiet_rank(d, "score", feature.col = "feature", label = "none")
    src <- attr(p, "source_data")
    tied <- src[src$feature %in% c("B", "C", "D"), ]
    # Positions 2, 3, 4; competition rank 2 for all three
    expect_setequal(tied$rank, 2:4)
    expect_true(all(tied$rank_tied == 2L))
    expect_equal(src$rank_tied[src$feature == "E"], 5L)
    # Ties broken by name when `value` is the ranked column
    expect_equal(tied$feature[order(tied$rank)], c("B", "C", "D"))
    # Ties broken by `value` first when it differs
    p_v <- quiet_rank(d, "score", value = "other", feature.col = "feature", label = "none")
    src_v <- attr(p_v, "source_data")
    tied_v <- src_v[src_v$feature %in% c("B", "C", "D"), ]
    expect_equal(tied_v$feature[order(tied_v$rank)], c("B", "D", "C"))
    # Same result whatever the row order
    p_rev <- quiet_rank(d[rev(seq_len(nrow(d))), ], "score", feature.col = "feature", label = "none")
    src_rev <- attr(p_rev, "source_data")
    expect_equal(src_rev[order(src_rev$rank), c("feature", "rank", "rank_tied")],
                 src[order(src$rank), c("feature", "rank", "rank_tied")],
                 ignore_attr = TRUE)
})

test_that("split.by ranks within each level", {
    p <- quiet_rank(rank_fixture(), "score", feature.col = "feature", split.by = "group", label = "none")
    src <- attr(p, "source_data")
    expect_equal(sort(src$rank[src$split == "x"]), 1:3)
    expect_equal(sort(src$rank[src$split == "y"]), 1:4)
    expect_equal(src$feature[src$split == "y" & src$rank == 1], "D")
})

test_that("rows missing rank.by or value are dropped with a message", {
    d <- rank_fixture()
    d$score[2] <- NA
    expect_message(p <- plot_rank(d, "score", feature.col = "feature", label = "none"), "1 row")
    expect_false("B" %in% attr(p, "source_data")$feature)
})

test_that("label.side auto labels both ends of a signed metric, the first ranks otherwise", {
    d <- rank_fixture()
    p <- quiet_rank(d, "score", feature.col = "feature", label.n = 1)
    src <- attr(p, "source_data")
    expect_setequal(src$feature[src$labelled], c("A", "G"))
    d$score <- abs(d$score) + 1
    p2 <- quiet_rank(d, "score", feature.col = "feature", label.n = 1)
    src2 <- attr(p2, "source_data")
    expect_equal(src2$feature[src2$labelled], "A")
    p3 <- quiet_rank(d, "score", feature.col = "feature", label.n = 1, label.side = "last")
    src3 <- attr(p3, "source_data")
    expect_equal(src3$feature[src3$labelled], src3$feature[src3$rank == 7])
})

test_that("ties at the label.n cutoff are reported", {
    expect_message(
        plot_rank(rank_fixture(), "score", feature.col = "feature", label.n = 2, label.side = "first"),
        "tied with a labelled one"
    )
})

test_that("selected labels: matched exactly, missing ones warned", {
    expect_warning(
        p <- quiet_rank(rank_fixture(), "score", feature.col = "feature", label = "selected",
                        label.features = c("C", "Z")),
        "not found"
    )
    src <- attr(p, "source_data")
    expect_equal(src$feature[src$labelled], "C")
    expect_error(plot_rank(rank_fixture(), "score", feature.col = "feature", label = "selected"),
                 "label.features")
})

test_that("inputs are checked", {
    d <- rank_fixture()
    expect_error(plot_rank(d, "group", feature.col = "feature"), "numeric")
    expect_error(plot_rank(d, "nope", feature.col = "feature"), "column")
    expect_error(plot_rank(d, "score"), "row names")
    expect_error(plot_rank(d, "score", feature.col = "feature", label.n = -1), "label.n")
})

test_that("colours: numeric, discrete and the labelled highlight", {
    d <- rank_fixture()
    p <- quiet_rank(d, "score", feature.col = "feature", color.by = "other", label = "none")
    expect_true("color" %in% names(attr(p, "source_data")))
    expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
    p2 <- quiet_rank(d, "score", feature.col = "feature", color.by = "group", label = "none")
    expect_s3_class(ggplot2::ggplot_build(p2), "ggplot_built")
    expect_error(plot_rank(d, "score", feature.col = "feature", cols = c("red", "blue", "green")), "two colours")
})

test_that("labels are placed inside the panel, clear of each other", {
    set.seed(1)
    d <- data.frame(feature = paste0("Feature_", 1:300), score = rnorm(300))
    p <- quiet_rank(d, "score", feature.col = "feature", label.n = 5, bg = "white")
    built <- ggplot2::ggplot_build(p)
    # The text layer is the last; its labels are the 10 features at the ends
    txt <- built$data[[length(built$data)]]
    expect_equal(nrow(txt), 10)
    lims <- built$layout$panel_params[[1]]
    expect_true(all(txt$x > lims$x.range[1] & txt$x < lims$x.range[2]))
    expect_true(all(txt$y > lims$y.range[1] & txt$y < lims$y.range[2]))
})

test_that("rank_place_labels keeps boxes apart and off the points", {
    # Ten labels around one corner of a falling curve
    x <- seq(0.1, 0.3, length.out = 10)
    y <- seq(1.6, 1.2, length.out = 10)
    obs_x <- c(x, seq(0.3, 2.4, length.out = 200))
    obs_y <- c(y, seq(1.2, 0.2, length.out = 200))
    res <- rank_place_labels(
        x = x, y = y, w = rep(0.4, 10), h = rep(0.1, 10),
        side = rep(1, 10), vside = rep(1, 10), order = 1:10,
        obs_x = obs_x, obs_y = obs_y, obs_r = rep(0.02, 210), self = 1:10,
        panel.width = 2.5, panel.height = 1.8
    )
    expect_true(all(res$level <= 3))
    # No two boxes overlap
    for (a in 1:9) {
        for (b in (a + 1):10) {
            sep <- abs(res$cx[a] - res$cx[b]) >= 0.4 | abs(res$cy[a] - res$cy[b]) >= 0.1
            expect_true(sep)
        }
    }
})

test_that("wrap_name breaks long names after separators", {
    expect_equal(wrap_name("EPITHELIAL_MESENCHYMAL_TRANSITION", 25), "EPITHELIAL_MESENCHYMAL_\nTRANSITION")
    expect_equal(wrap_name("short", 25), "short")
    expect_equal(wrap_name("3-hydroxy-3-methylglutaryl-CoA", 15), "3-hydroxy-3-\nmethylglutaryl-\nCoA")
    expect_equal(wrap_name("anything", NULL), "anything")
})

test_that("segments_cross ignores touching ends and crossings near the start", {
    expect_true(segments_cross(0, 0, 1, 1, 0, 1, 1, 0))
    expect_false(segments_cross(0, 0, 1, 1, 0, 0, 1, 0))
    expect_false(segments_cross(0, 0, 1, 1, 0, 0.1, 0.1, 0, min.dist = 0.2))
})
