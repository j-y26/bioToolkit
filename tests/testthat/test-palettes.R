# Palettes (roadmap item F4)

test_that("every palette is well formed", {
    pals <- list_palettes()
    expect_equal(anyDuplicated(tolower(pals$name)), 0L)
    expect_true(all(pals$type %in% c("qualitative", "sequential", "diverging")))
    expect_true(all(nzchar(pals$source)))
    for (nm in pals$name) {
        cols <- bt_palette(nm)
        expect_true(all(grepl("^#[0-9A-F]{6}$", cols)), info = nm)
        # Group colours must all be distinct
        if (pals$type[pals$name == nm] == "qualitative") expect_equal(anyDuplicated(cols), 0L, info = nm)
    }
})

test_that("scanpy-20 is the user's reference list", {
    expect_equal(bt_palette("scanpy-20")[c(1, 8, 11, 20)], c("#1F77B4", "#7F7F7F", "#AEC7E8", "#9EDAE5"))
    expect_length(bt_palette("scanpy-20"), 20)
})

test_that("colours can be picked, reordered and named by index", {
    cols <- bt_palette("scanpy-20")[c(3, 1)]
    expect_equal(cols, c("#2CA02C", "#1F77B4"))
    expect_equal(bt_palette("SCANPY-20"), bt_palette("scanpy-20"))
    expect_equal(bt_palette("npg", 3, reverse = TRUE), rev(bt_palette("npg")[1:3]))
})

test_that("qualitative palettes are never interpolated; continuous ones are", {
    expect_equal(bt_palette("npg", 4), bt_palette("npg")[1:4])
    expect_error(bt_palette("npg", 11), "has 10 colours")
    expect_length(bt_palette("tol-burd", 101), 101)
    # An odd number of colours keeps the centre of a diverging palette
    expect_equal(bt_palette("tol-burd", 9)[5], "#F7F7F7")
    expect_error(bt_palette("npgg"), "Unknown palette")
})

test_that("default colours: scanpy-20 up to 20 groups, scanpy-102 beyond", {
    expect_equal(default_colours(5), bt_palette("scanpy-20", 5))
    expect_equal(default_colours(20), bt_palette("scanpy-20"))
    expect_equal(default_colours(30), bt_palette("scanpy-102", 30))
    expect_message(default_colours(110), "More than 102")
})

test_that("show_palettes numbers each colour, wraps long palettes and ramps continuous ones", {
    p <- show_palettes(c("scanpy-20", "scanpy-102", "viridis"), ncol = 20)
    expect_s3_class(p, "ggplot")
    d <- p$data
    s20 <- d[d$palette == "scanpy-20", ]
    expect_equal(s20$label, as.character(1:20))
    expect_equal(s20$fill, bt_palette("scanpy-20"))
    expect_equal(max(d$y[d$palette == "scanpy-102"]), 6)
    expect_true(all(d$label[d$palette == "viridis"] == ""))
    # Not a fixed-size plot: it fills the plotting window
    expect_null(p$theme$panel.widths)
    expect_error(show_palettes("nope"), "Unknown palette")
})

test_that("show_palettes draws colours and lists of palettes and colours", {
    cell_cols <- c(T = "#1F77B4", B = "steelblue", NK = "red")
    d <- show_palettes(cell_cols)$data
    expect_equal(levels(d$palette), "cell_cols")
    expect_equal(d$fill, c("#1F77B4", "#4682B4", "#FF0000"))
    expect_equal(d$label, c("1", "2", "3"))
    expect_equal(levels(show_palettes(c("red", "blue"))$data$palette), "colours")
    d <- show_palettes(list(mine = cell_cols, "npg", c("gold", "navy"), ramp = viridisLite::mako(200)))$data
    expect_equal(levels(d$palette), c("mine", "npg", "3", "ramp"))
    expect_equal(d$fill[d$palette == "npg"], bt_palette("npg"))
    expect_true(all(d$label[d$palette == "ramp"] == ""))
    expect_error(show_palettes(c("npg", "red")), "mixes")
    expect_error(show_palettes(c("npg", "nope")), "nope")
    expect_error(show_palettes(list(c("red", "zz"))), "Invalid colour")
    expect_error(show_palettes(list(1:3)), "not a palette name")
    expect_warning(p <- show_palettes(names = "npg"), class = "deprecatedWarning")
    expect_equal(levels(p$data$palette), "npg")
})

test_that("check_palette flags colours that are hard to tell apart", {
    skip_if_not_installed("colorspace")
    expect_silent(res <- check_palette("okabe-ito"))
    expect_equal(nrow(res), 0)
    # The red and green of D3 look alike with deuteranopia (Delta E about 7)
    expect_warning(res <- check_palette(c("#D62728", "#2CA02C")), "hard to tell apart")
    expect_true(any(res$vision %in% c("deutan", "protan")))
})

test_that("single-hue palettes span the same lightness, light to dark", {
    seq_pals <- grep("^seq-", list_palettes("sequential")$name, value = TRUE)
    expect_setequal(seq_pals, paste0("seq-", c("purple", "blue", "teal", "green", "orange", "red", "magenta")))
    lum <- function(x) grDevices::convertColor(t(grDevices::col2rgb(x)) / 255, from = "sRGB", to = "Luv")[, 1]
    for (nm in seq_pals) {
        l <- lum(bt_palette(nm))
        expect_true(all(diff(l) < 0), info = nm)
        # Tinted, never white, at the light end; coloured, never black, at the dark end
        expect_true(l[1] <= 80.5 && l[1] > 75, info = nm)
        expect_true(l[length(l)] >= 29.5 && l[length(l)] < 36, info = nm)
    }
})

test_that("any palette name with \"_r\" is that palette reversed, without being stored", {
    expect_equal(bt_palette("viridis_r"), rev(bt_palette("viridis")))
    expect_equal(bt_palette("plasma_r"), rev(bt_palette("plasma")))
    expect_equal(bt_palette("NPG_R"), rev(bt_palette("npg")))
    expect_equal(bt_palette("seq-blue_r", 3), rev(bt_palette("seq-blue", 3)))
    expect_false(any(grepl("_r$", list_palettes()$name)))
    expect_null(find_palette("bioToolkit"))
    expect_null(find_palette("nope_r"))
    expect_error(bt_palette("nope_r"), "Unknown palette")
    skip_if_not_installed("RColorBrewer")
    expect_equal(resolve_continuous_palette("BuPu_r", FALSE), rev(RColorBrewer::brewer.pal(9, "BuPu")))
})
