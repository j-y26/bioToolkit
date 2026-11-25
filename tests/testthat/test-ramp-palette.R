# Single-hue ramps from given colours

test_that("a single colour gives a light-to-dark vector in its hue", {
    skip_if_not_installed("colorspace")
    ramp <- ramp_palette("#1F77B4")
    expect_type(ramp, "character")
    expect_length(ramp, 9)
    hcl <- colorspace::coords(methods::as(colorspace::hex2RGB(ramp), "polarLUV"))
    # Lightness falls along the ramp, from about 96 to 30
    expect_true(all(diff(hcl[, "L"]) < 0))
    expect_equal(unname(hcl[c(1, 9), "L"]), c(96, 30), tolerance = 0.5)
    # Hue stays that of the input colour (blue, H about 245)
    input_h <- colorspace::coords(methods::as(colorspace::hex2RGB("#1F77B4"), "polarLUV"))[, "H"]
    expect_true(all(abs(hcl[-1, "H"] - input_h) < 2))
})

test_that("vectors and lists give named lists of ramps", {
    skip_if_not_installed("colorspace")
    named <- ramp_palette(c(Fibroblast = "#1F77B4", Endothelial = "firebrick"), n = 5)
    expect_named(named, c("Fibroblast", "Endothelial"))
    expect_true(all(lengths(named) == 5))
    unnamed <- ramp_palette(c("#1F77B4", "red"))
    expect_named(unnamed, c("#1F77B4", "red"))
    from_list <- ramp_palette(list(a = "#1F77B4", b = "red"))
    expect_named(from_list, c("a", "b"))
    expect_equal(from_list$a, ramp_palette("#1F77B4"))
})

test_that("grey stays grey, reverse flips, bad input stops", {
    skip_if_not_installed("colorspace")
    grey <- ramp_palette("grey50", n = 3)
    rgb <- grDevices::col2rgb(grey)
    expect_true(all(apply(rgb, 2, function(x) diff(range(x)) <= 1)))
    expect_equal(ramp_palette("red", reverse = TRUE), rev(ramp_palette("red")))
    expect_error(ramp_palette("notacolour"), "Not a colour")
    expect_error(ramp_palette("red", n = 1), "`n`")
    expect_error(ramp_palette("red", peak = 1), "`peak`")
})
