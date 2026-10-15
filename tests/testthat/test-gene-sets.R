# format_gene_sets(): list <-> GMT <-> GeneSetCollection

test_that("format_gene_sets converts a list to GMT lines and back", {
    skip_if_not_installed("GSEABase")
    sets <- list(SET_A = c("CD3E", "CD3D"), SET_B = c("MS4A1", "CD79A", "CD19"))
    gmt <- format_gene_sets(sets, to = "gmt")
    expect_equal(gmt, c("SET_A\t\tCD3E\tCD3D", "SET_B\t\tMS4A1\tCD79A\tCD19"))
    f <- tempfile(fileext = ".gmt")
    format_gene_sets(sets, to = "gmt", file = f)
    expect_equal(format_gene_sets(f, to = "list"), sets)
    expect_s4_class(format_gene_sets(sets), "GeneSetCollection")
})

test_that("format_gene_sets keeps a gene listed twice in a set once", {
    skip_if_not_installed("GSEABase")
    expect_equal(format_gene_sets(list(S = c("A", "B", "A")), to = "gmt"), "S\t\tA\tB")
})
