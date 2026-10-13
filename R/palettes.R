# Colour palettes (roadmap item F4)
#
# Design
# ------
# Kept deliberately simple:
# * The palettes are a plain named list, `palette_list()`. Each entry has its
#   colours, its type and its source.
# * Users get colours with `bt_palette("scanpy-20")`, a plain vector of hex
#   codes, and pick, reorder and name them with base R:
#       cols <- bt_palette("scanpy-20")[c(1, 3, 2, 10)]
#       names(cols) <- c("T cells", "B cells", "Monocytes", "NK cells")
#   `show_palettes()` shows every colour with its index, for exactly this.
# * Plot functions have good defaults and one argument, `cols`, to change
#   them:
#     groups      `cols = NULL` gives the default palette; a palette name or
#                 a vector of colours (named by group, or in level order)
#                 replaces it. Default: "scanpy-20" up to 20 groups,
#                 "scanpy-102" beyond.
#     continuous  `cols = NULL` gives the plot's default: "plasma_r" for
#                 expression on embeddings, "viridis_r" for unscaled dot
#                 plots, "seq-purple" (etc.) for significance, "blue-magenta"
#                 for diverging values; a palette name, or colours from low
#                 to high, e.g. c(low, high) or c(low, mid, high).
# * Any palette name with "_r" appended is that palette reversed
#   (find_palette()); reversed palettes are not stored separately.
#
# Adding a palette
# ----------------
# Add an entry to `palette_list()` with its source. Copy hex values from the
# primary source, or from a package that reproduces it, never retype them from
# memory. tests/testthat/test-palettes.R checks every entry.
#
# Sources
# -------
# Colour values are taken from published palettes and credited in each
# entry's `source`. Hex codes are data, not code; no code is copied from the
# packages they were checked against (ggsci, khroma, scanpy, grDevices,
# viridisLite).


#' The palettes of the package
#'
#' @return A named list; each entry has `colours` (hex codes), `type`
#'         ("qualitative", "sequential" or "diverging") and `source`
#'
#' @keywords internal
palette_list <- function() {
    # Colours are stored as "#RRGGBB": viridisLite adds an alpha channel
    pal <- function(colours, type, source) {
        hex <- grDevices::rgb(t(grDevices::col2rgb(colours)), maxColorValue = 255)
        list(colours = hex, type = type, source = source)
    }
    tol <- "Paul Tol (2021), Colour schemes, SRON/EPS/TN/09-002, https://personal.sron.nl/~pault/"
    ggsci <- function(journal) {
        paste0("ggsci (N. Xiao), after the colours of ", journal, " figures; not an official journal palette")
    }
    viridis <- "viridisLite (Garnier et al.); matplotlib colour maps (van der Walt & Smith)"
    list(
        # Qualitative: for groups ------------------------------------------------
        `scanpy-20` = pal(
            c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
              "#7f7f7f", "#bcbd22", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896",
              "#c5b0d5", "#c49c94", "#f7b6d2", "#c7c7c7", "#dbdb8d", "#9edae5"),
            "qualitative",
            paste("matplotlib 'tab10' (the colours scanpy uses for up to 10 groups) followed by",
                  "the 10 lighter 'tab20' tones; D3 category10 / category20 (M. Bostock).",
                  "Default for up to 20 groups")
        ),
        `scanpy-20-cb` = pal(
            c("#1F77B4", "#FF7F0E", "#279E68", "#D62728", "#AA40FC", "#8C564B", "#E377C2",
              "#B5BD61", "#17BECF", "#AEC7E8", "#FFBB78", "#98DF8A", "#FF9896", "#C5B0D5",
              "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5", "#AD494A", "#8C6D31"),
            "qualitative",
            paste("scanpy default_20 for 11-20 groups (Wolf et al. 2018, Genome Biol 19:15):",
                  "tab20 without greys, three colours adjusted for colour-blindness")
        ),
        `scanpy-102` = pal(
            c("#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
              "#FFDBE5", "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF",
              "#997D87", "#5A0007", "#809693", "#6A3A4C", "#1B4400", "#4FC601", "#3B5DFF",
              "#4A3B53", "#FF2F80", "#61615A", "#BA0900", "#6B7900", "#00C2A0", "#FFAA92",
              "#FF90C9", "#B903AA", "#D16100", "#DDEFFF", "#000035", "#7B4F4B", "#A1C299",
              "#300018", "#0AA6D8", "#013349", "#00846F", "#372101", "#FFB500", "#C2FFED",
              "#A079BF", "#CC0744", "#C0B9B2", "#C2FF99", "#001E09", "#00489C", "#6F0062",
              "#0CBD66", "#EEC3FF", "#456D75", "#B77B68", "#7A87A1", "#788D66", "#885578",
              "#FAD09F", "#FF8A9A", "#D157A0", "#BEC459", "#456648", "#0086ED", "#886F4C",
              "#34362D", "#B4A8BD", "#00A6AA", "#452C2C", "#636375", "#A3C8C9", "#FF913F",
              "#938A81", "#575329", "#00FECF", "#B05B6F", "#8CD0FF", "#3B9700", "#04F757",
              "#C8A1A1", "#1E6E00", "#7900D7", "#A77500", "#6367A9", "#A05837", "#6B002C",
              "#772600", "#D790FF", "#9B9700", "#549E79", "#FFF69F", "#201625", "#72418F",
              "#BC23FF", "#99ADC0", "#3A2465", "#922329", "#5B4534", "#FDE8DC", "#404E55",
              "#0089A3", "#CB7E98", "#A4E804", "#324E72"),
            "qualitative",
            "scanpy default_102 ('godsnot_102', maximally distinct colours, black removed). Default above 20 groups"
        ),
        `okabe-ito` = pal(
            # Black, first in the original, is moved last so that it is used only
            # for an eighth group: black points hide outlines and labels
            c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000"),
            "qualitative",
            "Okabe & Ito (2008), Color Universal Design; Wong (2011) Nat Methods 8:441. Colour-blind safe"
        ),
        # Tol's grey for missing data is left out of his palettes: grey is kept
        # for cells that are not coloured (`na.value`)
        `tol-bright` = pal(c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377"),
                           "qualitative", paste(tol, "Colour-blind safe", sep = ". ")),
        `tol-vibrant` = pal(c("#EE7733", "#0077BB", "#33BBEE", "#EE3377", "#CC3311", "#009988"),
                            "qualitative", paste(tol, "Colour-blind safe", sep = ". ")),
        `tol-muted` = pal(c("#CC6677", "#332288", "#DDCC77", "#117733", "#88CCEE", "#882255",
                            "#44AA99", "#999933", "#AA4499"),
                          "qualitative", paste(tol, "Colour-blind safe", sep = ". ")),
        `tol-light` = pal(c("#77AADD", "#EE8866", "#EEDD88", "#FFAABB", "#99DDFF", "#44BB99",
                            "#BBCC33", "#AAAA00"),
                          "qualitative", paste(tol, "Colour-blind safe; light, for annotation strips", sep = ". ")),
        `tol-high-contrast` = pal(c("#004488", "#DDAA33", "#BB5566"),
                                  "qualitative", paste(tol, "Colour-blind safe; distinct in greyscale", sep = ". ")),
        npg = pal(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
                    "#91D1C2", "#DC0000", "#7E6148", "#B09C85"),
                  "qualitative", ggsci("Nature Publishing Group")),
        aaas = pal(c("#3B4992", "#EE0000", "#008B45", "#631879", "#008280", "#BB0021",
                     "#5F559B", "#A20056", "#808180", "#1B1919"),
                   "qualitative", ggsci("Science (AAAS)")),
        nejm = pal(c("#BC3C29", "#0072B5", "#E18727", "#20854E", "#7876B1", "#6F99AD",
                     "#FFDC91", "#EE4C97"),
                   "qualitative", ggsci("New England Journal of Medicine")),
        lancet = pal(c("#00468B", "#ED0000", "#42B540", "#0099B4", "#925E9F", "#FDAF91",
                       "#AD002A", "#ADB6B6", "#1B1919"),
                     "qualitative", ggsci("The Lancet")),
        jama = pal(c("#374E55", "#DF8F44", "#00A1D5", "#B24745", "#79AF97", "#6A6599", "#80796B"),
                   "qualitative", ggsci("JAMA")),
        jco = pal(c("#0073C2", "#EFC000", "#868686", "#CD534C", "#7AA6DC", "#003C67",
                    "#8F7700", "#3B3B3B", "#A73030", "#4A6990"),
                  "qualitative", ggsci("Journal of Clinical Oncology")),
        `tableau-10` = pal(grDevices::palette.colors(palette = "Tableau 10"), "qualitative",
                           "Tableau 10 (Tableau Software), as in grDevices::palette.colors()"),
        `polychrome-36` = pal(grDevices::palette.colors(palette = "Polychrome 36"), "qualitative",
                              "Coombes et al. (2019) J Stat Softw 90(1), as in grDevices::palette.colors(); Seurat's 'polychrome'"),
        alphabet = pal(grDevices::palette.colors(palette = "Alphabet"), "qualitative",
                       "Green-Armytage (2010), a colour alphabet, as in grDevices::palette.colors(); Seurat's 'alphabet'"),
        paired = pal(grDevices::palette.colors(palette = "Paired"), "qualitative",
                     "ColorBrewer 'Paired' (C. Brewer), as in grDevices::palette.colors(); light / dark pairs"),
        # Sequential: low to high values -------------------------------------------
        # Single-hue ramps (pale tint to deep colour) for any continuous,
        # one-directional value: significance, expression, scores. All are
        # trimmed to the same lightness range (see hcl_single_hue()), so that
        # the same lightness means the same value whatever the hue: panels of
        # one figure drawn in different hues stay comparable
        `seq-purple` = pal(hcl_single_hue("Purples 3"), "sequential", hcl_source("Purples 3")),
        `seq-blue` = pal(hcl_single_hue("Blues 3"), "sequential", hcl_source("Blues 3")),
        `seq-teal` = pal(hcl_single_hue("Mint"), "sequential", hcl_source("Mint")),
        `seq-green` = pal(hcl_single_hue("Greens 3"), "sequential", hcl_source("Greens 3")),
        `seq-orange` = pal(hcl_single_hue("Oranges"), "sequential", hcl_source("Oranges")),
        `seq-red` = pal(hcl_single_hue("Reds 3"), "sequential", hcl_source("Reds 3")),
        `seq-magenta` = pal(hcl_single_hue("Magenta"), "sequential", hcl_source("Magenta")),
        viridis = pal(viridisLite::viridis(256), "sequential", viridis),
        magma = pal(viridisLite::magma(256), "sequential", viridis),
        inferno = pal(viridisLite::inferno(256), "sequential", viridis),
        plasma = pal(viridisLite::plasma(256), "sequential", viridis),
        cividis = pal(viridisLite::cividis(256), "sequential", "Nu\u00f1ez et al. (2018) PLoS ONE 13:e0199239; via viridisLite"),
        mako = pal(viridisLite::mako(256), "sequential", "seaborn, via viridisLite"),
        rocket = pal(viridisLite::rocket(256), "sequential", "seaborn, via viridisLite"),
        `tol-iridescent` = pal(c("#FEFBE9", "#FCF7D5", "#F5F3C1", "#EAF0B5", "#DDECBF", "#D0E7CA",
                                 "#C2E3D2", "#B5DDD8", "#A8D8DC", "#9BD2E1", "#8DCBE4", "#81C4E7",
                                 "#7BBCE7", "#7EB2E4", "#88A5DD", "#9398D2", "#9B8AC4", "#9D7DB2",
                                 "#9A709E", "#906388", "#805770", "#684957", "#46353A"),
                               "sequential", tol),
        `tol-ylorbr` = pal(c("#FFFFE5", "#FFF7BC", "#FEE391", "#FEC44F", "#FB9A29", "#EC7014",
                             "#CC4C02", "#993404", "#662506"),
                           "sequential", tol),
        # Diverging: values on both sides of a centre ------------------------------
        `blue-magenta` = pal(c("#005C9B", "#0075A6", "#3D92B7", "#88B1CA", "#BECFDB", "#EBEBEB",
                               "#E0C6D1", "#D19CB4", "#C16F98", "#AE3C7B", "#9B0062"),
                             "diverging",
                             "bioToolkit; the sides stay distinct with red-green colour-blindness. Default diverging"),
        `blue-purple` = pal(c("#0B3C7D", "#175196", "#326DB0", "#6B97CA", "#ACC2DF", "#F0EEF2",
                              "#D7BCD8", "#BA84BB", "#964097", "#752379", "#561A5E"),
                            "diverging",
                            "bioToolkit; truer purple, but both sides look blue with red-green colour-blindness"),
        `tol-sunset` = pal(c("#364B9A", "#4A7BB7", "#6EA6CD", "#98CAE1", "#C2E4EF", "#EAECCC",
                             "#FEDA8B", "#FDB366", "#F67E4B", "#DD3D2D", "#A50026"),
                           "diverging", tol),
        `tol-burd` = pal(c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7", "#FDDBC7",
                           "#F4A582", "#D6604D", "#B2182B"),
                         "diverging", paste(tol, "(same colours as ColorBrewer RdBu, reversed)")),
        `tol-prgn` = pal(c("#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7", "#D9F0D3",
                           "#ACD39E", "#5AAE61", "#1B7837"),
                         "diverging", paste(tol, "(same colours as ColorBrewer PRGn)"))
    )
}


#' A single-hue sequential palette, trimmed to a lightness range
#'
#' Takes one of base R's HCL sequential palettes (\code{grDevices::hcl.colors()},
#' Zeileis et al. 2020) and keeps the part between lightness `light` and
#' `dark` (CIE L*, 0 = black, 100 = white). The untrimmed palettes run from
#' near white, which vanishes on a white or transparent background, to near
#' black, which hides the hue; L* 80 is a clear tint (dots of the lowest
#' values stay visible on white) and L* 30 still clearly
#' coloured. Every palette trimmed this way spans the same lightness, so that
#' different hues encode the same value with the same lightness.
#'
#' @param name Name of an \code{hcl.colors()} sequential palette
#'
#' @param n Number of colours, default is 9
#'
#' @param light,dark Lightness (L*) of the light and dark ends, default 80
#'                   and 30
#'
#' @return `n` colours from light to dark, taken from the palette (not
#'         interpolated)
#'
#' @keywords internal
hcl_single_hue <- function(name,
                           n = 9,
                           light = 80,
                           dark = 30) {
    # Light to dark, finely sampled; lightness L* of each colour (CIELUV L*,
    # the L of HCL)
    full <- grDevices::hcl.colors(501, name, rev = TRUE)
    lum <- grDevices::convertColor(t(grDevices::col2rgb(full)) / 255, from = "sRGB", to = "Luv")[, 1]
    keep <- full[lum <= light & lum >= dark]
    keep[round(seq(1, length(keep), length.out = n))]
}


#' Source text of a palette trimmed by hcl_single_hue()
#'
#' @param name Name of the \code{hcl.colors()} palette
#'
#' @keywords internal
hcl_source <- function(name) {
    paste0("grDevices::hcl.colors(\"", name, "\") (Zeileis et al. 2020, J Stat Softw 96(1)), ",
           "trimmed to lightness L* 80 to 30")
}


#' Look up a palette by name, ignoring case
#'
#' A name ending in "_r" (matplotlib's convention) is the palette of the
#' name without it, reversed: "viridis_r", "plasma_r", "npg_r"... Reversed
#' palettes are not stored, so every palette has one.
#'
#' @param name Palette name
#'
#' @return The palette entry (colours reversed for an "_r" name), or NULL
#'         when there is none
#'
#' @keywords internal
find_palette <- function(name) {
    if (!(is.character(name) && length(name) == 1 && !is.na(name))) return(NULL)
    pals <- palette_list()
    reverse <- grepl("_r$", name, ignore.case = TRUE)
    base <- if (reverse) sub("_r$", "", name, ignore.case = TRUE) else name
    i <- match(tolower(base), tolower(names(pals)))
    if (is.na(i)) return(NULL)
    entry <- pals[[i]]
    if (reverse) entry$colours <- rev(entry$colours)
    entry
}


#' Whether strings are valid R colours
#'
#' @param x A character vector
#'
#' @return A logical vector
#'
#' @keywords internal
is_colour <- function(x) {
    vapply(x, function(col) {
        !is.na(col) && !inherits(try(grDevices::col2rgb(col), silent = TRUE), "try-error")
    }, logical(1), USE.NAMES = FALSE)
}


#' @title Colours of a palette
#'
#' @description Returns the colours of a palette of the package as a plain
#'              vector of hex codes, to pick from, reorder, name and pass to
#'              the `cols` argument of any plot function (or to
#'              \code{scale_fill_manual()}, ComplexHeatmap, base graphics...).
#'              Use \code{\link{show_palettes}} to see the colours with their
#'              index.
#'
#' @param name Palette name (case is ignored), see \code{\link{list_palettes}}.
#'             Add "_r" to any name for the palette reversed, e.g.
#'             "viridis_r" (yellow to dark purple) or "npg_r"
#'
#' @param n Number of colours. If NULL (default), all colours of the palette.
#'          Qualitative palettes return their first `n` colours and are never
#'          interpolated, since a colour in between two group colours would
#'          read as an intermediate group; sequential and diverging palettes
#'          are interpolated to `n` colours
#'
#' @param reverse Whether to reverse the order of the colours, default is FALSE
#'
#' @return A character vector of hex colours
#'
#' @examples
#' bt_palette("scanpy-20")
#'
#' # Pick and reorder colours by index, then name them by cell type
#' cols <- bt_palette("scanpy-20")[c(1, 3, 2, 10)]
#' names(cols) <- c("T cells", "B cells", "Monocytes", "NK cells")
#' cols
#'
#' # Low, mid and high colours of a diverging scale
#' bt_palette("tol-burd")[c(1, 5, 9)]
#'
#' @seealso \code{\link{show_palettes}}, \code{\link{list_palettes}}
#'
#' @export
#'
bt_palette <- function(name,
                       n = NULL,
                       reverse = FALSE) {
    entry <- find_palette(name)
    if (is.null(entry)) {
        stop("Unknown palette '", paste(name, collapse = " "), "'. Available: ",
             paste(names(palette_list()), collapse = ", "), call. = FALSE)
    }
    cols <- entry$colours
    if (!is.null(n)) {
        if (!(is.numeric(n) && length(n) == 1 && !is.na(n) && n >= 0)) {
            stop("`n` must be a single non-negative number", call. = FALSE)
        }
        if (entry$type == "qualitative") {
            if (n > length(cols)) {
                stop("Palette '", name, "' has ", length(cols), " colours but ", n,
                     " were asked for; use a longer palette such as 'scanpy-102'", call. = FALSE)
            }
            cols <- cols[seq_len(n)]
        } else {
            # Interpolated in CIELAB, as by ggplot2's gradient scales
            cols <- toupper(scales::gradient_n_pal(cols)(seq(0, 1, length.out = n)))
        }
    }
    if (reverse) rev(cols) else cols
}


#' @title List the palettes of the package
#'
#' @param type Types to list: any of "qualitative" (for groups), "sequential"
#'             and "diverging". If NULL (default), all
#'
#' @return A data frame with the palette `name`, its `type`, its number of
#'         colours `n` and its `source`, for citing
#'
#' @examples
#' list_palettes("qualitative")
#'
#' @export
#'
list_palettes <- function(type = NULL) {
    pals <- palette_list()
    if (!is.null(type)) {
        type <- match.arg(type, c("qualitative", "sequential", "diverging"), several.ok = TRUE)
        pals <- Filter(function(p) p$type %in% type, pals)
    }
    data.frame(
        name = names(pals),
        type = vapply(pals, `[[`, "", "type"),
        n = vapply(pals, function(p) length(p$colours), integer(1)),
        source = vapply(pals, `[[`, "", "source"),
        row.names = NULL
    )
}


#' @title Show palettes with the index of each colour
#'
#' @description Draws each palette as a row of swatches numbered by their
#'              index, so that colours can be picked and reordered with
#'              `bt_palette(name)[indices]`. Long palettes wrap onto several
#'              rows. Palettes defined by many colours (the viridis family)
#'              are drawn as a smooth ramp instead. Names ending in "_r" show
#'              a palette reversed. The plot
#'              fills the plotting window; it is meant for browsing, not for
#'              publication.
#'
#' @param names Palette names. If NULL (default), all palettes of `type`
#'
#' @param type Types to show when `names` is NULL: any of "qualitative",
#'             "sequential" and "diverging". If NULL (default), all
#'
#' @param ncol Number of swatches per row, default is 20
#'
#' @return A ggplot object
#'
#' @examples
#' show_palettes("scanpy-20")
#' show_palettes(type = "qualitative")
#'
#' @export
#'
show_palettes <- function(names = NULL,
                          type = NULL,
                          ncol = 20) {
    if (is.null(names)) names <- list_palettes(type)$name
    rows <- lapply(names, function(nm) {
        cols <- bt_palette(nm)
        if (length(cols) > 102) {
            # A ramp of 100 tiles spanning the width of a row, without indices
            cols <- bt_palette(nm, 100)
            return(data.frame(palette = nm, x = 0.5 + (seq_along(cols) - 0.5) * ncol / 100,
                              y = 1, width = ncol / 100, fill = cols, label = ""))
        }
        i <- seq_along(cols)
        data.frame(palette = nm, x = (i - 1) %% ncol + 1, y = (i - 1) %/% ncol + 1,
                   width = 1, fill = cols, label = as.character(i))
    })
    df <- do.call(rbind, rows)
    df$palette <- factor(df$palette, levels = unique(df$palette))
    # Index numbers in white on dark colours and in black on light ones
    rgb <- grDevices::col2rgb(df$fill)
    df$text <- ifelse(colSums(rgb * c(0.299, 0.587, 0.114)) < 140, "white", "black")

    ggplot(df, aes(x = .data$x, y = .data$y)) +
        geom_tile(aes(fill = I(.data$fill), width = .data$width), height = 0.9) +
        geom_text(aes(label = .data$label, colour = I(.data$text)), size = 3) +
        # One band per palette, as tall as its number of rows of swatches
        facet_grid(palette ~ ., scales = "free_y", space = "free_y", switch = "y") +
        scale_x_continuous(limits = c(0.5, ncol + 0.5), expand = expansion(0)) +
        scale_y_reverse(expand = expansion(add = 0.05)) +
        theme_void() +
        theme(
            strip.text.y.left = element_text(angle = 0, hjust = 1, margin = margin(r = 6)),
            panel.spacing = unit(4, "pt"),
            plot.background = element_rect(fill = "white", colour = NA),
            plot.margin = margin(6, 6, 6, 6)
        )
}


#' @title Check that colours stay distinct with colour vision deficiency
#'
#' @description Simulates deuteranopia and protanopia (red-green colour
#'              blindness, together about 8\% of men of European descent) and,
#'              optionally, tritanopia, with \pkg{colorspace}, and flags pairs
#'              of colours that become hard to tell apart. Useful after picking
#'              colours for cell types.
#'
#'              Distances are Euclidean distances in CIELAB (Delta E 1976).
#'              About 2.3 is just noticeable side by side; colours of small
#'              points or thin lines need far more to be told apart, hence the
#'              default of 10. Colours already close with normal vision are
#'              reported too.
#'
#' @param cols A palette name or a vector of colours
#'
#' @param min.distance Pairs closer than this are flagged, default is 10
#'
#' @param cvd Deficiencies to simulate, any of "deutan" and "protan"
#'            (default, both) and "tritan"
#'
#' @return A data frame of the flagged pairs (`colour1`, `colour2`, `vision`,
#'         `distance`), invisibly. A warning lists them when there are any
#'
#' @examples
#' \dontrun{
#' check_palette("okabe-ito")
#' check_palette(bt_palette("scanpy-20")[1:8])
#' }
#'
#' @export
#'
check_palette <- function(cols,
                          min.distance = 10,
                          cvd = c("deutan", "protan")) {
    check_installed("colorspace", reason = "to simulate colour vision deficiency")
    cvd <- match.arg(cvd, c("deutan", "protan", "tritan"), several.ok = TRUE)
    if (length(cols) == 1 && !is.null(find_palette(cols))) cols <- bt_palette(cols)
    bad <- !is_colour(cols)
    if (any(bad)) stop("Invalid colour(s): ", paste(cols[bad], collapse = ", "), call. = FALSE)
    hex <- grDevices::rgb(t(grDevices::col2rgb(cols)), maxColorValue = 255)

    simulate <- list(
        normal = identity,
        deutan = colorspace::deutan,
        protan = colorspace::protan,
        tritan = colorspace::tritan
    )[c("normal", cvd)]
    pairs <- utils::combn(length(hex), 2)
    res <- do.call(rbind, lapply(names(simulate), function(v) {
        lab <- colorspace::coords(methods::as(colorspace::hex2RGB(simulate[[v]](hex)), "LAB"))
        d <- sqrt(colSums((t(lab[pairs[1, ], , drop = FALSE]) - t(lab[pairs[2, ], , drop = FALSE]))^2))
        data.frame(colour1 = cols[pairs[1, ]], colour2 = cols[pairs[2, ]], vision = v, distance = d)
    }))
    res <- res[res$distance < min.distance, , drop = FALSE]
    rownames(res) <- NULL
    if (nrow(res)) {
        shown <- utils::head(res, 10)
        warning(nrow(res), " pair(s) of colours are hard to tell apart:\n",
                paste0("  ", shown$colour1, " / ", shown$colour2, " (", shown$vision,
                       ", Delta E ", round(shown$distance, 1), ")", collapse = "\n"),
                if (nrow(res) > 10) "\n  ...", call. = FALSE)
    }
    invisible(res)
}


# Colours used by the plot functions --------------------------------------------

#' Default colours for `n` groups
#'
#' "scanpy-20" for up to 20 groups, "scanpy-102" for up to 102, and evenly
#' spaced hues beyond.
#'
#' @param n Number of groups
#'
#' @return `n` colours
#'
#' @keywords internal
default_colours <- function(n) {
    if (n <= 20) return(bt_palette("scanpy-20", n))
    if (n <= 102) return(bt_palette("scanpy-102", n))
    message("More than 102 groups: using evenly spaced hues, which are hard to tell apart")
    scales::hue_pal()(n)
}


#' Colours for the levels of a grouping variable
#'
#' Every level gets a colour, including levels with no cells, so that a group
#' gets the same colour in a plot of a subset of cells as in a plot of all of
#' them.
#'
#' @param cols NULL for the default colours; a qualitative palette name;
#'             colours named by group; or colours in the order of the group
#'             levels
#'
#' @param group A factor
#'
#' @param arg Name of the argument, for messages
#'
#' @return A list with `group` and `cols`, named by group
#'
#' @keywords internal
group_colours <- function(cols,
                          group,
                          arg = "cols") {
    lv <- levels(group)
    if (length(cols) == 1 && !is.null(entry <- find_palette(cols))) {
        if (entry$type != "qualitative") {
            stop("`", arg, "` names a ", entry$type, " palette; give a qualitative one ",
                 "(see list_palettes(\"qualitative\"))", call. = FALSE)
        }
        cols <- bt_palette(cols, length(lv))
    }
    if (is.null(cols)) {
        cols <- stats::setNames(default_colours(length(lv)), lv)
    } else if (is.null(names(cols))) {
        if (length(cols) < length(lv)) {
            stop("`", arg, "` has ", length(cols), " colours but there are ",
                 length(lv), " groups", call. = FALSE)
        }
        cols <- stats::setNames(cols[seq_along(lv)], lv)
    } else if (length(missing_cols <- setdiff(lv, names(cols)))) {
        stop("No colour given for group(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
    }
    list(group = group, cols = cols)
}


#' Resolve the `cols` argument of a continuous colour scale into colours
#'
#' Shared by every continuous scale of the package (feature plots, dot plots,
#' heatmaps, density differences), so that they all accept the same values.
#'
#' @param cols NULL for the default (`default`, or "blue-magenta" for a
#'             diverging scale); a sequential or diverging palette name of the
#'             package, or an RColorBrewer palette name; or two or more
#'             colours from low to high, e.g. `c(low, high)` or
#'             `c(low, mid, high)`
#'
#' @param diverging Whether the scale is diverging, which sets the default
#'
#' @param default Default sequential palette, "plasma_r" unless the plot
#'                sets another
#'
#' @return A vector of colours, from low to high
#'
#' @keywords internal
resolve_continuous_palette <- function(cols,
                                       diverging,
                                       default = "plasma_r") {
    if (is.null(cols)) cols <- if (diverging) "blue-magenta" else default
    if (length(cols) == 1) {
        entry <- find_palette(cols)
        if (!is.null(entry)) {
            if (entry$type == "qualitative") {
                stop("Palette '", cols, "' is qualitative, for groups; give a sequential or ",
                     "diverging palette (see list_palettes())", call. = FALSE)
            }
            return(entry$colours)
        }
        # RColorBrewer names, also with "_r" for reversed (e.g. "BuPu_r")
        brewer <- sub("_r$", "", cols)
        if (requireNamespace("RColorBrewer", quietly = TRUE) &&
            brewer %in% rownames(RColorBrewer::brewer.pal.info)) {
            n <- RColorBrewer::brewer.pal.info[brewer, "maxcolors"]
            out <- RColorBrewer::brewer.pal(n, brewer)
            return(if (brewer != cols) rev(out) else out)
        }
        stop("Unknown palette '", cols, "'; give a palette name (see list_palettes()), ",
             "an RColorBrewer palette name, or two or more colours from low to high", call. = FALSE)
    }
    ok <- is_colour(cols)
    if (!all(ok)) stop("Invalid colour(s) in `cols`: ", paste(cols[!ok], collapse = ", "), call. = FALSE)
    unname(cols)
}
