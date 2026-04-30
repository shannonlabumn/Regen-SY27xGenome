setwd("~/Desktop/Harpreet/")
## ==========================================================
## R-gene ideograms (large PNGs) using RIdeogram
## ==========================================================

## ---------- CONFIG ----------
chrom_sizes_path <- "chrom.sizes"
gff_files <- c(
  "RSY_shortid_braker_aa_CNB.gff3",
  "RSY_shortid_braker_aa_CNL.gff3",
  "RSY_shortid_braker_aa_NBARC.gff3",
  "RSY_shortid_braker_aa_NLR.gff3",
  "RSY_shortid_braker_aa_RNB.gff3",
  "RSY_shortid_braker_aa_RNL.gff3",
  "RSY_shortid_braker_aa_RxNB.gff3",
  "RSY_shortid_braker_aa_RxNL.gff3",
  "RSY_shortid_braker_aa_TNB.gff3",
  "RSY_shortid_braker_aa_TNL.gff3"
)

window_bp <- 100000L          # bin size (use 250000L-500000L for thicker bars)
normalize_per_chr <- TRUE      # scale counts to 0..1 per chromosome for visual comparability

## Figure size / resolution (PNG)
png_width_in  <- 14
png_height_in <- 10
png_dpi       <- 600
keep_svg      <- FALSE         # set TRUE if you also want to keep SVGs

## Ideogram layout (vignette params)
ideo_width <- 280              # widen drawing region (vignette default is ~170)
legend_x   <- ideo_width - 10  # Legend position X
legend_y   <- 35               # Legend position Y

## ---------- PACKAGES ----------
if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table")
if (!requireNamespace("RIdeogram",  quietly = TRUE)) install.packages("RIdeogram")
if (!requireNamespace("rtracklayer", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("rtracklayer")
}
# helpers for conversion / panel
if (!requireNamespace("rsvg", quietly = TRUE))   try(install.packages("rsvg"),   silent = TRUE)
if (!requireNamespace("magick", quietly = TRUE)) try(install.packages("magick"), silent = TRUE)

library(data.table)
library(rtracklayer)
library(RIdeogram)

## ---------- Helper: SVG -> big PNG (uses convertSVG per vignette; falls back to rsvg) ----------
save_png <- function(svg_path, width_in = png_width_in, height_in = png_height_in, dpi = png_dpi, keep_svg = FALSE) {
  base <- sub("\\.svg$", "", basename(svg_path))
  out_dir <- dirname(svg_path)
  png_path <- file.path(out_dir, paste0(base, ".png"))
  ok <- FALSE
  
  # Prefer the vignette's converter (lets us set size/dpi)
  try({
    RIdeogram::convertSVG(svg = svg_path, file = file.path(out_dir, base),
                          device = "png", width = width_in, height = height_in, dpi = dpi)
    ok <- file.exists(png_path)
  }, silent = TRUE)
  
  # Fallback: rsvg (pixels = inches * dpi)
  if (!ok && requireNamespace("rsvg", quietly = TRUE)) {
    try({
      rsvg::rsvg_png(svg_path, png_path,
                     width  = as.integer(width_in  * dpi),
                     height = as.integer(height_in * dpi))
      ok <- file.exists(png_path)
    }, silent = TRUE)
  }
  
  if (!ok) warning("PNG conversion failed: ", svg_path)
  if (ok && !keep_svg) try(unlink(svg_path), silent = TRUE)
  invisible(png_path)
}

## ---------- CHROM ORDER ----------
chrom_order <- sprintf("Chr%s.%s", rep(1:8, each = 4), rep(1:4, 8))

## ---------- READ chrom.sizes & KARYOTYPE (3-column style is fine per vignette) ----------
cs <- fread(chrom_sizes_path, col.names = c("chr", "size"))
cs <- cs[chr %in% chrom_order]
cs[, chr := factor(chr, levels = chrom_order)]
setorder(cs, chr)

# RIdeogram accepts 3 or 5 columns; we keep 4 with a Name column for labels
karyo <- data.frame(
  Chr   = as.character(cs$chr),
  Start = as.integer(0L),
  End   = as.integer(cs$size),
  Name  = as.character(cs$chr),
  check.names = FALSE
)
stopifnot(nrow(karyo) == 32L)

## ---------- BUILD FIXED BINS ----------
bins <- rbindlist(lapply(seq_len(nrow(cs)), function(i) {
  chr <- as.character(cs$chr[i]); L <- as.integer(cs$size[i])
  starts <- as.integer(seq(0L, max(0L, L - 1L), by = window_bp))
  data.table(chr = chr, start = starts, end = pmin(starts + window_bp, L))
}))

## ---------- GFF3 -> gene spans (robust) ----------
read_genes_from_gff <- function(path) {
  gr  <- import(path)  # GRanges (1-based closed)
  typ <- as.character(mcols(gr)$type)
  if ("gene" %in% typ) {
    gr <- gr[typ == "gene"]
    dt <- data.table(chr = as.character(seqnames(gr)),
                     start = as.integer(start(gr)) - 1L,
                     end   = as.integer(end(gr)))
  } else if ("mRNA" %in% typ) {
    gr <- gr[typ == "mRNA"]
    dt <- data.table(chr = as.character(seqnames(gr)),
                     start = as.integer(start(gr)) - 1L,
                     end   = as.integer(end(gr)))
  } else if (!is.null(mcols(gr)$Parent)) {
    tmp <- data.table(chr   = as.character(seqnames(gr)),
                      start = as.integer(start(gr)) - 1L,
                      end   = as.integer(end(gr)),
                      parent= as.character(mcols(gr)$Parent))
    dt <- tmp[, .(start = min(start), end = max(end)), by = .(chr, parent)][, parent := NULL][]
  } else {
    dt <- data.table(chr = as.character(seqnames(gr)),
                     start = as.integer(start(gr)) - 1L,
                     end   = as.integer(end(gr)))
  }
  dt <- dt[chr %in% chrom_order & end > start]
  dt[]
}

message("Reading R-gene GFF3s ...")
r_gene_lists <- lapply(gff_files, read_genes_from_gff)
track_names  <- sub("\\.gff3$", "", basename(gff_files))
track_names  <- sub("^RSY_shortid_braker_aa_", "", track_names)
names(r_gene_lists) <- track_names

## ---------- Bin counts (exact overlap via foverlaps) ----------
count_overlaps <- function(features_dt, bins_dt) {
  if (nrow(features_dt) == 0L) {
    out <- copy(bins_dt); out[, value := 0L]; return(out[])
  }
  b <- copy(bins_dt); f <- copy(features_dt)
  setkey(b, chr, start, end); setkey(f, chr, start, end)
  ol <- foverlaps(f, b, by.x = c("chr","start","end"), by.y = c("chr","start","end"), nomatch = 0L)
  dens <- ol[, .N, by = .(chr, start, end)]
  out <- merge(b, dens, by = c("chr","start","end"), all.x = TRUE)
  out[is.na(N), N := 0L]
  setnames(out, "N", "value")
  out[]
}

dens_total <- count_overlaps(rbindlist(r_gene_lists, use.names = TRUE, fill = TRUE), bins)
dens_list  <- lapply(r_gene_lists, count_overlaps, bins_dt = bins)

if (normalize_per_chr) {
  dens_total[, value := ifelse(max(value) > 0, value / max(value), 0), by = chr]
  dens_list <- lapply(dens_list, function(d) { d[, value := ifelse(max(value) > 0, value / max(value), 0), by = chr]; d })
}

## ---------- Rename to RIdeogram's required columns ----------
to_rideo <- function(d) {
  setDT(d)
  setnames(d, c("chr","start","end","value"), c("Chr","Start","End","Value"), skip_absent = TRUE)
  d[, `:=`(Chr = as.character(Chr), Start = as.integer(Start), End = as.integer(End))]
  d[Chr %in% karyo$Chr][]
}
dens_total <- to_rideo(dens_total)
dens_list  <- lapply(dens_list, to_rideo)

## ---------- Palettes (overlaid / label heatmap) ----------
pal_total <- colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(100)
pal_label <- colorRampPalette(c("#ffffe5", "#fdae6b", "#e6550d"))(100)

## ==========================================================
## FIGURE 1 — Total R-gene density (big, crisp PNG)
## ==========================================================
message("Drawing total-density ideogram ...")
svg1 <- "Rgenes_total_density.svg"
ideogram(
  karyotype = karyo,          # (Chr, Start, End, Name)
  overlaid  = dens_total,     # (Chr, Start, End, Value)
  colorset1 = pal_total,
  width     = ideo_width,     # wider region per vignette
  Lx        = legend_x,
  Ly        = legend_y,
  output    = svg1
)
save_png(svg1, keep_svg = keep_svg)

## ==========================================================
## FIGURES 2..11 — Total + one family as label heatmap (10 PNGs)
## ==========================================================
message("Drawing family heatmaps ...")
for (nm in names(dens_list)) {
  lab <- dens_list[[nm]]
  svg_out <- paste0("Rgenes_total_plus_", nm, ".svg")
  ideogram(
    karyotype = karyo,
    overlaid  = dens_total,
    colorset1 = pal_total,
    label     = lab,
    label_type= "heatmap",
    colorset2 = pal_label,
    width     = ideo_width,
    Lx        = legend_x,
    Ly        = legend_y,
    output    = svg_out
  )
  save_png(svg_out, keep_svg = keep_svg)
}

## ==========================================================
## Optional: make a 2×5 panel of the 10 family PNGs
## ==========================================================
panel_out <- "Rgenes_families_panel_2x5.png"
if (requireNamespace("magick", quietly = TRUE)) {
  message("Assembling 2x5 panel ...")
  family_pngs <- file.path(getwd(), paste0("Rgenes_total_plus_", names(dens_list), ".png"))
  if (all(file.exists(family_pngs))) {
    imgs <- magick::image_read(family_pngs)
    # simple montage; adjust geometry if you want larger margins
    panel <- magick::image_montage(imgs, tile = "5x2", geometry = "+10+10")
    magick::image_write(panel, panel_out, density = paste0(png_dpi, "x", png_dpi))
  } else {
    message("Panel skipped: not all family PNGs were found.")
  }
} else {
  message("Install 'magick' to auto-build a 2x5 panel.")
}

message("Done. PNGs are in: ", getwd(),
        "\n - Rgenes_total_density.png",
        "\n - Rgenes_total_plus_<FAMILY>.png (10 files)",
        "\n - Optional panel: ", panel_out)

