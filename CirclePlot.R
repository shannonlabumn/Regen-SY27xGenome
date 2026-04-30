setwd("~/Desktop/Harpreet/")
library(circlize)
library(data.table)
library(ComplexHeatmap)
library(grid)
library(viridisLite)
library(RColorBrewer)
library(colorspace)
library(dplyr)


chrom <- fread("chrom.sizes", col.names = c("chr", "size"))
chrom$chr <- sub("(?i)^Chr", "", as.character(chrom$chr), perl = TRUE)

# Optional plotting/output
out_pdf <- "circlize_genome_32chr.pdf"
plot_width  <- 10
plot_height <- 10

## ---------- 2) PARAMETERS ----------
# Synteny filtering
min_mapq   <- 20
min_block  <- 10e3   # min aligned length on each side (10 kb)
max_links  <- 3000   # cap to avoid overplotting; decrease if too busy

# Track heights
ideogram_th   <- 0.06
gd_th         <- 0.08
repeat_th     <- 0.06

# Label size
chr_label_cex <- 0.55
axis_label    <- FALSE    # tick labels off (crowded genomes)
axis_cex      <- 0.35

## ---------- 3) READ & STANDARDIZE INPUTS ----------
# Helper to standardize column names case-insensitively
std_names <- function(dt, mapping) {
  ln <- tolower(names(dt))
  for (i in seq_along(mapping)) {
    w <- which(ln == tolower(names(mapping)[i]))
    if (length(w) == 1) setnames(dt, w, mapping[[i]])
  }
  dt
}

# --- Karyotype ---
karyo <- read.table("chrom.sizes", header = FALSE, sep="\t")
karyotype <- data.frame(Chr=karyo$V1, Start=rep(0, nrow(karyo)), End=karyo$V2)
kar <- karyotype
kar <- std_names(kar, c(Chr="Chr", Start="Start", End="End"))
stopifnot(all(c("Chr","Start","End") %in% names(kar)))
setorder(kar, Chr, Start)

# Chromosome lengths
setDT(kar)        
chr_len <- kar[, .(length = max(End, na.rm=TRUE)), by=Chr]
# If you want a specific order, relevel here (e.g., Chr1..Chr32)
# chr_len <- chr_len[order(as.integer(sub(".*?(\\d+)$","\\1", Chr)))]
sectors <- chr_len$Chr

# --- Gene density (per 100 kb) ---
gd <- fread("genes_density_100Kb.bed", col.names = c("Chr", "Start", "End", "Value"))
gd <- gd[, c(1:4)]
gd <- std_names(gd, c(Chr="Chr", Start="Start", End="End", Value="Value"))
stopifnot(all(c("Chr","Start","End","Value") %in% names(gd)))
gd[, `:=`(Start = as.numeric(Start), End = as.numeric(End), Value = as.numeric(Value))]
gd <- gd[Chr %in% sectors]
setorder(gd, Chr, Start)

# Repeats bed files
rep <- fread("LINE_100kb.cov")
rep <- rep[, -c(5,6,7)]
colnames(rep) = c("Chr", "Start", "End", "Value")
rep <- std_names(rep, c(Chr="Chr", Start="Start", End="End", Value="Value"))
stopifnot(all(c("Chr","Start","End","Value") %in% names(rep)))
rep[, `:=`(Start = as.numeric(Start), End = as.numeric(End), Value = as.numeric(Value))]
rep <- rep[Chr %in% sectors]
setorder(rep, Chr, Start)


# 12 reds from very light to deep red
reds12 <- colorRampPalette(brewer.pal(9, "Reds"))(12)

# Use quantile breaks to spread colors across your data
probs <- seq(0, 1, length.out = length(reds12))
brks  <- as.numeric(quantile(gd$Value, probs = probs, na.rm = TRUE))

# Ensure strictly increasing breaks (fallback to linear if many ties)
brks <- unique(brks)
if (length(brks) != length(reds12)) {
  rng  <- range(gd$Value, na.rm = TRUE)
  brks <- seq(rng[1], rng[2], length.out = length(reds12))
}

gd_col_fun <- circlize::colorRamp2(brks, reds12)

chr_col8 <- setNames(qualitative_hcl(8, palette = "Dark 3"),
                     paste0("Chr", 1:8))
circos.clear()
circos.par(start.degree= 90, gap.after = c(rep(2, nrow(chrom)-1), 8),
           track.margin = c(0.002, 0.002),
           cell.padding = c(0.00, 0.00, 0.00, 0.00))  # bigger gap after last chrom (degrees)

# Optional: draw outer track as "ideogram" ticks
circos.initialize(factors = sectors,
                  xlim = cbind(0, chr_len$length[match(sectors, chr_len$Chr)]))

## ---------- 6) OUTER CHR TRACK (labels + ticks) ----------
circos.track(
  ylim = c(0,1), track.height = ideogram_th, bg.border = NA,
  panel.fun = function(x, y) {
    si   <- CELL_META$sector.index        # e.g., "Chr3.2"
    xlim <- CELL_META$xlim
    
    root <- sub("\\..*$", "", si)         # "Chr3.2" -> "Chr3"
    fill_col <- chr_col8[root]            # repeat same color for .1.. .4
    
    circos.rect(xlim[1], 0, xlim[2], 1, col = fill_col, border = "#bdbdbd")
    
    # ticks & labels (keep whatever you had)
    # ticks + labels at ALL major ticks (Mb)
    maj  <- pretty(xlim, n = 4)
    maj  <- maj[maj >= xlim[1] & maj <= xlim[2]]        # keep inside sector
    labs <- sprintf("%.0f", maj/1e6)                    # or paste0(...," Mb")
    
    circos.axis(
      major.at = maj,
      labels   = labs,
      labels.cex = 0.45,
      minor.ticks = 4,
      major.tick.length = mm_y(1.6)
    )
    
    circos.text(mean(xlim), 1 + mm_y(2.2), si,
                cex = chr_label_cex, facing = "bending.inside", adj = c(0.5, -1.5))
  }
)



## --- Gene density heatmap as rectangles (robust) ---
# Assumes: sectors, chr_len, gd (Chr, Start, End, Value), gd_col_fun already defined

# 1) Prechecks & ordering
gd_df <- as.data.frame(gd[, .(Chr, Start, End, Value)])
gd_df$Chr   <- as.character(gd_df$Chr)
stopifnot(all(unique(gd_df$Chr) %in% sectors))
gd_df$Chr   <- factor(gd_df$Chr, levels = sectors)   # match sector order
gd_df$Start <- as.numeric(gd_df$Start)
gd_df$End   <- as.numeric(gd_df$End)
gd_df$Value <- as.numeric(gd_df$Value)

# 2) Zero padding & set a healthy track height
circos.par(cell.padding = c(0,0,0,0), track.margin = c(0.001, 0.001))
gd_th <- 0.12   # increase to 0.15–0.18 if you want it bolder

circos.genomicTrackPlotRegion(
  gd_df,
  ylim = c(0, 1),
  track.height = gd_th,
  bg.border = NA,
  panel.fun = function(region, value, ...) {
    bed  <- cbind(region, value)  # Chr Start End Value
    lo   <- brks[1]; hi <- brks[length(brks)]
    vals <- pmin(hi, pmax(lo, bed$Value))   # clip to scale
    cols <- gd_col_fun(vals)
    cols[!is.finite(bed$Value)] <- "white" # show missing as light gray
    circos.genomicRect(bed, ybottom = 0, ytop = 1, col = cols, border = NA)
  }
)


circos.par(cell.padding = c(0,0,0,0), track.margin = c(0.050, 0.002))

hi <- as.numeric(quantile(rep$Value, 0.99, na.rm=TRUE))

circos.par(track.margin = c(0.04, 0.002))
circos.genomicTrackPlotRegion(
  rep, ylim = c(0, hi), track.height = 0.15, bg.border = NA,
  panel.fun = function(region, value, ...) {
    bed <- cbind(region, value)  # has Value
    circos.genomicRect(
      bed,
      ybottom = 0,
      ytop    = pmin(bed$Value, hi),   # bars up to clipped max
      col     = "#191970",             # or "midnightblue" etc.
      border  = NA
    )
  }
)

# Repeats bed files
rep <- fread("LTR_100kb.cov")
rep <- rep[, -c(5,6,7)]
colnames(rep) = c("Chr", "Start", "End", "Value")
rep <- std_names(rep, c(Chr="Chr", Start="Start", End="End", Value="Value"))
stopifnot(all(c("Chr","Start","End","Value") %in% names(rep)))
rep[, `:=`(Start = as.numeric(Start), End = as.numeric(End), Value = as.numeric(Value))]
rep <- rep[Chr %in% sectors]
setorder(rep, Chr, Start)

circos.par(track.margin = c(0.050, 0.002))

hi <- as.numeric(quantile(rep$Value, 0.99, na.rm=TRUE))

circos.par(track.margin = c(0.04, 0.002))
circos.genomicTrackPlotRegion(
  rep, ylim = c(0, hi), track.height = 0.15, bg.border = NA,
  panel.fun = function(region, value, ...) {
    bed <- cbind(region, value)  # has Value
    circos.genomicRect(
      bed,
      ybottom = 0,
      ytop    = pmin(bed$Value, hi),   # bars up to clipped max
      col     = "#8B0000",             # or "midnightblue" etc.
      border  = NA
    )
  }
)
# Repeats bed files
rep <- fread("repeat_total_100Kb.cov")
rep <- rep[, -c(5,6,7)]
colnames(rep) = c("Chr", "Start", "End", "Value")
rep <- std_names(rep, c(Chr="Chr", Start="Start", End="End", Value="Value"))
stopifnot(all(c("Chr","Start","End","Value") %in% names(rep)))
rep[, `:=`(Start = as.numeric(Start), End = as.numeric(End), Value = as.numeric(Value))]
rep <- rep[Chr %in% sectors]
setorder(rep, Chr, Start)

circos.par(track.margin = c(0.050, 0.002))

hi <- as.numeric(quantile(rep$Value, 0.99, na.rm=TRUE))

circos.par(track.margin = c(0.04, 0.002))
circos.genomicTrackPlotRegion(
  rep, ylim = c(0, hi), track.height = 0.15, bg.border = NA,
  panel.fun = function(region, value, ...) {
    bed <- cbind(region, value)  # has Value
    circos.genomicRect(
      bed,
      ybottom = 0,
      ytop    = pmin(bed$Value, hi),   # bars up to clipped max
      col     = "darkgreen",             # or "midnightblue" etc.
      border  = NA
    )
  }
)

## ---------- SYNTENY LINKS FROM BEDPE (Chr1 Start1 End1 Chr2 Start2 End2) ----------

bedpe_file <- "synteny_all.bedpe"  # << set your path
min_block  <- 1e2                 # keep links ≥10 kb on both sides
max_links  <- 10000                 # downsample if too dense
alpha_link <- 0.1                  # transparency 0..1

# 1) Read
bp <- fread(bedpe_file)
setnames(bp, c("Chr1","Start1","End1","Chr2","Start2","End2"))  # matches your header
bp[, c("Start1","End1","Start2","End2") := lapply(.SD, as.numeric),
   .SDcols = c("Start1","End1","Start2","End2")]
bp[, `:=`(Chr1 = as.character(Chr1), Chr2 = as.character(Chr2))]

# 2) Keep only chromosomes present in the circle & valid intervals
bp <- bp[Chr1 %in% sectors & Chr2 %in% sectors]
bp <- bp[End1 > Start1 & End2 > Start2]

# 3) Clamp to chromosome lengths (avoids warnings if anything runs out-of-bounds)
chr2len <- setNames(chr_len$length, chr_len$Chr)
bp[, Start1 := pmax(0, pmin(Start1, chr2len[Chr1]))]
bp[, End1   := pmax(0, pmin(End1,   chr2len[Chr1]))]
bp[, Start2 := pmax(0, pmin(Start2, chr2len[Chr2]))]
bp[, End2   := pmax(0, pmin(End2,   chr2len[Chr2]))]
bp <- bp[End1 > Start1 & End2 > Start2]

# # 4) Filter tiny blocks
# bp <- bp[(End1 - Start1) >= min_block & (End2 - Start2) >= min_block]
# 
# # 5) De-duplicate symmetric pairs (A↔B) so each link is drawn once
# swap <- bp$Chr1 > bp$Chr2 | (bp$Chr1 == bp$Chr2 & bp$Start1 > bp$Start2) |
#   (bp$Chr1 == bp$Chr2 & bp$Start1 == bp$Start2 & bp$End1 > bp$End2)
# bp[swap, c("Chr1","Start1","End1","Chr2","Start2","End2") :=
#      .(Chr2, Start2, End2, Chr1, Start1, End1)]
# bp <- unique(bp, by = c("Chr1","Start1","End1","Chr2","Start2","End2"))

# 6) Coloring: length-based (log scale). Inter- vs intra-chrom get different hues.
len_min <- pmin(bp$End1 - bp$Start1, bp$End2 - bp$Start2)
is_intra <- bp$Chr1 == bp$Chr2

# two palettes: blue for inter-chrom, green for intra-chrom (adjust if you like)
len_brks  <- quantile(log1p(len_min), probs = c(0, .5, 1), na.rm = TRUE)
col_fun_I <- colorRamp2(len_brks, c("#dbe9f6", "#7bb6d6", "#1f77b4"))  # inter
col_fun_S <- colorRamp2(len_brks, c("#e2f1e8", "#7bc87c", "#2ca02c"))  # intra

link_cols <- ifelse(is_intra,
                    col_fun_S(log1p(len_min)),
                    col_fun_I(log1p(len_min)))
link_cols <- add_transparency(link_cols, alpha_link)

# # 7) Optional: only inter-chrom links
# bp <- bp[Chr1 != Chr2]; link_cols <- link_cols[bp$Chr1 != bp$Chr2]
# 
# # # 8) Downsample if huge
# if (nrow(bp) > max_links) {
#    set.seed(1)
#    keep <- sample.int(nrow(bp), max_links)
#    bp <- bp[keep]
#    link_cols <- link_cols[keep]
#  }

# 9) Draw (after you’ve added the inner tracks)
df1 <- data.frame(chr=factor(bp$Chr1, levels=sectors), start=bp$Start1, end=bp$End1)
df2 <- data.frame(chr=factor(bp$Chr2, levels=sectors), start=bp$Start2, end=bp$End2)

circos.genomicLink(df1, df2, col = link_cols, border = NA, lwd = 2)

################################################################
###############################################################

bedpe_file <- "synteny_all.bedpe"

# 0) (Optional) render to a fast raster device instead of PDF with alpha
# install.packages("ragg")
# ragg::agg_png("circos_synteny.png", width = 3000, height = 3000, res = 360)

bp <- fread(bedpe_file)
setnames(bp, c("Chr1","Start1","End1","Chr2","Start2","End2"))
bp[, c("Start1","End1","Start2","End2") := lapply(.SD, as.numeric),
   .SDcols = c("Start1","End1","Start2","End2")]
bp <- bp[Chr1 %in% sectors & Chr2 %in% sectors &
           End1 > Start1 & End2 > Start2]

# Clamp to chromosome bounds
chr2len <- setNames(chr_len$length, chr_len$Chr)
bp[, Start1 := pmax(0, pmin(Start1, chr2len[Chr1]))]
bp[, End1   := pmax(0, pmin(End1,   chr2len[Chr1]))]
bp[, Start2 := pmax(0, pmin(Start2, chr2len[Chr2]))]
bp[, End2   := pmax(0, pmin(End2,   chr2len[Chr2]))]
bp <- bp[End1 > Start1 & End2 > Start2]

# 1) Filter tiny blocks HARD (you had 1e2 = 100 bp; use ≥ 5e4 = 50 kb or ≥ 1e4 = 10 kb)
min_block <- 1e4
bp <- bp[(End1-Start1) >= min_block & (End2-Start2) >= min_block]

# 2) De-duplicate symmetric pairs (A↔B) so each link is drawn once
swap <- bp$Chr1 > bp$Chr2 | (bp$Chr1 == bp$Chr2 & bp$Start1 > bp$Start2) |
  (bp$Chr1 == bp$Chr2 & bp$Start1 == bp$Start2 & bp$End1 > bp$End2)
bp[swap, c("Chr1","Start1","End1","Chr2","Start2","End2") :=
     .(Chr2, Start2, End2, Chr1, Start1, End1)]
bp <- unique(bp, by = c("Chr1","Start1","End1","Chr2","Start2","End2"))

# 3) Keep the strongest signal: longest blocks per Chr1–Chr2 pair
bp[, len := pmin(End1-Start1, End2-Start2)]
bp[, pair := paste(Chr1, Chr2, sep="|")]

# cap per pair and overall (tune these)
k_per_pair <- 200
bp <- bp[order(-len), .SD[1:min(.N, k_per_pair)], by = pair]

max_links <- 10000
if (nrow(bp) > max_links) bp <- bp[order(-len)][1:max_links]

# 4) Color per chromosome (Chr1) and draw with thin lines, modest alpha
alpha_link <- 0.35     # heavy alpha on PDF is slow; modest alpha is fine
lwd_link   <- 0.4

# reuse the per-chromosome palette defined earlier:
chr_col <- setNames(qualitative_hcl(length(sectors), "Dark 3"), sectors)
link_cols <- add_transparency(chr_col[bp$Chr1], alpha_link)

df1 <- data.frame(chr=factor(bp$Chr1, levels=sectors), start=bp$Start1, end=bp$End1)
df2 <- data.frame(chr=factor(bp$Chr2, levels=sectors), start=bp$Start2, end=bp$End2)

circos.genomicLink(df1, df2, col = link_cols, border = NA, lwd = lwd_link)




