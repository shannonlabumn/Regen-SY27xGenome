#!/usr/bin/env Rscript
# 02_te_proximity_with_genes.R
# Measure the nearest repeat and Helitron distance for each annotated gene.
# Summarize all genes and the background excluding R-gene-overlapping loci.
# Optionally compare this background with an existing R-gene distance table.
# Distances with no target on the same sequence are NA, not zero or infinity.
# Requirements: data.table, GenomicRanges, IRanges, and Unix grep.
# Run from the project directory; relative paths resolve there.
# Usage: Rscript 02_te_proximity_with_genes.R

project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
if (!dir.exists(project_dir)) stop("PROJECT_DIR not found: ", project_dir)
setwd(project_dir)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
})

# STEP 1: Configure annotation paths for one consistent assembly.
# Scripts 01/02 and 03 used different general GFFs in the supplied analyses;
# separate defaults preserve that distinction. Set both to the same file only
# if that matches the intended analysis.
GENOME_GFF    <- "data/annotations/general_genes.gff3"
BRAKER_GFF    <- "data/annotations/rgene_models.gff3"
RGENE_CLASSES <- "data/annotations/rgene_classes.tsv"      
OUTDIR        <- "results/repeat_gene_association"
# Cache reduced repeat intervals; use a fresh OUTDIR if annotations change.
MASK_ANY_RDS  <- file.path(OUTDIR, "any_repeat_mask.rds")
# Helitron mask contains reduced DNA/Helitron intervals (strand ignored).
MASK_HEL_RDS  <- file.path(OUTDIR, "helitron_mask.rds")
TE_GFF        <- "data/repeats/EDTA.intact.gff3"
OUT_CSV       <- file.path(OUTDIR, "per_gene_genome_te_proximity.csv")
OUT_SUMMARY   <- file.path(OUTDIR, "genome_te_proximity_summary.txt")

# Optional externally generated R-gene distances; see README for columns.
RGENE_TE_PROXIMITY <- file.path("data/rgene_comparisons/per_gene_te_proximity.csv")

# Cumulative nearest-distance thresholds, in bp.
THRESHOLDS <- c(1000, 5000, 10000, 20000)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# STEP 2: Reuse cached intervals or build masks from EDTA classifications.
# reduce() merges overlapping and directly adjacent intervals, so resulting
# Helitron loci need not correspond one-to-one with original EDTA records.
if (file.exists(MASK_ANY_RDS) && file.exists(MASK_HEL_RDS)) {
  message("[1/5] Loading cached repeat masks...")
  any_repeat_mask <- readRDS(MASK_ANY_RDS)
  helitron_mask <- readRDS(MASK_HEL_RDS)
} else {
  message("[1/5] No cached masks found -- rebuilding from EDTA GFF3 (slow step)...")
  te_dt <- fread(
    cmd = paste("grep -v '^#'", shQuote(TE_GFF)),
    header = FALSE, sep = "\t", quote = "",
    col.names = c("seqid", "source", "type", "start", "end", "score", "strand",
                  "phase", "attributes")
  )
  te_dt[, classification := sub(".*classification=([^;]+).*", "\\1", attributes)]
  te_gr <- GRanges(seqnames = te_dt$seqid,
                   ranges = IRanges(start = te_dt$start, end = te_dt$end),
                   classification = te_dt$classification)
  any_repeat_mask <- reduce(te_gr, ignore.strand = TRUE)
  helitron_mask <- reduce(te_gr[te_gr$classification == "DNA/Helitron"],
                          ignore.strand = TRUE)
  saveRDS(any_repeat_mask, MASK_ANY_RDS)
  saveRDS(helitron_mask, MASK_HEL_RDS)
}

message("[2/5] Getting whole-gene coordinates for the complete annotation...")
# STEP 3: Read nine-column GFF3 annotations; coordinates are one-based inclusive.
# GFF3 inputs must contain feature rows only after comments (no FASTA section).
gff <- fread(
  cmd = paste("grep -v '^#'", shQuote(GENOME_GFF)),
  header = FALSE, sep = "\t", quote = "",
  col.names = c("seqid", "source", "type", "start", "end", "score", "strand",
                "phase", "attributes")
)
get_attr <- function(x, key) sub(paste0(".*", key, "=([^;]+).*"), "\\1", x)

# Select gene spans; remove duplicate records before counting genes.
genes <- gff[type == "gene"]
genes[, gene_id := get_attr(attributes, "ID")]

n_before_dedup <- nrow(genes)
genes <- unique(genes, by = c("gene_id", "seqid", "start", "end"))
n_dropped <- n_before_dedup - nrow(genes)
if (n_dropped > 0) {
  message(sprintf("  NOTE: dropped %d exact-duplicate 'gene' rows.", n_dropped))
}
# If a gene ID has multiple coordinate sets, retain the first and warn.
# Review these warnings before interpreting the outputs.
still_dup <- genes$gene_id[duplicated(genes$gene_id)]
if (length(still_dup) > 0) {
  warning(sprintf(
    "%d gene_id(s) map to >1 DIFFERENT coordinate set (e.g. %s) -- keeping only the first occurrence. Investigate before trusting those genes' distances.",
    length(unique(still_dup)), paste(head(unique(still_dup), 5), collapse = ", ")
  ))
  genes <- unique(genes, by = "gene_id")
}

genes_gr <- GRanges(seqnames = genes$seqid,
                    ranges = IRanges(start = genes$start, end = genes$end),
                    gene_id = genes$gene_id)
message(sprintf("  %d unique genes with coordinates.", length(genes_gr)))

message("[2b/5] Flagging general-annotation genes that overlap an R-gene locus...")
# Restrict R-gene models to IDs in column 1 of the headerless class table.
# Their sequence names must match the repeat and general-gene annotations.
rgene_target_ids <- unique(fread(RGENE_CLASSES, header = FALSE)$V1)
braker_gff <- fread(
  cmd = paste("grep -v '^#'", shQuote(BRAKER_GFF)),
  header = FALSE, sep = "\t", quote = "",
  col.names = c("seqid", "source", "type", "start", "end", "score", "strand",
                "phase", "attributes")
)
rgenes_raw <- braker_gff[type == "gene"]
rgenes_raw[, gene_id := get_attr(attributes, "ID")]
rgenes_raw <- rgenes_raw[gene_id %in% rgene_target_ids]
rgenes_raw <- unique(rgenes_raw, by = c("gene_id", "seqid", "start", "end"))
rgenes_raw <- unique(rgenes_raw, by = "gene_id")
rgene_loci_gr <- GRanges(seqnames = rgenes_raw$seqid,
                         ranges = IRanges(start = rgenes_raw$start, end = rgenes_raw$end))

# Flag general genes overlapping an R-gene locus by genomic coordinates.
ov_rgene <- findOverlaps(genes_gr, rgene_loci_gr, ignore.strand = TRUE)
overlaps_rgene_locus <- logical(length(genes_gr))
overlaps_rgene_locus[unique(queryHits(ov_rgene))] <- TRUE
message(sprintf("  %d general-annotation genes overlap at least one of the %d R-gene loci.",
                sum(overlaps_rgene_locus), length(rgene_loci_gr)))

message("[3/5] Computing distance to nearest repeat (0 = gene overlaps a repeat)...")
# STEP 4: Compute strand-independent nearest distances on the same sequence.
# Initialize with NA to preserve queries without any eligible nearest target.
# Zero distance includes overlap and direct adjacency.
dist_to_nearest <- function(query_gr, mask_gr) {
  result <- rep(NA_real_, length(query_gr))
  hits <- distanceToNearest(query_gr, mask_gr, ignore.strand = TRUE)
  if (length(hits) > 0) {
    result[queryHits(hits)] <- mcols(hits)$distance
  }
  result
}

dist_any <- dist_to_nearest(genes_gr, any_repeat_mask)
dist_hel <- dist_to_nearest(genes_gr, helitron_mask)

# STEP 5: Assemble per-gene metrics, preserving unmeasured values as NA.
per_gene <- data.table(
  gene_id = genes_gr$gene_id,
  chrom = as.character(seqnames(genes_gr)),
  start = start(genes_gr),
  end = end(genes_gr),
  overlaps_rgene_locus = overlaps_rgene_locus,
  dist_to_any_repeat = dist_any,
  dist_to_helitron = dist_hel
)

n_unmeasured <- sum(is.na(per_gene$dist_to_any_repeat))
if (n_unmeasured > 0) {
  message(sprintf("  NOTE: %d genes on scaffolds with no repeat annotation -- distance left NA.",
                  n_unmeasured))
}

fwrite(per_gene, OUT_CSV)

message("[4/5] Writing proximity summary...")
lines <- c()
lines <- c(lines, sprintf("Total genes in general annotation (%s): %d", GENOME_GFF, nrow(per_gene)))
lines <- c(lines, "")

# Report measurable denominators separately; NA values do not enter percentages.
report_thresholds <- function(dist_vec, label) {
  n_total <- length(dist_vec)
  ok <- !is.na(dist_vec)
  n_na <- sum(!ok)
  dist_vec <- dist_vec[ok]
  n <- length(dist_vec)
  out <- c(sprintf("--- %s ---", label))
  out <- c(out, sprintf("Genes with a measurable distance: %d / %d (excluded: %d)", n, n_total, n_na))
  if (n == 0) return(c(out, "  (no genes to report on -- skipping)"))
  out <- c(out, sprintf("Zero distance (overlapping or adjacent): %d (%.1f%%)",
                        sum(dist_vec == 0), 100 * mean(dist_vec == 0)))
  for (t in THRESHOLDS) {
    out <- c(out, sprintf("Within %d bp (cumulative): %d (%.1f%%)", t, sum(dist_vec <= t), 100 * mean(dist_vec <= t)))
  }
  out
}

lines <- c(lines, "==== ALL GENES, GENERAL ANNOTATION (general annotation, INCLUDES R-gene loci) ====")
lines <- c(lines, report_thresholds(per_gene$dist_to_any_repeat, "Any repeat"))
lines <- c(lines, "")
lines <- c(lines, report_thresholds(per_gene$dist_to_helitron, "Helitron"))
lines <- c(lines, "")

# Exclude R-gene-overlapping general models for the comparison background.
non_rgene <- per_gene[overlaps_rgene_locus == FALSE]
lines <- c(lines, "==== NON-R GENES ONLY (R-gene-overlapping loci excluded -- this is the correct background) ====")
lines <- c(lines, report_thresholds(non_rgene$dist_to_any_repeat, "Any repeat"))
lines <- c(lines, "")
lines <- c(lines, report_thresholds(non_rgene$dist_to_helitron, "Helitron"))
lines <- c(lines, "")

# Optional comparison: summarize external R-gene metrics and run Wilcoxon.
# These are unadjusted tests; spatial dependence is not modeled.
if (file.exists(RGENE_TE_PROXIMITY)) {
  rgene_prox <- fread(RGENE_TE_PROXIMITY)
  lines <- c(lines, sprintf("==== R GENES (canonical set, from %s, n=%d) ====",
                            RGENE_TE_PROXIMITY, nrow(rgene_prox)))
  lines <- c(lines, report_thresholds(rgene_prox$dist_to_any_repeat, "Any repeat"))
  lines <- c(lines, "")
  lines <- c(lines, report_thresholds(rgene_prox$dist_to_helitron, "Helitron"))
  lines <- c(lines, "")

  nonrgene_any <- non_rgene$dist_to_any_repeat[!is.na(non_rgene$dist_to_any_repeat)]
  rgene_any <- rgene_prox$dist_to_any_repeat[!is.na(rgene_prox$dist_to_any_repeat)]
  wt_any <- wilcox.test(rgene_any, nonrgene_any)
  lines <- c(lines, sprintf("Wilcoxon (distance to any repeat, R genes vs. non-R genes ONLY): W=%.3e, p=%.3e",
                            wt_any$statistic, wt_any$p.value))
  lines <- c(lines, sprintf("  medians: R genes=%.0f bp, non-R genes=%.0f bp",
                            median(rgene_any), median(nonrgene_any)))

  nonrgene_hel <- non_rgene$dist_to_helitron[!is.na(non_rgene$dist_to_helitron)]
  rgene_hel <- rgene_prox$dist_to_helitron[!is.na(rgene_prox$dist_to_helitron)]
  wt_hel <- wilcox.test(rgene_hel, nonrgene_hel)
  lines <- c(lines, sprintf("Wilcoxon (distance to Helitron, R genes vs. non-R genes ONLY): W=%.3e, p=%.3e",
                            wt_hel$statistic, wt_hel$p.value))
  lines <- c(lines, sprintf("  medians: R genes=%.0f bp, non-R genes=%.0f bp",
                            median(rgene_hel), median(nonrgene_hel)))
  lines <- c(lines, "(A smaller median distance = closer to a repeat -- check which group's median is smaller before interpreting the p-value's meaning.)")
} else {
  lines <- c(lines, sprintf(
    "NOTE: %s not found -- skipping the R-gene-vs-genome proximity comparison. Provide the optional R-gene proximity table; see README.md.",
    RGENE_TE_PROXIMITY))
}

# FINAL STEP: Save a human-readable summary and print it to the console.
writeLines(lines, OUT_SUMMARY)
writeLines(lines)
cat("\nWrote per-gene distances to", OUT_CSV, "\n")
cat("Wrote summary to", OUT_SUMMARY, "\n")
