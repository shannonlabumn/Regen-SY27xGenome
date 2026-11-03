#!/usr/bin/env Rscript
# 03_helitron_gene_proximity.R
# Measure the nearest general gene and canonical R gene for each reduced
# Helitron interval. The denominator is Helitron loci, not genes.
# MIN_CDS_LENGTH=100 retains the supplied general-gene filter. Its current
# implementation sums CDS rows across transcripts; it is not a longest-
# isoform or non-redundant CDS-length filter.
# The gene-span coverage ratio is descriptive, not a permutation-based null.
# Requirements: data.table, GenomicRanges, IRanges, and Unix grep.
# Run from the project directory; relative paths resolve there.
# Usage: Rscript 03_helitron_gene_proximity.R

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
GENOME_GFF     <- "data/annotations/general_genes_helitron_analysis.gff3"
BRAKER_GFF     <- "data/annotations/rgene_models.gff3"
RGENE_CLASSES  <- "data/annotations/rgene_classes.tsv"
OUTDIR         <- "results/repeat_gene_association"
# Helitron mask contains reduced DNA/Helitron intervals (strand ignored).
MASK_HEL_RDS   <- file.path(OUTDIR, "helitron_mask.rds")
TE_GFF         <- "data/repeats/EDTA.intact.gff3"
OUT_CSV        <- file.path(OUTDIR, "per_helitron_gene_proximity.csv")
OUT_SUMMARY    <- file.path(OUTDIR, "helitron_gene_proximity_summary.txt")

# Retain genes with summed CDS-row length >=100 bp; NULL disables the filter.
MIN_CDS_LENGTH <- 100

# Supply the total length of the SAME assembly used for every annotation.
# Example: sum column 2 of the genome FASTA .fai file, then export the value.
GENOME_SIZE_BP <- suppressWarnings(as.numeric(Sys.getenv("GENOME_SIZE_BP", unset = "")))
if (length(GENOME_SIZE_BP) != 1L || !is.finite(GENOME_SIZE_BP) || GENOME_SIZE_BP <= 0) {
  stop("Set GENOME_SIZE_BP to the verified total assembly length before running script 03")
}

# Cumulative nearest-distance thresholds, in bp.
THRESHOLDS <- c(1000, 5000, 10000, 20000)

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
get_attr <- function(x, key) sub(paste0(".*", key, "=([^;]+).*"), "\\1", x)

# STEP 2: Reuse cached intervals or build masks from EDTA classifications.
# reduce() merges overlapping and directly adjacent intervals, so resulting
# Helitron loci need not correspond one-to-one with original EDTA records.
if (file.exists(MASK_HEL_RDS)) {
  message("[1/5] Loading cached Helitron mask...")
  helitron_gr <- readRDS(MASK_HEL_RDS)
} else {
  message("[1/5] No cached Helitron mask found -- building from EDTA GFF3 (slow step)...")
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
  helitron_gr <- reduce(te_gr[te_gr$classification == "DNA/Helitron"], ignore.strand = TRUE)
  saveRDS(helitron_gr, MASK_HEL_RDS)
}
message(sprintf("  %d non-redundant Helitron loci (after collapsing overlapping/nested EDTA calls).",
                length(helitron_gr)))

message("[2/5] Getting whole-gene coordinates for the general annotation...")
# STEP 3: Read nine-column GFF3 annotations; coordinates are one-based inclusive.
# GFF3 inputs must contain feature rows only after comments (no FASTA section).
gff <- fread(
  cmd = paste("grep -v '^#'", shQuote(GENOME_GFF)),
  header = FALSE, sep = "\t", quote = "",
  col.names = c("seqid", "source", "type", "start", "end", "score", "strand",
                "phase", "attributes")
)
# Select gene spans; remove duplicate records before counting genes.
genes <- gff[type == "gene"]
genes[, gene_id := get_attr(attributes, "ID")]
genes <- unique(genes, by = c("gene_id", "seqid", "start", "end"))
# If a gene ID has multiple coordinate sets, retain the first and warn.
# Review these warnings before interpreting the outputs.
still_dup <- genes$gene_id[duplicated(genes$gene_id)]
if (length(still_dup) > 0) {
  warning(sprintf("%d gene_id(s) map to >1 different coordinate set -- keeping first occurrence.",
                  length(unique(still_dup))))
  genes <- unique(genes, by = "gene_id")
}

if (!is.null(MIN_CDS_LENGTH)) {
  message(sprintf("  MIN_CDS_LENGTH sensitivity filter active: restricting to genes with CDS >= %d bp.",
                  MIN_CDS_LENGTH))
  mrna <- gff[type == "mRNA"]
  mrna[, transcript_id := get_attr(attributes, "ID")]
  mrna[, gene_id := get_attr(attributes, "Parent")]
  tx2gene <- setNames(mrna$gene_id, mrna$transcript_id)
  cds <- gff[type == "CDS"]
  cds[, transcript_id := get_attr(attributes, "Parent")]
  cds[, gene_id := tx2gene[transcript_id]]
  cds_len <- cds[, .(cds_length = sum(end - start + 1)), by = gene_id]
  genes <- merge(genes, cds_len, by = "gene_id", all.x = TRUE)
  n_before <- nrow(genes)
  genes <- genes[!is.na(cds_length) & cds_length >= MIN_CDS_LENGTH]
  message(sprintf("  Kept %d / %d genes after the CDS-length filter.", nrow(genes), n_before))
}

all_genes_gr <- GRanges(seqnames = genes$seqid,
                        ranges = IRanges(start = genes$start, end = genes$end),
                        gene_id = genes$gene_id)
message(sprintf("  %d general-annotation genes used for the 'any gene' comparison.", length(all_genes_gr)))

genic_footprint_bp <- sum(width(reduce(all_genes_gr, ignore.strand = TRUE)))
pct_genome_genic <- 100 * genic_footprint_bp / GENOME_SIZE_BP
message(sprintf("  This gene set covers %.1f Mb (%.1f%% of the ~%.2f Gb genome) -- that's the 'by chance' overlap rate to compare against.",
                genic_footprint_bp / 1e6, pct_genome_genic, GENOME_SIZE_BP / 1e9))

message("[3/5] Getting whole-gene coordinates for the canonical R-gene set...")
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
rgenes_gr <- GRanges(seqnames = rgenes_raw$seqid,
                     ranges = IRanges(start = rgenes_raw$start, end = rgenes_raw$end),
                     gene_id = rgenes_raw$gene_id)
message(sprintf("  %d R genes with coordinates.", length(rgenes_gr)))

message("[4/5] Computing, per Helitron, distance to nearest gene (0 = overlaps a gene)...")
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

dist_any_gene <- dist_to_nearest(helitron_gr, all_genes_gr)
dist_rgene <- dist_to_nearest(helitron_gr, rgenes_gr)

# STEP 5: Export one row per reduced Helitron locus and both gene distances.
per_helitron <- data.table(
  helitron_id = sprintf("%s:%d-%d", as.character(seqnames(helitron_gr)), start(helitron_gr), end(helitron_gr)),
  chrom = as.character(seqnames(helitron_gr)),
  start = start(helitron_gr),
  end = end(helitron_gr),
  width = width(helitron_gr),
  dist_to_any_gene = dist_any_gene,
  dist_to_rgene = dist_rgene
)
fwrite(per_helitron, OUT_CSV)

message("[5/5] Writing summary...")
# Report measurable denominators separately; NA values do not enter percentages.
report_thresholds <- function(dist_vec, label) {
  n_total <- length(dist_vec)
  ok <- !is.na(dist_vec)
  n_na <- sum(!ok)
  dist_vec <- dist_vec[ok]
  n <- length(dist_vec)
  out <- c(sprintf("--- Distance from each Helitron to nearest %s ---", label))
  out <- c(out, sprintf("Helitrons with a measurable distance: %d / %d (excluded: %d, no %s on that scaffold at all)",
                        n, n_total, n_na, label))
  if (n == 0) return(c(out, "  (no Helitrons to report on -- skipping)"))
  out <- c(out, sprintf("Zero distance to a gene (overlapping or adjacent): %d (%.1f%%)",
                        sum(dist_vec == 0), 100 * mean(dist_vec == 0)))
  for (t in THRESHOLDS) {
    out <- c(out, sprintf("Within %d bp of a gene (cumulative): %d (%.1f%%)",
                          t, sum(dist_vec <= t), 100 * mean(dist_vec <= t)))
  }
  out <- c(out, sprintf("Median distance to nearest %s: %.0f bp", label, median(dist_vec)))
  out
}

lines <- c()
lines <- c(lines, sprintf("Total non-redundant Helitron loci: %d", length(helitron_gr)))
lines <- c(lines, sprintf("(from %s, intact-EDTA reference)", if (file.exists(MASK_HEL_RDS)) "cached mask" else TE_GFF))
lines <- c(lines, "")
lines <- c(lines, sprintf("General-annotation gene set used: %d genes%s, covering %.1f%% of the genome (this is the 'genic coverage baseline' rate).",
                          length(all_genes_gr),
                          if (!is.null(MIN_CDS_LENGTH)) sprintf(" (CDS >= %d bp filter applied)", MIN_CDS_LENGTH) else "",
                          pct_genome_genic))
lines <- c(lines, report_thresholds(per_helitron$dist_to_any_gene, "gene (general annotation)"))
observed_pct_overlap <- 100 * mean(per_helitron$dist_to_any_gene[!is.na(per_helitron$dist_to_any_gene)] == 0)
lines <- c(lines, sprintf("Descriptive ratio: observed zero-distance proportion %.1f%% vs. %.1f%% gene-span coverage baseline = %.2fx",
                          observed_pct_overlap, pct_genome_genic, observed_pct_overlap / pct_genome_genic))
lines <- c(lines, "")
lines <- c(lines, report_thresholds(per_helitron$dist_to_rgene, "R gene (canonical set)"))

# FINAL STEP: Save a human-readable summary and print it to the console.
writeLines(lines, OUT_SUMMARY)
writeLines(lines)
cat("\nWrote per-Helitron distances to", OUT_CSV, "\n")
cat("Wrote summary to", OUT_SUMMARY, "\n")
