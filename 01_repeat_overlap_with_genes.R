#!/usr/bin/env Rscript
# 01_repeat_overlap_with_genes.R
# Measure repeat overlap across whole-gene spans and reduced CDS footprints.
# General-annotation genes overlapping canonical R-gene loci are excluded
# from the background used for optional R-gene distribution comparisons.
# CDS features are preferred; an exon fallback includes UTRs and is not CDS.
# R-gene and general annotations have separate ID namespaces: compare loci
# by coordinates, not by matching gene IDs across these annotations.
# Requirements: data.table, GenomicRanges, IRanges, and Unix grep.
# Run from the project directory; relative paths resolve there.
# Usage: Rscript 01_repeat_overlap_with_genes.R

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
OUT_CSV       <- file.path(OUTDIR, "per_gene_genome_repeat_overlap.csv")
OUT_SUMMARY   <- file.path(OUTDIR, "genome_repeat_overlap_summary.txt")

# Optional externally generated comparison table; not produced by this set.
RGENE_CDS_OVERLAP      <- "data/rgene_comparisons/rgene_final_confidence_per_gene.csv"
RGENE_WHOLEGENE_OVERLAP <- "data/rgene_comparisons/rgene_confidence_per_gene.csv"      

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# STEP 2: Reuse cached intervals or build masks from EDTA classifications.
# reduce() merges overlapping and directly adjacent intervals, so resulting
# Helitron loci need not correspond one-to-one with original EDTA records.
if (file.exists(MASK_ANY_RDS) && file.exists(MASK_HEL_RDS)) {
  message("[1/6] Loading cached repeat masks...")
  any_repeat_mask <- readRDS(MASK_ANY_RDS)
  helitron_mask <- readRDS(MASK_HEL_RDS)
} else {
  message("[1/6] No cached masks found -- rebuilding from EDTA GFF3 (slow step)...")
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

message("[2/6] Parsing genome-wide GFF3...")
# STEP 3: Read nine-column GFF3 annotations; coordinates are one-based inclusive.
# GFF3 inputs must contain feature rows only after comments (no FASTA section).
gff <- fread(
  cmd = paste("grep -v '^#'", shQuote(GENOME_GFF)),
  header = FALSE, sep = "\t", quote = "",
  col.names = c("seqid", "source", "type", "start", "end", "score", "strand",
                "phase", "attributes")
)
message("  Feature type counts in this file:")
print(table(gff$type))

get_attr <- function(x, key) sub(paste0(".*", key, "=([^;]+).*"), "\\1", x)

has_cds <- "CDS" %in% gff$type
if (has_cds) {
  message("  'CDS' features found -- using CDS-only footprint (matches the R-gene methodology).")
} else {
  warning("No 'CDS' rows found in this GFF3 -- falling back to 'exon' footprint for the coding-only metric. This is NOT equivalent to true CDS (it will include UTRs), so treat pct_*_cdslike as approximate, not directly comparable to the R-gene CDS-only numbers.")
}

message("[3/6] Extracting whole-gene coordinates and de-duplicating...")
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
    "%d gene_id(s) map to >1 DIFFERENT coordinate set (e.g. %s) -- keeping only the first occurrence. This is a real non-unique-ID situation, not a simple file duplicate -- investigate before trusting those genes' numbers.",
    length(unique(still_dup)), paste(head(unique(still_dup), 5), collapse = ", ")
  ))
  genes <- unique(genes, by = "gene_id")
}
message(sprintf("  %d unique genes with coordinates.", nrow(genes)))

genes_gr <- GRanges(seqnames = genes$seqid,
                    ranges = IRanges(start = genes$start, end = genes$end),
                    gene_id = genes$gene_id)

message("[3b/6] Flagging general-annotation genes that overlap an R-gene locus...")
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
rgenes_raw <- unique(rgenes_raw, by = "gene_id")  # belt-and-suspenders, see 05 for the full check
rgene_loci_gr <- GRanges(seqnames = rgenes_raw$seqid,
                         ranges = IRanges(start = rgenes_raw$start, end = rgenes_raw$end))

# Flag general genes overlapping an R-gene locus by genomic coordinates.
ov_rgene <- findOverlaps(genes_gr, rgene_loci_gr, ignore.strand = TRUE)
overlaps_rgene_locus <- logical(length(genes_gr))
overlaps_rgene_locus[unique(queryHits(ov_rgene))] <- TRUE
message(sprintf(
  "  %d general-annotation genes overlap at least one of the %d R-gene loci (%d R-gene loci matched no general-annotation gene at all).",
  sum(overlaps_rgene_locus), length(rgene_loci_gr),
  length(rgene_loci_gr) - length(unique(subjectHits(ov_rgene)))
))

message("[4/6] Building per-gene coding (or exon-fallback) footprint...")
# Map CDS/exon Parent transcript IDs to mRNA Parent gene IDs.
# Assumes single-parent features and mRNA records; inspect other GFF layouts.
mrna <- gff[type == "mRNA"]
mrna[, transcript_id := get_attr(attributes, "ID")]
mrna[, gene_id := get_attr(attributes, "Parent")]
tx2gene <- setNames(mrna$gene_id, mrna$transcript_id)

cds_type <- if (has_cds) "CDS" else "exon"
cds <- gff[type == cds_type]
cds[, transcript_id := get_attr(attributes, "Parent")]
cds[, gene_id := tx2gene[transcript_id]]
cds <- cds[!is.na(gene_id)]

cds_gr <- GRanges(seqnames = cds$seqid,
                  ranges = IRanges(start = cds$start, end = cds$end),
                  gene_id = cds$gene_id)

# Reduce coding blocks per gene so shared bases across isoforms count once.
cds_by_gene <- split(cds_gr, cds_gr$gene_id)
cds_reduced <- reduce(cds_by_gene, ignore.strand = TRUE)
flat <- unlist(cds_reduced, use.names = FALSE)
mcols(flat)$gene_id <- rep(names(cds_reduced), times = elementNROWS(cds_reduced))

# STEP 4: Sum actual intersecting bases against the non-redundant repeat mask.
overlap_bp_per_interval <- function(query_gr, mask_gr) {
  hits <- findOverlaps(query_gr, mask_gr, ignore.strand = TRUE)
  result <- integer(length(query_gr))
  if (length(hits) > 0) {
    inter <- pintersect(query_gr[queryHits(hits)], mask_gr[subjectHits(hits)],
                        ignore.strand = TRUE)
    bp_by_query <- tapply(width(inter), queryHits(hits), sum)
    result[as.integer(names(bp_by_query))] <- bp_by_query
  }
  result
}

message("[5/6] Computing whole-gene-span and CDS-only overlap (this is the slow step at genome scale)...")
bp_any_wg <- overlap_bp_per_interval(genes_gr, any_repeat_mask)
bp_hel_wg <- overlap_bp_per_interval(genes_gr, helitron_mask)

interval_bp_any <- overlap_bp_per_interval(flat, any_repeat_mask)
interval_bp_hel <- overlap_bp_per_interval(flat, helitron_mask)
interval_dt <- data.table(
  gene_id = mcols(flat)$gene_id,
  block_length = width(flat),
  bp_any = interval_bp_any,
  bp_hel = interval_bp_hel
)
cds_per_gene <- interval_dt[, .(
  cds_length = sum(block_length),
  bp_overlap_any_repeat_cds = sum(bp_any),
  bp_overlap_helitron_cds = sum(bp_hel)
), by = gene_id]

# STEP 5: Assemble per-gene metrics, preserving unmeasured values as NA.
per_gene <- data.table(
  gene_id = genes_gr$gene_id,
  chrom = as.character(seqnames(genes_gr)),
  start = start(genes_gr),
  end = end(genes_gr),
  gene_length = width(genes_gr),
  overlaps_rgene_locus = overlaps_rgene_locus,
  bp_overlap_any_repeat_wholegene = bp_any_wg,
  bp_overlap_helitron_wholegene = bp_hel_wg
)
per_gene[, pct_any_repeat_wholegene := 100 * bp_overlap_any_repeat_wholegene / gene_length]
per_gene[, pct_helitron_wholegene := 100 * bp_overlap_helitron_wholegene / gene_length]

per_gene <- merge(per_gene, cds_per_gene, by = "gene_id", all.x = TRUE)
per_gene[, pct_any_repeat_cds := 100 * bp_overlap_any_repeat_cds / cds_length]
per_gene[, pct_helitron_cds := 100 * bp_overlap_helitron_cds / cds_length]

fwrite(per_gene, OUT_CSV)

message("[6/6] Writing summary...")
lines <- c()
lines <- c(lines, sprintf("Total genes in general annotation (%s): %d", GENOME_GFF, nrow(per_gene)))
lines <- c(lines, sprintf("Coding-footprint feature type used: %s%s", cds_type,
                          if (!has_cds) " (fallback -- no CDS rows in this file)" else ""))
lines <- c(lines, "")

n_no_cds <- sum(is.na(per_gene$pct_any_repeat_cds))
if (n_no_cds > 0) {
  lines <- c(lines, sprintf("NOTE: %d genes had no %s rows and got NA for the CDS-only metric (left NA, not 0).",
                            n_no_cds, cds_type))
  lines <- c(lines, "")
}

summarize_group <- function(dt, wg_col, cds_col, label, n_label = NULL) {
  out <- c(sprintf("--- %s (n=%d) ---", label, if (is.null(n_label)) nrow(dt) else n_label))
  wg_vals <- dt[[wg_col]][!is.na(dt[[wg_col]])]
  out <- c(out, sprintf("Whole-gene %%overlap any-repeat: median=%.2f%%, mean=%.2f%%, %%genes>0%%=%.1f%%, n=%d",
                        median(wg_vals), mean(wg_vals), 100 * mean(wg_vals > 0), length(wg_vals)))
  cds_vals <- dt[[cds_col]][!is.na(dt[[cds_col]])]
  out <- c(out, sprintf("CDS-only %%overlap any-repeat: median=%.2f%%, mean=%.2f%%, %%genes>0%%=%.1f%%, n_measured=%d",
                        median(cds_vals), mean(cds_vals), 100 * mean(cds_vals > 0), length(cds_vals)))
  out
}

lines <- c(lines, summarize_group(per_gene, "pct_any_repeat_wholegene", "pct_any_repeat_cds",
                                  "All genes, general annotation (general annotation, INCLUDES R-gene loci)"))
lines <- c(lines, "")

# Exclude R-gene-overlapping general models for the comparison background.
non_rgene <- per_gene[overlaps_rgene_locus == FALSE]
lines <- c(lines, summarize_group(non_rgene, "pct_any_repeat_wholegene", "pct_any_repeat_cds",
                                  "Non-R genes only, general annotation (R-gene-overlapping loci excluded -- this is the correct background)"))
lines <- c(lines, "")

# Optional comparison: summarize external R-gene metrics and run Wilcoxon.
# These are unadjusted tests; spatial dependence is not modeled.
if (file.exists(RGENE_CDS_OVERLAP)) {
  rgene_cds_dt <- fread(RGENE_CDS_OVERLAP)
  cds_vals_r <- rgene_cds_dt$pct_any_repeat_cds[!is.na(rgene_cds_dt$pct_any_repeat_cds)]
  lines <- c(lines, sprintf("--- R genes (BRAKER/FindPlantNLRs set, from %s, n=%d, n_measured=%d) ---",
                            RGENE_CDS_OVERLAP, nrow(rgene_cds_dt), length(cds_vals_r)))
  lines <- c(lines, sprintf("CDS-only %%overlap any-repeat: median=%.2f%%, mean=%.2f%%, %%genes>0%%=%.1f%%",
                            median(cds_vals_r), mean(cds_vals_r), 100 * mean(cds_vals_r > 0)))
  lines <- c(lines, "")

  nonrgene_cds_vals <- non_rgene$pct_any_repeat_cds[!is.na(non_rgene$pct_any_repeat_cds)]
  wt_cds <- wilcox.test(cds_vals_r, nonrgene_cds_vals)
  lines <- c(lines, sprintf("Wilcoxon (CDS-only %%overlap, R genes vs. non-R genes ONLY): W=%.3e, p=%.3e",
                            wt_cds$statistic, wt_cds$p.value))
  lines <- c(lines, sprintf("  medians: R genes=%.2f%%, non-R genes=%.2f%% (check which is larger before interpreting the direction of the difference)",
                            median(cds_vals_r), median(nonrgene_cds_vals)))
} else {
  lines <- c(lines, sprintf(
    "NOTE: %s not found in the working directory -- skipping the R-gene-vs-genome CDS-overlap comparison. Provide the optional R-gene CDS-overlap table; see README.md.",
    RGENE_CDS_OVERLAP))
}
lines <- c(lines, "")

# Optional comparison: summarize external R-gene metrics and run Wilcoxon.
# These are unadjusted tests; spatial dependence is not modeled.
if (file.exists(RGENE_WHOLEGENE_OVERLAP)) {
  rgene_wg_dt <- fread(RGENE_WHOLEGENE_OVERLAP)
  lines <- c(lines, sprintf(
    "--- R genes, whole-gene-span (from %s, n=%d -- externally generated comparison set; verify its membership) ---",
    RGENE_WHOLEGENE_OVERLAP, nrow(rgene_wg_dt)))
  nonrgene_wg_vals <- non_rgene$pct_any_repeat_wholegene[!is.na(non_rgene$pct_any_repeat_wholegene)]
  wt_wg <- wilcox.test(rgene_wg_dt$pct_any_repeat, nonrgene_wg_vals)
  lines <- c(lines, sprintf("Wilcoxon (whole-gene %%overlap, R genes vs. non-R genes ONLY): W=%.3e, p=%.3e",
                            wt_wg$statistic, wt_wg$p.value))
  lines <- c(lines, sprintf("  medians: R genes=%.2f%%, non-R genes=%.2f%%",
                            median(rgene_wg_dt$pct_any_repeat), median(nonrgene_wg_vals)))
} else {
  lines <- c(lines, sprintf(
    "NOTE: %s not found in the working directory -- skipping the whole-gene-span comparison.",
    RGENE_WHOLEGENE_OVERLAP))
}

# FINAL STEP: Save a human-readable summary and print it to the console.
writeLines(lines, OUT_SUMMARY)
writeLines(lines)
cat("\nWrote per-gene table to", OUT_CSV, "\n")
cat("Wrote summary to", OUT_SUMMARY, "\n")
