#!/usr/bin/env Rscript
# Workflow step 1/5: validate linkage-map metadata and write one-base SNP BED.
# Coordinates refer to the original genotyping reference, not RegenSY.

# ============================================================
# Prepare linkage-map marker coordinates from a CSV file
#
# Input:
#   CSV linkage-map file containing:
#   Marker.Name, LG, Ref.Chrom, Ref.Position,
#   Ref.Allele.x, Alt.Allele.x, Map.Position
#
# Outputs:
#   linkage_marker_metadata.tsv
#   linkage_markers.snp.bed
#   linkage_marker_summary.tsv
# ============================================================

library(dplyr)
library(readr)
library(stringr)

# ============================================================
# USER SETTINGS
# ============================================================
# Run all five scripts from the same project directory. Override PROJECT_DIR
# to use another working directory without editing these public scripts.
project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
if (!dir.exists(project_dir)) stop("PROJECT_DIR does not exist: ", project_dir)
setwd(project_dir)
# Linkage-map CSV file
map_file <- Sys.getenv("MAP_FILE", unset = "data/linkage_map.csv")

# ------------------------------------------------------------
# Chromosome-name format in the ORIGINAL reference FASTA
# ------------------------------------------------------------
#
# First check your original reference FASTA:
#
# grep '^>' original_reference_genome.fa | head
#
# Choose ONE:
#
# "marker"
#     Marker.Name begins with the exact chromosome name:
#     Chr1_745191 -> Chr1
#
# "refchrom"
#     Use Ref.Chrom exactly:
#     1 -> 1
#
# "add_Chr"
#     Add "Chr" to Ref.Chrom:
#     1 -> Chr1
#
chromosome_name_mode <- "marker"

# Optional additional prefix.
#
# Leave empty for FASTA headers such as:
# >Chr1
#
# Use an appropriate prefix if FASTA headers resemble:
# >ZhongmuNo.1_Chr1
#
# chromosome_prefix <- "ZhongmuNo.1_"
#
chromosome_prefix <- ""

# ============================================================
# CHECK INPUT
# ============================================================

if (!file.exists(map_file)) {
  stop(
    "The map file was not found: ",
    normalizePath(map_file, mustWork = FALSE)
  )
}

# ============================================================
# READ CSV FILE
# ============================================================

map_raw <- read_csv(
  file = map_file,
  na = c("", "NA", "NaN", ".", "null", "NULL"),
  show_col_types = FALSE,
  trim_ws = TRUE
)

cat("Rows read from CSV:", nrow(map_raw), "\n")
cat("Columns read from CSV:", ncol(map_raw), "\n\n")

cat("Column names found in the CSV:\n")
print(names(map_raw))

# ============================================================
# CHECK REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "Marker.Name",
  "LG",
  "Ref.Chrom",
  "Ref.Position",
  "Ref.Allele.x",
  "Alt.Allele.x",
  "Map.Position"
)

missing_columns <- setdiff(
  required_columns,
  names(map_raw)
)

if (length(missing_columns) > 0) {
  stop(
    "Required columns are missing: ",
    paste(missing_columns, collapse = ", "),
    "\n\nColumns found were:\n",
    paste(names(map_raw), collapse = ", ")
  )
}

# ============================================================
# CLEAN REQUIRED COLUMNS
# ============================================================

markers <- map_raw %>%
  transmute(
    Marker.Name = str_trim(as.character(Marker.Name)),

    LG = str_trim(as.character(LG)),

    Ref.Chrom = str_trim(as.character(Ref.Chrom)),

    Ref.Position = suppressWarnings(
      as.integer(as.numeric(Ref.Position))
    ),

    Ref.Allele = toupper(
      str_trim(as.character(Ref.Allele.x))
    ),

    Alt.Allele = toupper(
      str_trim(as.character(Alt.Allele.x))
    ),

    Map.Position = suppressWarnings(
      as.numeric(Map.Position)
    ),

    Dose.in.P1 = if ("Dose.in.P1" %in% names(map_raw)) {
      suppressWarnings(as.numeric(map_raw$Dose.in.P1))
    } else {
      NA_real_
    },

    Dose.in.P2 = if ("Dose.in.P2" %in% names(map_raw)) {
      suppressWarnings(as.numeric(map_raw$Dose.in.P2))
    } else {
      NA_real_
    }
  ) %>%
  filter(
    !is.na(Marker.Name),
    Marker.Name != "",
    !is.na(Ref.Chrom),
    Ref.Chrom != "",
    !is.na(Ref.Position),
    Ref.Position >= 1,
    !is.na(Map.Position)
  )

cat(
  "Rows retained after removing blank or invalid rows:",
  nrow(markers),
  "\n"
)

if (nrow(markers) == 0) {
  stop(
    "No valid markers remained after filtering. ",
    "Check the column contents and input delimiter."
  )
}

# ============================================================
# CREATE FASTA CHROMOSOME NAMES
# ============================================================

if (chromosome_name_mode == "marker") {

  markers <- markers %>%
    mutate(
      Fasta.Chrom = str_extract(
        Marker.Name,
        "^[^_]+"
      )
    )

} else if (chromosome_name_mode == "refchrom") {

  markers <- markers %>%
    mutate(
      Fasta.Chrom = Ref.Chrom
    )

} else if (chromosome_name_mode == "add_Chr") {

  markers <- markers %>%
    mutate(
      Fasta.Chrom = if_else(
        str_detect(
          Ref.Chrom,
          regex("^Chr", ignore_case = TRUE)
        ),
        Ref.Chrom,
        paste0("Chr", Ref.Chrom)
      )
    )

} else {

  stop(
    "Invalid chromosome_name_mode: ",
    chromosome_name_mode,
    "\nUse 'marker', 'refchrom' or 'add_Chr'."
  )
}

markers <- markers %>%
  mutate(
    Fasta.Chrom = paste0(
      chromosome_prefix,
      Fasta.Chrom
    )
  )

# ============================================================
# CHECK CHROMOSOME NAMES
# ============================================================

missing_chromosome_names <- markers %>%
  filter(
    is.na(Fasta.Chrom) |
      Fasta.Chrom == ""
  )

if (nrow(missing_chromosome_names) > 0) {
  stop(
    nrow(missing_chromosome_names),
    " markers do not have a valid FASTA chromosome name."
  )
}

# ============================================================
# CHECK FOR DUPLICATED MARKER NAMES
# ============================================================

duplicated_markers <- markers %>%
  count(
    Marker.Name,
    name = "Number.of.Rows"
  ) %>%
  filter(Number.of.Rows > 1)

if (nrow(duplicated_markers) > 0) {

  warning(
    nrow(duplicated_markers),
    " duplicated Marker.Name values were found. ",
    "Only the first occurrence will be retained."
  )

  write_tsv(
    duplicated_markers,
    "duplicated_marker_names.tsv"
  )
}

markers <- markers %>%
  distinct(
    Marker.Name,
    .keep_all = TRUE
  )

# ============================================================
# CREATE BED COORDINATES
# ============================================================
#
# Ref.Position is assumed to use one-based coordinates.
#
# BED uses zero-based, half-open coordinates.
#
# Example:
#     Ref.Position = 745191
#
# BED:
#     start = 745190
#     end   = 745191
#
# The Bash script will expand this one-base interval by
# 250 bp on both sides.
# ============================================================

markers <- markers %>%
  mutate(
    SNP.BED.Start = Ref.Position - 1L,
    SNP.BED.End = Ref.Position
  )

# ============================================================
# CHECK ALLELE VALUES
# ============================================================

valid_bases <- c("A", "C", "G", "T")

markers <- markers %>%
  mutate(
    Ref.Allele.Valid =
      !is.na(Ref.Allele) &
      Ref.Allele %in% valid_bases,

    Alt.Allele.Valid =
      !is.na(Alt.Allele) &
      Alt.Allele %in% valid_bases
  )

invalid_alleles <- markers %>%
  filter(
    !Ref.Allele.Valid |
      !Alt.Allele.Valid
  )

if (nrow(invalid_alleles) > 0) {

  warning(
    nrow(invalid_alleles),
    " markers contain missing or non-ACGT alleles. ",
    "They will remain in the BED file."
  )

  write_tsv(
    invalid_alleles,
    "markers_with_nonstandard_alleles.tsv",
    na = "NA"
  )
}

# ============================================================
# SORT MARKERS
# ============================================================

markers <- markers %>%
  mutate(
    LG.Numeric = suppressWarnings(
      as.numeric(LG)
    )
  ) %>%
  arrange(
    LG.Numeric,
    LG,
    Map.Position,
    Ref.Position,
    Marker.Name
  )

# ============================================================
# WRITE MARKER METADATA
# ============================================================

marker_metadata <- markers %>%
  select(
    Marker.Name,
    LG,
    Map.Position,
    Ref.Chrom,
    Fasta.Chrom,
    Ref.Position,
    Ref.Allele,
    Alt.Allele,
    Dose.in.P1,
    Dose.in.P2,
    SNP.BED.Start,
    SNP.BED.End,
    Ref.Allele.Valid,
    Alt.Allele.Valid
  )

write_tsv(
  marker_metadata,
  "linkage_marker_metadata.tsv",
  na = "NA"
)

# ============================================================
# WRITE ONE-BASE SNP BED FILE
# ============================================================

snp_bed <- markers %>%
  transmute(
    chromosome = Fasta.Chrom,
    start = SNP.BED.Start,
    end = SNP.BED.End,
    marker = Marker.Name
  )

write_tsv(
  snp_bed,
  "linkage_markers.snp.bed",
  col_names = FALSE,
  na = "NA"
)

# ============================================================
# WRITE SUMMARY
# ============================================================

summary_table <- tibble(
  Measurement = c(
    "Rows in original CSV file",
    "Usable unique markers",
    "Number of linkage groups",
    "Number of reference chromosomes",
    "Markers with valid reference allele",
    "Markers with valid alternate allele",
    "Minimum map position (cM)",
    "Maximum map position (cM)"
  ),

  Value = c(
    nrow(map_raw),
    nrow(markers),
    n_distinct(markers$LG),
    n_distinct(markers$Fasta.Chrom),
    sum(markers$Ref.Allele.Valid, na.rm = TRUE),
    sum(markers$Alt.Allele.Valid, na.rm = TRUE),
    min(markers$Map.Position, na.rm = TRUE),
    max(markers$Map.Position, na.rm = TRUE)
  )
)

write_tsv(
  summary_table,
  "linkage_marker_summary.tsv"
)

# ============================================================
# FINAL OUTPUT
# ============================================================

cat("\n============================================\n")
cat("LINKAGE-MAP PREPARATION COMPLETED\n")
cat("============================================\n\n")

print(summary_table)

cat("\nFASTA chromosome names placed in the BED file:\n")
print(sort(unique(markers$Fasta.Chrom)))

cat("\nFirst six BED records:\n")
print(head(snp_bed))

cat("\nFiles created:\n")
cat("  linkage_marker_metadata.tsv\n")
cat("  linkage_markers.snp.bed\n")
cat("  linkage_marker_summary.tsv\n")

if (file.exists("duplicated_marker_names.tsv")) {
  cat("  duplicated_marker_names.tsv\n")
}

if (file.exists("markers_with_nonstandard_alleles.tsv")) {
  cat("  markers_with_nonstandard_alleles.tsv\n")
}

cat("\nBefore running the Bash script, confirm that the BED\n")
