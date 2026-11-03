#!/usr/bin/env Rscript
# Workflow step 4/5: quality-filter alignments and project SNP coordinates.
# Requires outputs from steps 01, 02, and 03.

# ============================================================
# Filter marker alignments and project the original SNP position
# through each CIGAR string onto the RegenSY assembly.
#
# Inputs:
#   RegenSY_marker_alignments/
#       all_haplotypes.primary_positions.tsv
#
#   linkage_marker_metadata.tsv
#
# Outputs:
#   RegenSY_marker_alignments/
#       all_alignment_QC.tsv
#       rejected_alignments.tsv
#       valid_alignment_candidates.tsv
#       best_marker_positions.tsv
#       markers_without_valid_position.tsv
#       genetic_physical_map.tsv
#       genetic_physical_map_wide.tsv
#       marker_mapping_summary.tsv
# ============================================================

library(dplyr)
library(readr)
library(stringr)
library(tidyr)
library(ggplot2)

# ============================================================
# USER SETTINGS
# ============================================================
# Run all five scripts from the same project directory. Override PROJECT_DIR
# to use another working directory without editing these public scripts.
project_dir <- Sys.getenv("PROJECT_DIR", unset = getwd())
if (!dir.exists(project_dir)) stop("PROJECT_DIR does not exist: ", project_dir)
setwd(project_dir)
alignment_file <- file.path(
  "RegenSY_marker_alignments",
  "all_haplotypes.primary_positions.tsv"
)

metadata_file <- "linkage_marker_metadata.tsv"

output_directory <- "RegenSY_marker_alignments"

# Minimum mapping quality
minimum_mapq <- 20

# At least 80% of the 501-bp marker sequence must participate
# in the alignment. This removes short clipped alignments.
minimum_query_coverage <- 0.80

# Optional: require marker chromosome and RegenSY chromosome to agree.
# FALSE retains discordant chromosomes for comparison in the plots.
# For example, Chr1_745191 must align to Chr1.
require_matching_chromosome <- FALSE

# If RegenSY names have extra haplotype prefixes, chromosome
# normalization below may need to be adjusted.
# ============================================================

dir.create(
  output_directory,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# CHECK INPUTS
# ============================================================

if (!file.exists(alignment_file)) {
  stop("Alignment file not found: ", alignment_file)
}

if (!file.exists(metadata_file)) {
  stop("Metadata file not found: ", metadata_file)
}

# ============================================================
# READ FILES
# ============================================================

alignments <- read_tsv(
  alignment_file,
  show_col_types = FALSE
)

metadata <- read_tsv(
  metadata_file,
  show_col_types = FALSE
)

required_alignment_columns <- c(
  "Marker.Name",
  "Haplotype",
  "RegenSY.Sequence",
  "Alignment.Start",
  "MAPQ",
  "Strand",
  "CIGAR",
  "Query.Length"
)

missing_alignment_columns <- setdiff(
  required_alignment_columns,
  names(alignments)
)

if (length(missing_alignment_columns) > 0) {
  stop(
    "Missing alignment columns: ",
    paste(missing_alignment_columns, collapse = ", ")
  )
}

required_metadata_columns <- c(
  "Marker.Name",
  "LG",
  "Map.Position",
  "Ref.Position",
  "Ref.Allele",
  "Alt.Allele"
)

missing_metadata_columns <- setdiff(
  required_metadata_columns,
  names(metadata)
)

if (length(missing_metadata_columns) > 0) {
  stop(
    "Missing metadata columns: ",
    paste(missing_metadata_columns, collapse = ", ")
  )
}

# ============================================================
# CIGAR PARSING FUNCTIONS
# ============================================================

parse_cigar <- function(cigar) {

  if (
    is.na(cigar) ||
    cigar == "" ||
    cigar == "*"
  ) {
    return(
      tibble(
        Length = integer(),
        Operation = character()
      )
    )
  }

  lengths <- as.integer(
    unlist(str_extract_all(cigar, "[0-9]+"))
  )

  operations <- unlist(
    str_extract_all(cigar, "[MIDNSHP=X]")
  )

  if (length(lengths) != length(operations)) {
    stop("Could not parse CIGAR string: ", cigar)
  }

  tibble(
    Length = lengths,
    Operation = operations
  )
}

# Query bases aligned to the reference:
# M, I, = and X consume query and are part of the alignment.
# Soft-clipped bases are excluded.
calculate_aligned_query_bases <- function(cigar) {

  cigar_table <- parse_cigar(cigar)

  if (nrow(cigar_table) == 0) {
    return(NA_integer_)
  }

  sum(
    cigar_table$Length[
      cigar_table$Operation %in% c("M", "I", "=", "X")
    ]
  )
}

calculate_soft_clipped_bases <- function(cigar) {

  cigar_table <- parse_cigar(cigar)

  if (nrow(cigar_table) == 0) {
    return(NA_integer_)
  }

  sum(
    cigar_table$Length[
      cigar_table$Operation == "S"
    ]
  )
}

calculate_reference_span <- function(cigar) {

  cigar_table <- parse_cigar(cigar)

  if (nrow(cigar_table) == 0) {
    return(NA_integer_)
  }

  sum(
    cigar_table$Length[
      cigar_table$Operation %in% c(
        "M", "D", "N", "=", "X"
      )
    ]
  )
}

# Project one query base through the CIGAR string.
#
# Query position is relative to the sequence orientation represented
# by the CIGAR. For a reverse-strand alignment, convert the original
# query coordinate before projection.
project_query_base <- function(
    alignment_start,
    cigar,
    query_length,
    original_query_position,
    strand) {

  if (
    is.na(alignment_start) ||
    is.na(cigar) ||
    is.na(query_length) ||
    is.na(original_query_position) ||
    is.na(strand)
  ) {
    return(NA_integer_)
  }

  if (strand == "-") {
    target_query_position <-
      query_length - original_query_position + 1L
  } else {
    target_query_position <- original_query_position
  }

  cigar_table <- parse_cigar(cigar)

  if (nrow(cigar_table) == 0) {
    return(NA_integer_)
  }

  query_position <- 1L
  reference_position <- as.integer(alignment_start)

  for (i in seq_len(nrow(cigar_table))) {

    operation <- cigar_table$Operation[i]
    operation_length <- cigar_table$Length[i]

    if (operation %in% c("M", "=", "X")) {

      query_end <- query_position + operation_length - 1L

      if (
        target_query_position >= query_position &&
        target_query_position <= query_end
      ) {
        offset <- target_query_position - query_position

        return(
          as.integer(reference_position + offset)
        )
      }

      query_position <- query_position + operation_length
      reference_position <-
        reference_position + operation_length

    } else if (operation %in% c("I", "S")) {

      query_end <- query_position + operation_length - 1L

      # The target base is inserted relative to the reference
      # or soft clipped, so it has no reference coordinate.
      if (
        target_query_position >= query_position &&
        target_query_position <= query_end
      ) {
        return(NA_integer_)
      }

      query_position <- query_position + operation_length

    } else if (operation %in% c("D", "N")) {

      reference_position <-
        reference_position + operation_length

    } else if (operation %in% c("H", "P")) {

      # Hard clipping and padding consume neither sequence here.
      next
    }
  }

  NA_integer_
}

# ============================================================
# CHROMOSOME NORMALIZATION
# ============================================================

normalize_chromosome <- function(x) {

  x <- as.character(x)

  # Extract Chr followed by a chromosome number.
  extracted <- str_extract(
    x,
    regex("Chr0*[0-9]+", ignore_case = TRUE)
  )

  extracted <- str_replace(
    extracted,
    regex("^chr", ignore_case = TRUE),
    "Chr"
  )

  extracted <- str_replace(
    extracted,
    "^Chr0+",
    "Chr"
  )

  extracted
}

# ============================================================
# CALCULATE ALIGNMENT QC
# ============================================================

# Read the actual extracted interval to handle chromosome-end clipping.
# BED starts are zero-based; Ref.Position is one-based, so the SNP's
# one-based query position is Ref.Position - Flank.Start.
flank_bed_file <- "linkage_markers.flank250.bed"
if (!file.exists(flank_bed_file)) stop("Missing extracted flank BED: ", flank_bed_file)
flanks <- read_tsv(flank_bed_file, col_names = FALSE, show_col_types = FALSE)
if (ncol(flanks) < 4L) stop("Flank BED must contain at least four columns")
flank_offsets <- tibble(
  Marker.Name = as.character(flanks[[4]]),
  Flank.Start = as.integer(flanks[[2]])
) %>% left_join(metadata %>% select(Marker.Name, Ref.Position), by = "Marker.Name") %>%
  mutate(SNP.Query.Position = as.integer(Ref.Position) - Flank.Start) %>%
  select(Marker.Name, SNP.Query.Position)
if (anyDuplicated(flank_offsets$Marker.Name)) stop("Duplicate marker IDs in flank offsets")
alignments <- alignments %>% left_join(flank_offsets, by = "Marker.Name")
if (any(is.na(alignments$SNP.Query.Position) |
        alignments$SNP.Query.Position < 1L |
        alignments$SNP.Query.Position > as.integer(alignments$Query.Length))) {
  stop("Missing or out-of-range SNP query position; check metadata and flank BED")
}

alignment_qc <- alignments %>%
  mutate(
    Alignment.Start = as.integer(Alignment.Start),
    MAPQ = as.numeric(MAPQ),
    Query.Length = as.integer(Query.Length),

    Expected.Chromosome =
      normalize_chromosome(Marker.Name),

    Observed.Chromosome =
      normalize_chromosome(RegenSY.Sequence),

    Chromosome.Match =
      !is.na(Expected.Chromosome) &
      !is.na(Observed.Chromosome) &
      Expected.Chromosome == Observed.Chromosome,

    Aligned.Query.Bases = vapply(
      CIGAR,
      calculate_aligned_query_bases,
      integer(1)
    ),

    Soft.Clipped.Bases = vapply(
      CIGAR,
      calculate_soft_clipped_bases,
      integer(1)
    ),

    Reference.Span = vapply(
      CIGAR,
      calculate_reference_span,
      integer(1)
    ),

    Query.Coverage =
      Aligned.Query.Bases / Query.Length,

    # Use the original SNP offset, including shortened boundary flanks.
    Original.SNP.Query.Position =
      as.integer(SNP.Query.Position),

    RegenSY.Position = mapply(
      FUN = project_query_base,
      alignment_start = Alignment.Start,
      cigar = CIGAR,
      query_length = Query.Length,
      original_query_position =
        Original.SNP.Query.Position,
      strand = Strand,
      USE.NAMES = FALSE
    ),

    SNP.Covered =
      !is.na(RegenSY.Position),

    Pass.MAPQ =
      !is.na(MAPQ) &
      MAPQ >= minimum_mapq,

    Pass.Query.Coverage =
      !is.na(Query.Coverage) &
      Query.Coverage >= minimum_query_coverage,

    Pass.Chromosome =
      if (require_matching_chromosome) {
        Chromosome.Match
      } else {
        TRUE
      },

    Pass.All.Filters =
      Pass.MAPQ &
      Pass.Query.Coverage &
      Pass.Chromosome &
      SNP.Covered,

    Rejection.Reason = case_when(
      !Pass.MAPQ ~ "MAPQ below threshold",
      !Pass.Query.Coverage ~ "Insufficient query coverage",
      !Pass.Chromosome ~ "Chromosome mismatch",
      !SNP.Covered ~ "SNP not aligned to reference",
      TRUE ~ "Accepted"
    )
  )

write_tsv(
  alignment_qc,
  file.path(
    output_directory,
    "all_alignment_QC.tsv"
  ),
  na = "NA"
)

# ============================================================
# SAVE REJECTED ALIGNMENTS
# ============================================================

rejected_alignments <- alignment_qc %>%
  filter(!Pass.All.Filters) %>%
  arrange(
    Marker.Name,
    Haplotype,
    Rejection.Reason
  )

write_tsv(
  rejected_alignments,
  file.path(
    output_directory,
    "rejected_alignments.tsv"
  ),
  na = "NA"
)

# ============================================================
# VALID CANDIDATE ALIGNMENTS
# ============================================================

valid_candidates <- alignment_qc %>%
  filter(Pass.All.Filters) %>%
  arrange(
    Marker.Name,
    Haplotype,
    desc(Query.Coverage),
    desc(MAPQ),
    Soft.Clipped.Bases,
    RegenSY.Sequence,
    RegenSY.Position
  ) %>%
  group_by(
    Marker.Name,
    Haplotype
  ) %>%
  mutate(
    Number.Valid.Candidates =
      n(),

    Candidate.Rank =
      row_number()
  ) %>%
  ungroup()

write_tsv(
  valid_candidates,
  file.path(
    output_directory,
    "valid_alignment_candidates.tsv"
  ),
  na = "NA"
)

# ============================================================
# SELECT BEST ALIGNMENT PER MARKER × HAPLOTYPE
# ============================================================

best_positions <- valid_candidates %>%
  filter(Candidate.Rank == 1) %>%
  select(
    Marker.Name,
    Haplotype,
    RegenSY.Sequence,
    RegenSY.Position,
    Alignment.Start,
    Strand,
    CIGAR,
    MAPQ,
    Query.Length,
    Aligned.Query.Bases,
    Query.Coverage,
    Soft.Clipped.Bases,
    Reference.Span,
    Expected.Chromosome,
    Observed.Chromosome,
    Number.Valid.Candidates
  )

write_tsv(
  best_positions,
  file.path(
    output_directory,
    "best_marker_positions.tsv"
  ),
  na = "NA"
)

# ============================================================
# JOIN GENETIC AND PHYSICAL POSITIONS
# ============================================================

genetic_physical_map <- metadata %>%
  select(
    Marker.Name,
    LG,
    Map.Position,
    Ref.Chrom,
    Ref.Position,
    Ref.Allele,
    Alt.Allele
  ) %>%
  left_join(
    best_positions,
    by = "Marker.Name"
  ) %>%
  mutate(
    RegenSY.Position.Mb =
      RegenSY.Position / 1e6
  ) %>%
  arrange(
    suppressWarnings(as.numeric(LG)),
    Map.Position,
    Marker.Name,
    Haplotype
  )

write_tsv(
  genetic_physical_map,
  file.path(
    output_directory,
    "genetic_physical_map.tsv"
  ),
  na = "NA"
)

# ============================================================
# IDENTIFY MARKERS WITHOUT A VALID POSITION
# ============================================================

markers_without_valid_position <- metadata %>%
  distinct(Marker.Name, LG, Map.Position) %>%
  anti_join(
    best_positions %>%
      distinct(Marker.Name),
    by = "Marker.Name"
  )

write_tsv(
  markers_without_valid_position,
  file.path(
    output_directory,
    "markers_without_valid_position.tsv"
  ),
  na = "NA"
)

# ============================================================
# CREATE WIDE TABLE: ONE ROW PER MARKER
# ============================================================

wide_positions <- best_positions %>%
  select(
    Marker.Name,
    Haplotype,
    RegenSY.Sequence,
    RegenSY.Position,
    MAPQ,
    Query.Coverage,
    Number.Valid.Candidates
  ) %>%
  pivot_wider(
    names_from = Haplotype,
    values_from = c(
      RegenSY.Sequence,
      RegenSY.Position,
      MAPQ,
      Query.Coverage,
      Number.Valid.Candidates
    ),
    names_glue = "{Haplotype}_{.value}"
  )

genetic_physical_map_wide <- metadata %>%
  select(
    Marker.Name,
    LG,
    Map.Position,
    Ref.Chrom,
    Ref.Position,
    Ref.Allele,
    Alt.Allele
  ) %>%
  left_join(
    wide_positions,
    by = "Marker.Name"
  ) %>%
  arrange(
    suppressWarnings(as.numeric(LG)),
    Map.Position,
    Ref.Position
  )

write_tsv(
  genetic_physical_map_wide,
  file.path(
    output_directory,
    "genetic_physical_map_wide.tsv"
  ),
  na = "NA"
)

# ============================================================
# CREATE SUMMARY
# ============================================================

total_markers <- n_distinct(metadata$Marker.Name)

mapping_summary <- alignment_qc %>%
  group_by(Haplotype) %>%
  summarise(
    Total.Input.Alignments = n(),
    Unique.Markers.With.Alignment =
      n_distinct(Marker.Name),
    Alignments.Passing.Filters =
      sum(Pass.All.Filters, na.rm = TRUE),
    Markers.With.Valid.Position =
      n_distinct(
        Marker.Name[Pass.All.Filters]
      ),
    Rejected.Low.MAPQ =
      sum(!Pass.MAPQ, na.rm = TRUE),
    Rejected.Low.Coverage =
      sum(
        Pass.MAPQ &
          !Pass.Query.Coverage,
        na.rm = TRUE
      ),
    Rejected.Chromosome.Mismatch =
      sum(
        Pass.MAPQ &
          Pass.Query.Coverage &
          !Pass.Chromosome,
        na.rm = TRUE
      ),
    Rejected.SNP.Not.Covered =
      sum(
        Pass.MAPQ &
          Pass.Query.Coverage &
          Pass.Chromosome &
          !SNP.Covered,
        na.rm = TRUE
      ),
    .groups = "drop"
  ) %>%
  mutate(
    Total.Markers.In.Map = total_markers,
    Percent.Markers.With.Valid.Position =
      100 * Markers.With.Valid.Position /
      Total.Markers.In.Map
  )

write_tsv(
  mapping_summary,
  file.path(
    output_directory,
    "marker_mapping_summary.tsv"
  ),
  na = "NA"
)

# ============================================================
# PRELIMINARY GENETIC–PHYSICAL PLOT
# ============================================================

plot_data <- genetic_physical_map %>%
  filter(
    !is.na(Haplotype),
    !is.na(RegenSY.Position.Mb),
    !is.na(Map.Position)
  )

if (nrow(plot_data) > 0) {

  genetic_physical_plot <- ggplot(
    plot_data,
    aes(
      x = RegenSY.Position.Mb,
      y = Map.Position,
      color = Haplotype
    )
  ) +
    geom_point(
      alpha = 0.55,
      size = 0.8
    ) +
    facet_wrap(
      ~ LG,
      scales = "free_x",
      ncol = 2
    ) +
    labs(
      x = "Physical position in RegenSY27× (Mb)",
      y = "Genetic position (cM)",
      color = "Haplotype"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )

  ggsave(
    filename = file.path(
      output_directory,
      "RegenSY_genetic_vs_physical_positions.png"
    ),
    plot = genetic_physical_plot,
    width = 10,
    height = 12,
    dpi = 300
  )

  ggsave(
    filename = file.path(
      output_directory,
      "RegenSY_genetic_vs_physical_positions.pdf"
    ),
    plot = genetic_physical_plot,
    width = 10,
    height = 12
  )
}

# ============================================================
# DISPLAY SUMMARY
# ============================================================

cat("\n============================================\n")
cat("MARKER POSITION PROJECTION COMPLETED\n")
cat("============================================\n\n")

print(mapping_summary)

cat("\nTotal markers in linkage map:",
    total_markers, "\n")

cat(
  "Markers with at least one valid RegenSY placement:",
  n_distinct(best_positions$Marker.Name),
  "\n"
)

cat(
  "Markers without any valid RegenSY placement:",
  nrow(markers_without_valid_position),
  "\n"
)

cat("\nOutputs are in:\n")
