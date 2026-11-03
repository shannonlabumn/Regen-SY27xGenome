#!/usr/bin/env Rscript
# Workflow step 5/5: compare cumulative physical and genetic coordinates.
# Requires step 04 and the four FASTA indexes created by step 03.

# ============================================================
# RegenSY27x linkage-map versus physical-position figures
#
# Inputs:
#   RegenSY_marker_alignments/genetic_physical_map.tsv
#
#   RSYhap1noscf.fltd.fasta.fai
#   RSYhap2noscf.fltd.fasta.fai
#   RSYhap3noscf.fltd.fasta.fai
#   RSYhap4noscf.fltd.fasta.fai
#
# Main outputs:
#   1. One whole-genome figure per haplotype
#   2. One cumulative 2 x 2 grid for all haplotypes
#   3. Tables of chromosome-concordant and discordant markers
#
# Important:
#   Markers mapped to a different RegenSY chromosome are retained.
#   Their genetic coordinate is determined by the linkage group,
#   while their physical coordinate is determined by their
#   alignment to the RegenSY assembly.
# ============================================================

# ============================================================
# PACKAGES
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
map_file <- file.path(
  "RegenSY_marker_alignments",
  "genetic_physical_map.tsv"
)

haplotype_fai_files <- c(
  Hap1 = "data/genomes/RSYhap1noscf.fltd.fasta.fai",
  Hap2 = "data/genomes/RSYhap2noscf.fltd.fasta.fai",
  Hap3 = "data/genomes/RSYhap3noscf.fltd.fasta.fai",
  Hap4 = "data/genomes/RSYhap4noscf.fltd.fasta.fai"
)

output_directory <- file.path(
  "RegenSY_marker_alignments",
  "whole_genome_figures"
)

# Gap between consecutive physical chromosomes
chromosome_gap_mb <- 3

# Gap between consecutive linkage groups
genetic_gap_cm <- 10

# Figure settings
separate_figure_width <- 13
separate_figure_height <- 5.5

combined_figure_width <- 16
combined_figure_height <- 10

figure_dpi <- 400

point_size <- 0.9
point_alpha <- 0.70

expected_haplotypes <- c(
  "Hap1",
  "Hap2",
  "Hap3",
  "Hap4"
)

expected_chromosomes <- paste0(
  "Chr",
  1:8
)

placement_colors <- c(
  "Same chromosome" = "#355C9A",
  "Different chromosome" = "#D55E00",
  "Unknown" = "grey60"
)

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# FUNCTIONS
# ============================================================

normalize_chromosome <- function(x) {

  x <- str_trim(
    as.character(x)
  )

  ifelse(
    x %in% expected_chromosomes,
    x,
    NA_character_
  )
}

chromosome_number <- function(x) {

  suppressWarnings(
    as.integer(
      str_remove(
        as.character(x),
        "^Chr"
      )
    )
  )
}

read_fai <- function(
    fai_file,
    haplotype_name) {

  if (!file.exists(fai_file)) {
    stop(
      "FASTA index file was not found: ",
      fai_file
    )
  }

  read_tsv(
    file = fai_file,

    col_names = c(
      "RegenSY.Sequence",
      "Chromosome.Length",
      "Offset",
      "Line.Bases",
      "Line.Width"
    ),

    col_types = cols(
      RegenSY.Sequence = col_character(),
      Chromosome.Length = col_double(),
      Offset = col_double(),
      Line.Bases = col_double(),
      Line.Width = col_double()
    ),

    show_col_types = FALSE
  ) %>%
    transmute(
      Haplotype = haplotype_name,

      RegenSY.Sequence = str_trim(
        RegenSY.Sequence
      ),

      Chromosome.Length =
        as.numeric(Chromosome.Length)
    )
}

# ============================================================
# CHECK INPUT FILES
# ============================================================

if (!file.exists(map_file)) {
  stop(
    "Genetic-physical map was not found: ",
    map_file
  )
}

missing_fai_files <- haplotype_fai_files[
  !file.exists(haplotype_fai_files)
]

if (length(missing_fai_files) > 0) {
  stop(
    "Missing FASTA index files:\n",
    paste(
      missing_fai_files,
      collapse = "\n"
    )
  )
}

# ============================================================
# READ GENETIC-PHYSICAL MAP
# ============================================================

map <- read_tsv(
  map_file,
  show_col_types = FALSE
)

required_map_columns <- c(
  "Marker.Name",
  "LG",
  "Map.Position",
  "Haplotype",
  "RegenSY.Sequence",
  "RegenSY.Position"
)

missing_map_columns <- setdiff(
  required_map_columns,
  names(map)
)

if (length(missing_map_columns) > 0) {
  stop(
    "The genetic-physical map is missing columns: ",
    paste(
      missing_map_columns,
      collapse = ", "
    )
  )
}

cat(
  "\nRows read from genetic-physical map:",
  nrow(map),
  "\n"
)

# ============================================================
# READ FASTA INDEXES
# ============================================================

chromosome_lengths <- bind_rows(
  lapply(
    expected_haplotypes,
    function(haplotype_name) {

      read_fai(
        fai_file =
          haplotype_fai_files[[haplotype_name]],

        haplotype_name =
          haplotype_name
      )
    }
  )
)

# ============================================================
# CLEAN CHROMOSOME LENGTHS
# ============================================================

lengths_clean <- chromosome_lengths %>%
  mutate(
    Haplotype = str_trim(
      as.character(Haplotype)
    ),

    RegenSY.Sequence = str_trim(
      as.character(RegenSY.Sequence)
    ),

    Chromosome = normalize_chromosome(
      RegenSY.Sequence
    ),

    Chromosome.Number = chromosome_number(
      Chromosome
    ),

    Chromosome.Length =
      as.numeric(Chromosome.Length),

    Chromosome.Length.Mb =
      Chromosome.Length / 1e6
  ) %>%
  filter(
    Haplotype %in% expected_haplotypes,
    Chromosome %in% expected_chromosomes,
    Chromosome.Number %in% 1:8
  ) %>%
  arrange(
    match(
      Haplotype,
      expected_haplotypes
    ),
    Chromosome.Number
  )

# ============================================================
# VALIDATE CHROMOSOME LENGTHS
# ============================================================

duplicate_chromosomes <- lengths_clean %>%
  count(
    Haplotype,
    Chromosome,
    name = "Number.of.Records"
  ) %>%
  filter(
    Number.of.Records != 1
  )

if (nrow(duplicate_chromosomes) > 0) {

  print(duplicate_chromosomes)

  stop(
    "At least one chromosome has duplicate FASTA-index records."
  )
}

chromosome_length_check <- lengths_clean %>%
  count(
    Haplotype,
    name = "Number.of.Chromosomes"
  )

cat("\nChromosomes detected per haplotype:\n")
print(chromosome_length_check)

if (
  nrow(chromosome_length_check) != 4 ||
  any(
    chromosome_length_check$Number.of.Chromosomes != 8
  )
) {

  print(
    lengths_clean %>%
      select(
        Haplotype,
        RegenSY.Sequence,
        Chromosome,
        Chromosome.Number,
        Chromosome.Length.Mb
      ),
    n = Inf
  )

  stop(
    "Every haplotype must contain exactly eight chromosomes."
  )
}

chromosome_presence <- tidyr::expand_grid(
  Haplotype = expected_haplotypes,
  Chromosome = expected_chromosomes
) %>%
  dplyr::left_join(
    lengths_clean %>%
      dplyr::select(
        Haplotype,
        Chromosome,
        Chromosome.Length
      ),
    by = c(
      "Haplotype",
      "Chromosome"
    )
  )

missing_chromosomes <- chromosome_presence %>%
  filter(
    is.na(Chromosome.Length)
  )

if (nrow(missing_chromosomes) > 0) {

  print(missing_chromosomes)

  stop(
    "At least one expected haplotype/chromosome is missing."
  )
}

write_tsv(
  lengths_clean,
  file.path(
    output_directory,
    "chromosome_lengths_used_for_plotting.tsv"
  ),
  na = "NA"
)

# ============================================================
# CLEAN GENETIC-PHYSICAL MAP
# ============================================================

map_clean <- map %>%
  mutate(
    Marker.Name = str_trim(
      as.character(Marker.Name)
    ),

    Haplotype = str_trim(
      as.character(Haplotype)
    ),

    RegenSY.Sequence = str_trim(
      as.character(RegenSY.Sequence)
    ),

    Physical.Chromosome = normalize_chromosome(
      RegenSY.Sequence
    ),

    Physical.Chromosome.Number =
      chromosome_number(
        Physical.Chromosome
      ),

    LG = str_trim(
      as.character(LG)
    ),

    LG.Number = suppressWarnings(
      as.integer(LG)
    ),

    Genetic.Chromosome = if_else(
      !is.na(LG.Number) &
        LG.Number %in% 1:8,
      paste0("Chr", LG.Number),
      NA_character_
    ),

    Map.Position = suppressWarnings(
      as.numeric(Map.Position)
    ),

    RegenSY.Position = suppressWarnings(
      as.numeric(RegenSY.Position)
    )
  ) %>%
  filter(
    !is.na(Marker.Name),
    Marker.Name != "",
    Haplotype %in% expected_haplotypes,
    Physical.Chromosome %in% expected_chromosomes,
    Physical.Chromosome.Number %in% 1:8,
    LG.Number %in% 1:8,
    !is.na(Map.Position),
    !is.na(RegenSY.Position)
  ) %>%
  mutate(
    LG.Chromosome.Match =
      LG.Number ==
      Physical.Chromosome.Number,

    Marker.Placement = case_when(
      LG.Chromosome.Match ~
        "Same chromosome",

      !LG.Chromosome.Match ~
        "Different chromosome",

      TRUE ~
        "Unknown"
    ),

    Marker.Placement = factor(
      Marker.Placement,
      levels = c(
        "Same chromosome",
        "Different chromosome",
        "Unknown"
      )
    )
  )

cat(
  "\nUsable marker-haplotype records:",
  nrow(map_clean),
  "\n"
)

# ============================================================
# REPORT CHROMOSOME-CONCORDANT AND DISCORDANT MARKERS
# ============================================================

placement_summary <- map_clean %>%
  count(
    Haplotype,
    Marker.Placement,
    name = "Number.of.Markers"
  ) %>%
  arrange(
    match(
      Haplotype,
      expected_haplotypes
    ),
    Marker.Placement
  )

cat("\nMarker-placement summary:\n")
print(placement_summary, n = Inf)

write_tsv(
  placement_summary,
  file.path(
    output_directory,
    "marker_placement_summary.tsv"
  )
)

chromosome_concordant_markers <- map_clean %>%
  filter(
    LG.Chromosome.Match
  ) %>%
  arrange(
    Haplotype,
    LG.Number,
    Map.Position
  )

chromosome_discordant_markers <- map_clean %>%
  filter(
    !LG.Chromosome.Match
  ) %>%
  arrange(
    Haplotype,
    LG.Number,
    Physical.Chromosome.Number,
    Map.Position
  )

write_tsv(
  chromosome_concordant_markers,
  file.path(
    output_directory,
    "markers_mapped_to_same_chromosome.tsv"
  ),
  na = "NA"
)

write_tsv(
  chromosome_discordant_markers,
  file.path(
    output_directory,
    "markers_mapped_to_different_chromosomes.tsv"
  ),
  na = "NA"
)

cat(
  "\nChromosome-discordant marker records:",
  nrow(chromosome_discordant_markers),
  "\n"
)

# ============================================================
# PHYSICAL CHROMOSOME LAYOUT
# ============================================================

chromosome_layout <- lengths_clean %>%
  group_by(
    Haplotype
  ) %>%
  arrange(
    Chromosome.Number,
    .by_group = TRUE
  ) %>%
  mutate(
    Chromosome.Offset.Mb = lag(
      cumsum(
        Chromosome.Length.Mb +
          chromosome_gap_mb
      ),
      default = 0
    ),

    Chromosome.Start.Mb =
      Chromosome.Offset.Mb,

    Chromosome.End.Mb =
      Chromosome.Offset.Mb +
      Chromosome.Length.Mb,

    Chromosome.Midpoint.Mb =
      Chromosome.Start.Mb +
      Chromosome.Length.Mb / 2
  ) %>%
  ungroup()

write_tsv(
  chromosome_layout,
  file.path(
    output_directory,
    "cumulative_physical_chromosome_layout.tsv"
  ),
  na = "NA"
)

cat("\nPhysical chromosome layout:\n")

print(
  chromosome_layout %>%
    select(
      Haplotype,
      Chromosome,
      Chromosome.Number,
      Chromosome.Length.Mb,
      Chromosome.Start.Mb,
      Chromosome.End.Mb
    ),
  n = Inf
)

# ============================================================
# JOIN MARKERS TO PHYSICAL CHROMOSOME LAYOUT
# ============================================================

plot_data <- map_clean %>%
  left_join(
    chromosome_layout %>%
      select(
        Haplotype,
        Chromosome,
        Chromosome.Length.Mb,
        Chromosome.Offset.Mb,
        Chromosome.Start.Mb,
        Chromosome.End.Mb,
        Chromosome.Midpoint.Mb
      ),

    by = c(
      "Haplotype",
      "Physical.Chromosome" = "Chromosome"
    )
  ) %>%
  mutate(
    Physical.Position.Mb =
      RegenSY.Position / 1e6,

    Genome.Position.Mb =
      Chromosome.Offset.Mb +
      Physical.Position.Mb
  )

missing_physical_layout <- plot_data %>%
  filter(
    is.na(Chromosome.Offset.Mb)
  )

if (nrow(missing_physical_layout) > 0) {

  write_tsv(
    missing_physical_layout,
    file.path(
      output_directory,
      "markers_missing_physical_layout.tsv"
    ),
    na = "NA"
  )

  stop(
    nrow(missing_physical_layout),
    " markers failed to join to the physical chromosome layout."
  )
}

positions_outside_chromosomes <- plot_data %>%
  filter(
    Physical.Position.Mb < 0 |
      Physical.Position.Mb >
      Chromosome.Length.Mb
  )

if (nrow(positions_outside_chromosomes) > 0) {

  write_tsv(
    positions_outside_chromosomes,
    file.path(
      output_directory,
      "positions_outside_chromosome_lengths.tsv"
    ),
    na = "NA"
  )

  stop(
    nrow(positions_outside_chromosomes),
    " markers have coordinates outside their physical chromosome."
  )
}

# ============================================================
# GENETIC LINKAGE-GROUP LAYOUT
# ============================================================

genetic_layout <- map_clean %>%
  select(
    Marker.Name,
    LG.Number,
    Genetic.Chromosome,
    Map.Position
  ) %>%
  distinct() %>%
  group_by(
    LG.Number,
    Genetic.Chromosome
  ) %>%
  summarise(
    Genetic.Minimum.cM =
      min(Map.Position, na.rm = TRUE),

    Genetic.Maximum.cM =
      max(Map.Position, na.rm = TRUE),

    Genetic.Length.cM =
      Genetic.Maximum.cM -
      Genetic.Minimum.cM,

    Number.of.Markers =
      n_distinct(Marker.Name),

    .groups = "drop"
  ) %>%
  arrange(
    LG.Number
  ) %>%
  mutate(
    Genetic.Offset.cM = lag(
      cumsum(
        Genetic.Length.cM +
          genetic_gap_cm
      ),
      default = 0
    ),

    Genetic.Start.cM =
      Genetic.Offset.cM,

    Genetic.End.cM =
      Genetic.Offset.cM +
      Genetic.Length.cM,

    Genetic.Midpoint.cM =
      Genetic.Start.cM +
      Genetic.Length.cM / 2
  )

if (nrow(genetic_layout) != 8) {

  print(genetic_layout, n = Inf)

  stop(
    "The genetic layout must contain exactly eight linkage groups."
  )
}

write_tsv(
  genetic_layout,
  file.path(
    output_directory,
    "cumulative_genetic_map_layout.tsv"
  ),
  na = "NA"
)

cat("\nGenetic linkage-group layout:\n")
print(genetic_layout, n = Inf)

# ============================================================
# ADD CUMULATIVE GENETIC COORDINATES
# ============================================================

cumulative_plot_data <- plot_data %>%
  left_join(
    genetic_layout %>%
      select(
        LG.Number,
        Genetic.Minimum.cM,
        Genetic.Maximum.cM,
        Genetic.Length.cM,
        Genetic.Offset.cM,
        Genetic.Start.cM,
        Genetic.End.cM,
        Genetic.Midpoint.cM
      ),

    by = "LG.Number"
  ) %>%
  mutate(
    Genetic.Position.Cumulative.cM =
      Genetic.Offset.cM +
      Map.Position -
      Genetic.Minimum.cM
  )

missing_genetic_layout <- cumulative_plot_data %>%
  filter(
    is.na(Genetic.Offset.cM) |
      is.na(Genetic.Position.Cumulative.cM)
  )

if (nrow(missing_genetic_layout) > 0) {

  write_tsv(
    missing_genetic_layout,
    file.path(
      output_directory,
      "markers_missing_genetic_layout.tsv"
    ),
    na = "NA"
  )

  stop(
    nrow(missing_genetic_layout),
    " markers failed to join to the genetic layout."
  )
}

# ============================================================
# MARKER COUNTS BY PHYSICAL CHROMOSOME
# ============================================================

observed_marker_counts <- plot_data %>%
  count(
    Haplotype,
    Physical.Chromosome,
    Physical.Chromosome.Number,
    name = "Number.of.Markers"
  )

marker_counts <- chromosome_layout %>%
  select(
    Haplotype,
    Physical.Chromosome = Chromosome,
    Physical.Chromosome.Number = Chromosome.Number
  ) %>%
  left_join(
    observed_marker_counts,
    by = c(
      "Haplotype",
      "Physical.Chromosome",
      "Physical.Chromosome.Number"
    )
  ) %>%
  mutate(
    Number.of.Markers = replace_na(
      Number.of.Markers,
      0L
    )
  ) %>%
  arrange(
    match(
      Haplotype,
      expected_haplotypes
    ),
    Physical.Chromosome.Number
  )

if (nrow(marker_counts) != 32) {
  stop(
    "Marker-count table must contain 32 rows."
  )
}

write_tsv(
  marker_counts,
  file.path(
    output_directory,
    "marker_counts_by_haplotype_physical_chromosome.tsv"
  )
)

cat("\nMarker counts by physical chromosome:\n")
print(marker_counts, n = Inf)

# ============================================================
# SAVE PLOTTING DATA
# ============================================================

write_tsv(
  plot_data,
  file.path(
    output_directory,
    "whole_genome_plotting_data.tsv"
  ),
  na = "NA"
)

write_tsv(
  cumulative_plot_data,
  file.path(
    output_directory,
    "cumulative_genetic_physical_plotting_data.tsv"
  ),
  na = "NA"
)

# ============================================================
# SEPARATE FIGURE FUNCTION
#
# x = cumulative physical position
# y = within-linkage-group cM position
# ============================================================

plot_one_haplotype <- function(haplotype_name) {

  hap_data <- plot_data %>%
    filter(
      Haplotype == haplotype_name
    )

  hap_layout <- chromosome_layout %>%
    filter(
      Haplotype == haplotype_name
    ) %>%
    arrange(
      Chromosome.Number
    )

  if (nrow(hap_layout) != 8) {
    stop(
      haplotype_name,
      " does not have eight chromosome-layout records."
    )
  }

  if (nrow(hap_data) == 0) {
    stop(
      "No markers are available for ",
      haplotype_name
    )
  }

  physical_boundaries <- hap_layout[
    2:nrow(hap_layout),
    ,
    drop = FALSE
  ]

  hap_x_max <- max(
    hap_layout$Chromosome.End.Mb,
    na.rm = TRUE
  ) * 1.04

  ggplot(
    hap_data,
    aes(
      x = Genome.Position.Mb,
      y = Map.Position,
      color = Marker.Placement
    )
  ) +

    geom_vline(
      data = physical_boundaries,
      aes(
        xintercept = Chromosome.Start.Mb
      ),
      inherit.aes = FALSE,
      color = "grey60",
      linewidth = 0.4
    ) +

    geom_point(
      size = point_size,
      alpha = point_alpha
    ) +

    scale_color_manual(
      values = placement_colors,
      breaks = c(
        "Same chromosome",
        "Different chromosome"
      ),
      name = "Marker placement",
      drop = FALSE
    ) +

    scale_x_continuous(
      breaks =
        hap_layout$Chromosome.Midpoint.Mb,

      labels =
        hap_layout$Chromosome,

      limits = c(
        0,
        hap_x_max
      ),

      expand = expansion(
        mult = c(0.002, 0.005)
      )
    ) +

    scale_y_continuous(
      expand = expansion(
        mult = c(0.02, 0.05)
      )
    ) +

    labs(
      title = paste0(
        haplotype_name,
        ": linkage-map and RegenSY27x physical positions"
      ),

      subtitle = paste0(
        "Orange points indicate markers assigned to a different ",
        "linkage-group chromosome"
      ),

      x = "RegenSY27x chromosome",

      y = "Genetic position (cM)"
    ) +

    theme_bw(
      base_size = 12
    ) +

    theme(
      panel.grid.minor =
        element_blank(),

      panel.grid.major.x =
        element_blank(),

      axis.text.x =
        element_text(
          face = "bold",
          size = 11
        ),

      axis.title =
        element_text(
          face = "bold"
        ),

      plot.title =
        element_text(
          face = "bold",
          size = 16
        ),

      plot.subtitle =
        element_text(
          size = 9
        ),

      legend.position =
        "bottom"
    )
}

# ============================================================
# GENERATE SEPARATE HAPLOTYPE FIGURES
# ============================================================

haplotype_plots <- list()

for (hap in expected_haplotypes) {

  message(
    "Generating separate figure for ",
    hap
  )

  p <- plot_one_haplotype(
    hap
  )

  haplotype_plots[[hap]] <- p

  ggsave(
    filename = file.path(
      output_directory,
      paste0(
        hap,
        "_whole_genome_genetic_physical_map.png"
      )
    ),
    plot = p,
    width = separate_figure_width,
    height = separate_figure_height,
    units = "in",
    dpi = figure_dpi,
    bg = "white"
  )

  ggsave(
    filename = file.path(
      output_directory,
      paste0(
        hap,
        "_whole_genome_genetic_physical_map.pdf"
      )
    ),
    plot = p,
    width = separate_figure_width,
    height = separate_figure_height,
    units = "in",
    device = "pdf"
  )
}

# ============================================================
# DATA FOR CUMULATIVE FOUR-HAPLOTYPE GRID
# ============================================================

physical_boundaries_all <- chromosome_layout %>%
  filter(
    Chromosome.Number > 1
  )

genetic_boundaries <- genetic_layout %>%
  filter(
    LG.Number > 1
  )

chromosome_labels <- chromosome_layout %>%
  select(
    Haplotype,
    Chromosome,
    Chromosome.Number,
    Chromosome.Midpoint.Mb
  )

global_physical_max <- max(
  chromosome_layout$Chromosome.End.Mb,
  na.rm = TRUE
)

# Extra room prevents Chr8 labels from being clipped
physical_axis_max <- global_physical_max * 1.04

global_genetic_max <- max(
  genetic_layout$Genetic.End.cM,
  na.rm = TRUE
)

genetic_axis_max <- global_genetic_max * 1.02

# ============================================================
# CUMULATIVE 2 x 2 GRID
#
# x = cumulative RegenSY physical position
# y = cumulative linkage-map position
# ============================================================

haplotype_labels <- c(
  Hap1 = "Haplotype 1",
  Hap2 = "Haplotype 2",
  Hap3 = "Haplotype 3",
  Hap4 = "Haplotype 4"
)

cumulative_map_figure <- ggplot(
  cumulative_plot_data,
  aes(
    x = Genome.Position.Mb,
    y = Genetic.Position.Cumulative.cM,
    color = Marker.Placement
  )
) +

  geom_vline(
    data = physical_boundaries_all,
    aes(
      xintercept = Chromosome.Start.Mb
    ),
    inherit.aes = FALSE,
    color = "grey65",
    linewidth = 0.35
  ) +

  geom_hline(
    data = genetic_boundaries,
    aes(
      yintercept = Genetic.Start.cM
    ),
    inherit.aes = FALSE,
    color = "grey65",
    linewidth = 0.35
  ) +

  geom_point(
    size = point_size,
    alpha = point_alpha
  ) +

  geom_text(
    data = chromosome_labels,
    aes(
      x = Chromosome.Midpoint.Mb,
      y = Inf,
      label = Chromosome
    ),
    inherit.aes = FALSE,
    vjust = 1.25,
    hjust = 0.5,
    fontface = "bold",
    size = 3.2,
    color = "grey20"
  ) +

  facet_wrap(
    ~ Haplotype,
    ncol = 2,
    scales = "fixed",
    labeller = labeller(
      Haplotype = haplotype_labels
    )
  ) +

  scale_color_manual(
    values = placement_colors,
    breaks = c(
      "Same chromosome",
      "Different chromosome"
    ),
    name = "Marker placement",
    drop = FALSE
  ) +

  scale_x_continuous(
    limits = c(
      0,
      physical_axis_max
    ),

    expand = expansion(
      mult = c(0.002, 0.005)
    )
  ) +

  scale_y_continuous(
    limits = c(
      0,
      genetic_axis_max
    ),

    expand = expansion(
      mult = c(0.002, 0.01)
    )
  ) +

  labs(
    # title = paste0(
    #   "Collinearity between the linkage map and ",
    #   "RegenSY27x haplotype assemblies"
    # ),
    #
    # subtitle = paste0(
    #   "Physical and genetic positions are cumulative; ",
    #   "orange points indicate chromosome-discordant placements"
    # ),

    x = "Cumulative physical position in RegenSY27x (Mb)",

    y = "Cumulative genetic-map position (cM)"
  ) +

  theme_bw(
    base_size = 12
  ) +

  theme(
    panel.grid.minor =
      element_blank(),

    panel.grid.major =
      element_blank(),

    strip.background =
      element_rect(
        fill = "grey90",
        color = "grey40",
        linewidth = 0.5
      ),

    strip.text =
      element_text(
        face = "bold",
        size = 12
      ),

    axis.title =
      element_text(
        face = "bold"
      ),

    axis.text =
      element_text(
        size = 9
      ),

    plot.title =
      element_text(
        face = "bold",
        size = 16
      ),

    plot.subtitle =
      element_text(
        size = 10
      ),

    legend.position =
      "bottom",

    panel.spacing =
      grid::unit(
        1,
        "lines"
      )
  )

# ============================================================
# SAVE CUMULATIVE GRID
# ============================================================

ggsave(
  filename = file.path(
    output_directory,
    "all_haplotypes_cumulative_genetic_physical_map.png"
  ),
  plot = cumulative_map_figure,
  width = combined_figure_width,
  height = combined_figure_height,
  units = "in",
  dpi = figure_dpi,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "all_haplotypes_cumulative_genetic_physical_map.pdf"
  ),
  plot = cumulative_map_figure,
  width = combined_figure_width,
  height = combined_figure_height,
  units = "in",
  device = "pdf"
)

# ============================================================
# FINAL FILE CHECK
# ============================================================

expected_separate_png <- file.path(
  output_directory,
  paste0(
    expected_haplotypes,
    "_whole_genome_genetic_physical_map.png"
  )
)

expected_separate_pdf <- file.path(
  output_directory,
  paste0(
    expected_haplotypes,
    "_whole_genome_genetic_physical_map.pdf"
  )
)

combined_png <- file.path(
  output_directory,
  "all_haplotypes_cumulative_genetic_physical_map.png"
)

combined_pdf <- file.path(
  output_directory,
  "all_haplotypes_cumulative_genetic_physical_map.pdf"
)

expected_files <- c(
  expected_separate_png,
  expected_separate_pdf,
  combined_png,
  combined_pdf
)

if (!all(file.exists(expected_files))) {

  missing_outputs <- expected_files[
    !file.exists(expected_files)
  ]

  stop(
    "Expected figure files were not created:\n",
    paste(
      missing_outputs,
      collapse = "\n"
    )
  )
}

# ============================================================
# FINAL REPORT
# ============================================================

cat("\n============================================\n")
cat("SCRIPT 05 COMPLETED SUCCESSFULLY\n")
cat("============================================\n\n")

cat("Output directory:\n")
cat(output_directory, "\n\n")

cat("Marker placement summary:\n")
print(placement_summary, n = Inf)

cat("\nSeparate haplotype figures:\n")
cat(
  paste(
    expected_separate_png,
    collapse = "\n"
  ),
  "\n"
)

cat("\nCombined cumulative figure:\n")
cat(combined_png, "\n")
cat(combined_pdf, "\n")

cat("\nDiscordant-marker table:\n")

cat(
  file.path(
    output_directory,
    "markers_mapped_to_different_chromosomes.tsv"
  ),
  "\n"
)