# Linkage-map versus genome-assembly comparison

Five numbered scripts prepare linkage-map SNPs, extract their flanking sequences,
align them to four Regen-SY27x haplotype assemblies, project SNP positions through
CIGAR strings, and create genetic–physical collinearity (Marey) plots.

## Requirements

- R with dplyr, readr, stringr, tidyr, and ggplot2.
- Bash, BWA with BWA-MEM support, samtools, bedtools, and standard Unix tools.
- The original workflow specified BWA 0.7.4, samtools 1.21, and bedtools 2.29.2.
  Confirm that the installed BWA provides `bwa mem`; version labels and local
  builds should be verified before reporting the software used.
- Activate software before execution. SLURM directives are examples; adjust
  resource requests to your cluster. Scripts also run directly with Bash.

## Project layout

Run from the directory containing these scripts. Place the linkage-map CSV at
`data/linkage_map.csv`, the original genotyping reference at
`data/genomes/original_reference.fasta`, and the four target assemblies at
`data/genomes/RSYhap1noscf.fltd.fasta` through
`data/genomes/RSYhap4noscf.fltd.fasta`.

The CSV must contain Marker.Name, LG, Ref.Chrom, Ref.Position, Ref.Allele.x,
Alt.Allele.x, and Map.Position. Ref.Position must be a one-based coordinate on
exactly the same reference assembly/version used for sequence extraction.
Configure chromosome_name_mode and chromosome_prefix in step 01 to match the
original FASTA sequence names. Review the chromosome-name parsing in steps 04
and 05 if your target assemblies use another naming convention.

## Run order

```bash
Rscript 01_prepare_linkage_markers.R
bash 02_extract_linkage_marker_sequences.sh
bash 03_align_markers_to_RegenSY_haplotypes.sh
Rscript 04_filter_and_project_RegenSY_positions.R
Rscript 05_whole_genome_marey_plots.R
```

Wait for each step to finish successfully before starting the next. On SLURM,
submit steps 02 and 03 using sbatch instead of bash and wait for job completion.
Do not submit dependent steps simultaneously without job dependencies.

PROJECT_DIR optionally sets the working directory for every script. MAP_FILE
sets the CSV path in step 01; ORIGINAL_REFERENCE sets the source FASTA in step
02; HAP1 through HAP4 set target FASTA paths in step 03. If target FASTA paths
are changed, update haplotype_fai_files in step 05 to match their .fai paths.
Relative paths are resolved within PROJECT_DIR. Default outputs are written in
that directory and RegenSY_marker_alignments/. Rerunning can overwrite outputs;
use a separate project directory to preserve a previous run.

## What each script does

| Script | Main inputs | Main outputs |
| --- | --- | --- |
| 01 | Linkage-map CSV | Marker metadata, one-base SNP BED, marker summary |
| 02 | SNP BED and original reference | Clipped flank BED, flank FASTA, length and extraction QC |
| 03 | Flank FASTA and four assemblies | Sorted BAMs, alignment QC, combined primary positions, FASTA indexes |
| 04 | Metadata, flank BED, primary positions | Alignment QC, rejected records, projected SNP positions, genetic_physical_map.tsv |
| 05 | Genetic–physical table and four FASTA indexes | Per-haplotype and combined figures, plotting data, concordance tables |

## Coordinates and interpretation

BED coordinates are zero-based, half-open; reference SNP positions and SAM
alignment starts are one-based. Step 02 extracts 250 bases on each side of the
SNP (501 bases away from chromosome ends), clipping at sequence boundaries.
Step 04 uses Ref.Position minus the extracted BED start to locate the SNP in
its actual query sequence, including shortened flanks. It then accounts for
strand orientation and CIGAR operations to project the SNP onto each assembly.
This corrects the supplied script's midpoint assumption for shortened flanks.

Step 04 retains the supplied thresholds: MAPQ >=20, query coverage >=0.80, and
require_matching_chromosome = FALSE. Primary mapped records are evaluated;
MAPQ filtering does not prove a marker has only one genomic match. Inspect the
QC outputs before interpreting discordance as a biological rearrangement.

Step 05 retains chromosome-discordant markers and uses cumulative physical and
genetic axes with display gaps. The plots describe collinearity; they are not a
separate statistical test of rearrangements or a formal recombination-rate fit.

## Changes made for public release

Personal filesystem paths and email addresses were removed. Inputs use a shared
portable layout, environment-specific module loads are commented examples, and
step comments explain inputs, filters, coordinate conventions, and outputs.
Step 03 uses portable arithmetic to inspect the reverse-strand flag and splits
pipeline threads across BWA and sorting. Step 04 corrects boundary-flank SNP
offsets as described above. Other analytical thresholds and plotting settings
are retained from the supplied scripts.

## Validation

The Bash scripts passed syntax checks. Full execution and R parsing were not
performed here because R and the input genomic datasets were unavailable.
