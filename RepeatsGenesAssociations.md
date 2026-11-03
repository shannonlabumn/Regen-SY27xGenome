# Repeat and gene association in Regen-SY27x

Three R scripts quantify repeat overlap and nearest distances from complementary
perspectives. Upload this folder as `repeat_gene_association/` in your repository.

## Requirements and inputs

R packages: data.table, GenomicRanges, and IRanges. The scripts also call Unix
grep through data.table::fread. Install packages before execution.

Place inputs at the relative paths below, or edit the configuration blocks:

| Input path | Description |
| --- | --- |
| data/annotations/general_genes.gff3 | General gene annotation for scripts 01 and 02 |
| data/annotations/general_genes_helitron_analysis.gff3 | General gene annotation for script 03 |
| data/annotations/rgene_models.gff3 | R-gene source models, including gene features |
| data/annotations/rgene_classes.tsv | Headerless table with canonical R-gene IDs in column 1 |
| data/repeats/EDTA.intact.gff3 | Intact-repeat GFF3 with classification attributes |

All annotations must use the same assembly coordinates and sequence names.
The two general-annotation filenames preserve the different input GFFs in the
supplied scripts. They may point to the same annotation if appropriate; do not
silently substitute one annotation version for another when reproducing results.
GFF3 files must have nine tab-separated feature columns, no embedded FASTA
section, and the expected gene/mRNA/CDS Parent hierarchy. Multiple parents and
alternative transcript feature types require adapting the parser.

## Run

Run from this folder, with `data/` inside it:

```bash
Rscript 01_repeat_overlap_with_genes.R
Rscript 02_te_proximity_with_genes.R
# Set the verified total length of the assembly in base pairs before step 03:
export GENOME_SIZE_BP=$(awk '{total += $2} END {printf "%.0f\n", total}' data/genomes/assembly.fasta.fai)
Rscript 03_helitron_gene_proximity.R
```

The FASTA index in the example must cover the complete assembly represented by
the annotations, including all haplotypes and scaffolds used in the analysis.
Alternatively, set GENOME_SIZE_BP to its verified numeric total directly.
PROJECT_DIR can override the working directory for all scripts.

## Analyses and outputs

Outputs and cached masks are written to `results/repeat_gene_association/`.

| Script | Question | Main outputs |
| --- | --- | --- |
| 01_repeat_overlap_with_genes.R | What proportion of each gene span and coding footprint intersects repeats? | per_gene_genome_repeat_overlap.csv; genome_repeat_overlap_summary.txt |
| 02_te_proximity_with_genes.R | How far is each gene from the nearest repeat or Helitron? | per_gene_genome_te_proximity.csv; genome_te_proximity_summary.txt |
| 03_helitron_gene_proximity.R | How far is each reduced Helitron locus from a general gene or R gene? | per_helitron_gene_proximity.csv; helitron_gene_proximity_summary.txt |

Scripts can build the repeat masks from EDTA if caches are absent. The Helitron
classification is matched exactly as `DNA/Helitron`. Reduced masks merge both
overlapping and directly adjacent intervals, ignoring strand. They represent
non-redundant intervals rather than individual original EDTA calls. Existing
caches are not validated against the current inputs: remove them or choose a
fresh OUTDIR when changing the repeat annotation or assembly.

## Optional R-gene comparisons

These three scripts do not generate the following external tables. Missing
files cause only their comparison sections to be skipped.

| Optional path under data/rgene_comparisons/ | Required metric columns |
| --- | --- |
| rgene_final_confidence_per_gene.csv | pct_any_repeat_cds |
| rgene_confidence_per_gene.csv | pct_any_repeat |
| per_gene_te_proximity.csv | dist_to_any_repeat; dist_to_helitron |

Verify each table's gene membership, annotation version, assembly, repeat mask,
and metric definition before comparing it. The supplied whole-gene comparison
source used a smaller R-gene set than the canonical list; the released scripts
report the supplied table's actual row count rather than assuming membership.
General-annotation models overlapping canonical R-gene loci are excluded from
the comparison background by coordinates, since the annotations use different
ID namespaces. This exclusion does not establish statistical independence:
Wilcoxon tests remain unadjusted and do not model genomic spatial dependence.

## Metric definitions and limits

- GFF3 coordinates are one-based inclusive. Whole-gene spans include introns.
- Script 01 reduces coding intervals within each gene before measuring overlap.
  If CDS features are absent, it uses exons and warns; these include UTRs and
  must not be interpreted as equivalent to CDS-only measurements.
- Nearest-distance results follow GenomicRanges conventions. Zero includes
  overlapping and directly adjacent intervals; it is not an exact overlap count.
- Queries without targets on the same sequence retain NA. Proximity percentages
  use measurable queries as their denominator, reported alongside total counts.
- Distance summaries are cumulative at 1, 5, 10, and 20 kb.
- Script 03 retains MIN_CDS_LENGTH = 100. Its supplied implementation sums CDS
  rows across transcripts, which can count shared bases more than once. This is
  not the reduced coding-footprint length used by script 01 or a longest-isoform
  length. Set NULL to disable the filter, or adapt it for another definition.
- Duplicate gene rows are removed; multiple distinct loci sharing a gene ID may
  be reduced to the first occurrence. Review warnings and check ID uniqueness.
- The script 03 ratio compares the zero-distance proportion of measurable
  Helitron loci with general-gene-span genome coverage. This is descriptive:
  interval length, sequence availability, adjacency, and spatial structure mean
  it is not a calibrated random-overlap expectation or formal enrichment test.

## Public-release changes and validation

Personal paths were removed, scripts were numbered consistently, historical
comments were replaced with current step descriptions, and output directories
are created recursively. Core overlap and distance algorithms were retained.
The inconsistent hard-coded genome-size expression in script 03 was replaced
by the required GENOME_SIZE_BP setting. Distance-zero labels were clarified to
include adjacency; the distance calculations themselves were not changed.

Files were checked for personal paths and reviewed for input/output consistency.
R parsing and full execution were not performed because R and the genomic input
datasets were unavailable in the preparation environment.
