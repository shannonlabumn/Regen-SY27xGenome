#!/bin/bash -l

#SBATCH --job-name=extract_map_markers
#SBATCH --output=extract_map_markers_%j.out
#SBATCH --error=extract_map_markers_%j.err
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

# Workflow step 2/5: extract marker flanks from the original reference.
# Activate samtools and bedtools before execution; module examples are optional.
set -euo pipefail
# Keep intermediate filenames consistent across all five steps.
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
cd -- "$PROJECT_DIR"

# ============================================================
# Extract ±250-bp flanking sequences around linkage-map SNPs
#
# Inputs:
#   linkage_markers.snp.bed
#   Original reference FASTA used to define Ref.Position
#
# Outputs:
#   linkage_markers.flank250.bed
#   linkage_markers.flank250.fa
#   linkage_marker_sequence_lengths.tsv
#   linkage_marker_extraction_summary.txt
# ============================================================

# ============================================================
# USER SETTINGS
# ============================================================

# IMPORTANT:
# This must be the same reference assembly/version used to
# generate Ref.Chrom and Ref.Position in the linkage map.
ORIGINAL_REFERENCE="${ORIGINAL_REFERENCE:-data/genomes/original_reference.fasta}"

SNP_BED="linkage_markers.snp.bed"

FLANK_SIZE=250

OUTPUT_BED="linkage_markers.flank250.bed"
OUTPUT_FASTA="linkage_markers.flank250.fa"

# ============================================================
# LOAD SOFTWARE
# ============================================================

# Modify these module names if your cluster uses different names.
 # Example environment module: module load samtools/1.21
# Example environment module: module load bedtools/2.29.2

# ============================================================
# CHECK SOFTWARE
# ============================================================

if ! command -v samtools > /dev/null 2>&1
then
    echo "ERROR: samtools was not found." >&2
    exit 1
fi

if ! command -v bedtools > /dev/null 2>&1
then
    echo "ERROR: bedtools was not found." >&2
    exit 1
fi

echo "samtools:"
samtools --version | head -n 1

echo "bedtools:"
bedtools --version

# ============================================================
# CHECK INPUT FILES
# ============================================================

if [[ ! -s "${ORIGINAL_REFERENCE}" ]]
then
    echo "ERROR: Original reference FASTA not found or empty:" >&2
    echo "${ORIGINAL_REFERENCE}" >&2
    exit 1
fi

if [[ ! -s "${SNP_BED}" ]]
then
    echo "ERROR: SNP BED file not found or empty:" >&2
    echo "${SNP_BED}" >&2
    exit 1
fi

# ============================================================
# VALIDATE BED FORMAT
# ============================================================

awk '
BEGIN {
    error_count = 0
}
NF < 4 {
    print "Invalid BED row with fewer than four columns at line " NR > "/dev/stderr"
    error_count++
}
$2 !~ /^[0-9]+$/ {
    print "Invalid BED start coordinate at line " NR ": " $2 > "/dev/stderr"
    error_count++
}
$3 !~ /^[0-9]+$/ {
    print "Invalid BED end coordinate at line " NR ": " $3 > "/dev/stderr"
    error_count++
}
$2 >= $3 {
    print "BED start is not smaller than end at line " NR > "/dev/stderr"
    error_count++
}
END {
    if (error_count > 0) {
        exit 1
    }
}
' "${SNP_BED}"

# ============================================================
# INDEX ORIGINAL REFERENCE
# ============================================================

if [[ ! -s "${ORIGINAL_REFERENCE}.fai" ]]
then
    echo "Creating FASTA index..."

    samtools faidx "${ORIGINAL_REFERENCE}"
else
    echo "Existing FASTA index found:"
    echo "${ORIGINAL_REFERENCE}.fai"
fi

# Create bedtools genome file:
# chromosome<TAB>chromosome_length
cut -f1,2 "${ORIGINAL_REFERENCE}.fai" \
    > original_reference.genome

# ============================================================
# CHECK FOR DUPLICATED FASTA SEQUENCE NAMES
# ============================================================

cut -f1 "${ORIGINAL_REFERENCE}.fai" |
    sort |
    uniq -d \
    > duplicated_reference_sequence_names.txt

if [[ -s duplicated_reference_sequence_names.txt ]]
then
    echo "ERROR: Duplicated reference sequence names were found:" >&2
    cat duplicated_reference_sequence_names.txt >&2
    exit 1
fi

# ============================================================
# CHECK BED CHROMOSOME NAMES AGAINST FASTA
# ============================================================

cut -f1 "${SNP_BED}" |
    sort -u \
    > bed_chromosomes.txt

cut -f1 original_reference.genome |
    sort -u \
    > fasta_chromosomes.txt

comm -23 \
    bed_chromosomes.txt \
    fasta_chromosomes.txt \
    > chromosomes_missing_from_reference.txt

if [[ -s chromosomes_missing_from_reference.txt ]]
then
    echo "ERROR: These BED chromosome names are absent from the FASTA:" >&2
    cat chromosomes_missing_from_reference.txt >&2

    echo "" >&2
    echo "Chromosome names in the BED file:" >&2
    cat bed_chromosomes.txt >&2

    echo "" >&2
    echo "First chromosome names in the FASTA:" >&2
    head -n 30 fasta_chromosomes.txt >&2

    echo "" >&2
    echo "Change chromosome_name_mode in the R script." >&2

    exit 1
fi

echo "All BED chromosome names were found in the reference FASTA."

# ============================================================
# SORT ONE-BASE SNP BED FILE
# ============================================================

bedtools sort \
    -i "${SNP_BED}" \
    -g original_reference.genome \
    > linkage_markers.snp.sorted.bed

# ============================================================
# EXTEND EACH SNP BY ±250 BP
# ============================================================
#
# Input SNP interval:
#   one base
#
# With -b 250:
#   total expected length = 501 bp
#
# bedtools clips intervals at chromosome boundaries.
# ============================================================

bedtools slop \
    -i linkage_markers.snp.sorted.bed \
    -g original_reference.genome \
    -b "${FLANK_SIZE}" \
    > "${OUTPUT_BED}"

# ============================================================
# CHECK FOR DUPLICATED MARKER NAMES
# ============================================================

cut -f4 "${OUTPUT_BED}" |
    sort |
    uniq -d \
    > duplicated_marker_names_in_bed.txt

if [[ -s duplicated_marker_names_in_bed.txt ]]
then
    echo "ERROR: Duplicated marker names were found in the BED file:" >&2
    head duplicated_marker_names_in_bed.txt >&2
    exit 1
fi

# ============================================================
# EXTRACT MARKER SEQUENCES
# ============================================================
#
# -nameOnly:
#   FASTA headers contain Marker.Name only.
#
# Example:
#   >Chr1_745191
# ============================================================

bedtools getfasta \
    -fi "${ORIGINAL_REFERENCE}" \
    -bed "${OUTPUT_BED}" \
    -nameOnly \
    -fo "${OUTPUT_FASTA}"

# ============================================================
# COUNT BED AND FASTA RECORDS
# ============================================================

BED_COUNT=$(wc -l < "${OUTPUT_BED}")
FASTA_COUNT=$(grep -c '^>' "${OUTPUT_FASTA}")

echo "BED marker count:   ${BED_COUNT}"
echo "FASTA record count: ${FASTA_COUNT}"

if [[ "${BED_COUNT}" -ne "${FASTA_COUNT}" ]]
then
    echo "ERROR: BED and FASTA record counts differ." >&2
    exit 1
fi

# ============================================================
# CALCULATE FASTA SEQUENCE LENGTHS
# ============================================================

awk '
BEGIN {
    OFS = "\t"
}
function print_record() {
    if (header != "") {
        print header, length(sequence)
    }
}
/^>/ {
    print_record()
    header = substr($0, 2)
    sequence = ""
    next
}
{
    gsub(/[[:space:]]/, "", $0)
    sequence = sequence $0
}
END {
    print_record()
}
' "${OUTPUT_FASTA}" \
> linkage_marker_sequence_lengths.tsv

# ============================================================
# IDENTIFY NONSTANDARD SEQUENCE LENGTHS
# ============================================================

EXPECTED_LENGTH=$((2 * FLANK_SIZE + 1))

awk -v expected="${EXPECTED_LENGTH}" '
BEGIN {
    OFS = "\t"
}
$2 != expected {
    print $1, $2
}
' linkage_marker_sequence_lengths.tsv \
> markers_with_shortened_flanks.tsv

SHORTENED_COUNT=$(wc -l < markers_with_shortened_flanks.tsv)

# Shortened sequences should generally only represent markers
# within 250 bp of the beginning or end of a chromosome.
if [[ "${SHORTENED_COUNT}" -gt 0 ]]
then
    echo "NOTICE: ${SHORTENED_COUNT} markers have sequences shorter than ${EXPECTED_LENGTH} bp."
    echo "These may be close to chromosome ends."
fi

# ============================================================
# CHECK FOR EMPTY SEQUENCES
# ============================================================

EMPTY_COUNT=$(
    awk '$2 == 0 {count++} END {print count+0}' \
        linkage_marker_sequence_lengths.tsv
)

if [[ "${EMPTY_COUNT}" -gt 0 ]]
then
    echo "ERROR: ${EMPTY_COUNT} empty sequences were extracted." >&2
    exit 1
fi

# ============================================================
# CREATE EXTRACTION SUMMARY
# ============================================================

{
    echo "Linkage-marker sequence extraction"
    echo
    echo "Original reference: ${ORIGINAL_REFERENCE}"
    echo "Input SNP BED: ${SNP_BED}"
    echo "Flank size per side: ${FLANK_SIZE} bp"
    echo "Expected sequence length: ${EXPECTED_LENGTH} bp"
    echo "Markers in final BED: ${BED_COUNT}"
    echo "Sequences in final FASTA: ${FASTA_COUNT}"
    echo "Markers with shortened flanks: ${SHORTENED_COUNT}"
    echo "Empty sequences: ${EMPTY_COUNT}"
    echo
    echo "Output BED: ${OUTPUT_BED}"
    echo "Output FASTA: ${OUTPUT_FASTA}"
} > linkage_marker_extraction_summary.txt

# ============================================================
# DISPLAY FINAL SUMMARY
# ============================================================

echo
echo "============================================"
echo "MARKER SEQUENCE EXTRACTION COMPLETED"
echo "============================================"
echo

cat linkage_marker_extraction_summary.txt

echo
echo "Sequence-length distribution:"

cut -f2 linkage_marker_sequence_lengths.tsv |
    sort -n |
    uniq -c

echo
echo "First five BED records:"
head -n 5 "${OUTPUT_BED}"

echo
echo "First FASTA headers:"
grep '^>' "${OUTPUT_FASTA}" |
    head -n 5

echo
echo "Files created:"
echo "  ${OUTPUT_BED}"
echo "  ${OUTPUT_FASTA}"
echo "  linkage_marker_sequence_lengths.tsv"
echo "  linkage_marker_extraction_summary.txt"
echo "  markers_with_shortened_flanks.tsv"
