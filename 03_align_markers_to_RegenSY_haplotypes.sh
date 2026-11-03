#!/bin/bash -l
#SBATCH --job-name=map_markers_RegenSY
#SBATCH --output=map_markers_RegenSY_%j.out
#SBATCH --error=map_markers_RegenSY_%j.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G

# Workflow step 3/5: align the marker flanks independently to four assemblies.
# Requirements: BWA with the mem command, samtools, and awk on PATH.
# Activate software before running; cluster module names are site-specific.
# Example: module load bwa/0.7.4 samtools/1.21
set -euo pipefail

# 1. Set project paths. All five scripts share the same working directory.
PROJECT_DIR="${PROJECT_DIR:-$PWD}"
cd -- "$PROJECT_DIR"
HAP1="${HAP1:-data/genomes/RSYhap1noscf.fltd.fasta}"
HAP2="${HAP2:-data/genomes/RSYhap2noscf.fltd.fasta}"
HAP3="${HAP3:-data/genomes/RSYhap3noscf.fltd.fasta}"
HAP4="${HAP4:-data/genomes/RSYhap4noscf.fltd.fasta}"
MARKERS="${MARKERS:-linkage_markers.flank250.fa}"
OUTDIR="RegenSY_marker_alignments"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
# BWA and samtools sort run in a pipeline and share the CPU allocation.
SORT_THREADS=1
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo 'Invalid THREADS' >&2; exit 1; }
BWA_THREADS=$((THREADS > 2 ? THREADS - 2 : 1))

# 2. Fail before alignment if a dependency or input is missing.
for program in bwa samtools awk; do
    command -v "$program" >/dev/null 2>&1 || { echo "Missing executable: $program" >&2; exit 1; }
done
for file in "$HAP1" "$HAP2" "$HAP3" "$HAP4" "$MARKERS"; do
    [[ -s "$file" ]] || { echo "Missing or empty input: $file" >&2; exit 1; }
done
mkdir -p "$OUTDIR"

align_haplotype() {
    local hap_name="$1" reference="$2"
    local prefix="${OUTDIR}/${hap_name}"
    local primary_count mapq20_count suffix need_index=0
    echo "Processing ${hap_name}: ${reference}"

    # 3. Build reference indexes beside each FASTA when needed.
    # Reference folders must be writable for this step.
    for suffix in amb ann bwt pac sa; do
        [[ -s "${reference}.${suffix}" ]] || need_index=1
    done
    if (( need_index )); then bwa index "$reference"; fi
    if [[ ! -s "${reference}.fai" ]]; then samtools faidx "$reference"; fi

    # 4. Retain all reported matches in a coordinate-sorted BAM.
    # -a allows secondary matches; primary-only filtering occurs below.
    # This workflow does not require a marker to have exactly one genomic match.
    bwa mem -t "$BWA_THREADS" -a "$reference" "$MARKERS" |
        samtools sort -@ "$SORT_THREADS" -o "${prefix}.all.sorted.bam"
    samtools index "${prefix}.all.sorted.bam"
    samtools flagstat -@ "$THREADS" "${prefix}.all.sorted.bam" > "${prefix}.flagstat.txt"

    # 5. Exclude unmapped (4), secondary (256), and supplementary (2048)
    # records: 4 + 256 + 2048 = 2308. Keep a separate MAPQ >=20 subset.
    samtools view -@ "$THREADS" -b -F 2308 "${prefix}.all.sorted.bam" > "${prefix}.primary_mapped.bam"
    samtools index "${prefix}.primary_mapped.bam"
    samtools view -@ "$THREADS" -b -q 20 "${prefix}.primary_mapped.bam" > "${prefix}.primary_MAPQ20.bam"
    samtools index "${prefix}.primary_MAPQ20.bam"

    # 6. Export primary records for CIGAR-based SNP projection in step 04.
    # SAM POS is one-based. It is the alignment start, not the SNP coordinate.
    # Test flag bit 16 arithmetically for portability across awk versions.
    samtools view "${prefix}.primary_mapped.bam" |
        awk -v hap="$hap_name" 'BEGIN {OFS="\t"}
        {strand = int($2/16)%2 ? "-" : "+";
         print $1,hap,$3,$4,$5,strand,$6,length($10)}' > "${prefix}.primary_positions.noheader.tsv"
    {
        printf 'Marker.Name\tHaplotype\tRegenSY.Sequence\tAlignment.Start\tMAPQ\tStrand\tCIGAR\tQuery.Length\n'
        cat "${prefix}.primary_positions.noheader.tsv"
    } > "${prefix}.primary_positions.tsv"
    primary_count=$(samtools view -c "${prefix}.primary_mapped.bam")
    mapq20_count=$(samtools view -c "${prefix}.primary_MAPQ20.bam")
    echo "${hap_name}: primary mapped = ${primary_count}; MAPQ >=20 = ${mapq20_count}"
}

# 7. Analyze assemblies sequentially, retaining independent output files.
align_haplotype Hap1 "$HAP1"
align_haplotype Hap2 "$HAP2"
align_haplotype Hap3 "$HAP3"
align_haplotype Hap4 "$HAP4"

# 8. Combine primary positions with a single header for step 04.
{
    printf 'Marker.Name\tHaplotype\tRegenSY.Sequence\tAlignment.Start\tMAPQ\tStrand\tCIGAR\tQuery.Length\n'
    for hap in Hap1 Hap2 Hap3 Hap4; do
        cat "${OUTDIR}/${hap}.primary_positions.noheader.tsv"
    done
} > "${OUTDIR}/all_haplotypes.primary_positions.tsv"

# 9. Summarize primary and high-quality alignment counts by haplotype.
{
    printf 'Haplotype\tPrimary.Mapped\tPrimary.MAPQ20\n'
    for hap in Hap1 Hap2 Hap3 Hap4; do
        primary_count=$(samtools view -c "${OUTDIR}/${hap}.primary_mapped.bam")
        mapq20_count=$(samtools view -c "${OUTDIR}/${hap}.primary_MAPQ20.bam")
        printf '%s\t%s\t%s\n' "$hap" "$primary_count" "$mapq20_count"
    done
} > "${OUTDIR}/alignment_summary.tsv"
cat "${OUTDIR}/alignment_summary.tsv"
