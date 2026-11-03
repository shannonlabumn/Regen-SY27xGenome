#!/bin/bash -l
#SBATCH --job-name=cafe5
#SBATCH --output=cafe5_%j.out
#SBATCH --time=50:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G


# Treat unset variables as errors and propagate pipeline failures.
# Model failures are handled explicitly so the cutoff scan can continue.
# No global errexit is enabled; inspect model statuses and logs after a run.
set -uo pipefail



# CAFE5 gene-family analyses for alfalfa haplotypes 1-4.
# Requirements: Bash >=4, Python 3, CAFE5, and standard Unix utilities.
# Activate your own software environment before running this script.
# Usage: bash cafe5.sh  OR  sbatch cafe5.sh
# Optional: INPUT_ROOT=/path/to/inputs OUTPUT_ROOT=/path/to/results bash cafe5.sh
#
# Default layout (relative to this script):
# data/orthofinder/haplotype1/Orthogroups/Orthogroups.GeneCount.tsv
# data/orthofinder/haplotype1/Species_Tree/SpeciesTree_dated_Mya.nwk
# Repeat this layout for haplotype2, haplotype3, and haplotype4.
# Each dataset must contain one haplotype per haplotype-resolved genome,
# plus the other species included in the analysis. This script does not
# select haplotypes from a combined OrthoFinder dataset.
# Dated trees must be prepared separately; tip names must match count headers.
# Outputs include prepared counts, model logs, parameter estimates, cutoff
# scans, incremental candidate tests, and per-haplotype/master summaries.
# The large-family workflow scans the cutoffs below, chooses the largest
# successful cutoff, and tests remaining families incrementally. It does
# not impose a fixed final difference cutoff (e.g., 300).

# STEP 1: Configure portable paths and execution settings.
# Default paths are anchored to this script, independent of the launch folder.
# Environment variables override defaults without editing the public script.
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
INPUT_ROOT="${INPUT_ROOT:-${SCRIPT_DIR}/data/orthofinder}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${SCRIPT_DIR}/results/cafe5}"
# Use the requested thread count, then the SLURM allocation, then 32.
# When submitting with SLURM, keep THREADS within the allocated CPU count.
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-32}}"
CAFE5_BIN="${CAFE5_BIN:-cafe5}"
PYTHON_BIN="${PYTHON_BIN:-python3}"
# Classify by the largest count across all species: <100 is small; >=100 is large.
LARGE_FAMILY_SIZE=100
# FORCE=1 archives existing haplotype folders before rerunning them.
FORCE="${FORCE:-0}"
HAPLOTYPES=(1 2 3 4)

# STEP 2: Define candidate limits on max(counts) - min(counts).
# These are absolute copy-number differences, not percentages or P-value cutoffs.
# Keep the limits in ascending order: each successful run replaces the best limit.
# Difference thresholds for large-family scans.
CUTOFFS=(
    100
    110
    120
    125
    130
    140
    150
    175
    200
    225
    250
    275
    300
    350
    400
    500
    600
    700
    800
)


# ============================================================
# 3. ORTHOFINDER DIRECTORIES
# ============================================================

# Map each independent haplotype analysis to its count table and dated tree.
# Matching haplotype labels across genomes identify datasets, not additional
# independent lineages within the same CAFE5 tree.
declare -A ORTHOFINDER_DIR TREE_FILE
GENECOUNT_RELATIVE="Orthogroups/Orthogroups.GeneCount.tsv"
for HAP in "${HAPLOTYPES[@]}"; do
    ORTHOFINDER_DIR[$HAP]="${INPUT_ROOT}/haplotype${HAP}"
    TREE_FILE[$HAP]="${ORTHOFINDER_DIR[$HAP]}/Species_Tree/SpeciesTree_dated_Mya.nwk"
done
# Edit the mapping above if your input folders or dated-tree filenames differ.

# ============================================================
# 5. HELPER FUNCTIONS
# ============================================================

# Return a readable time for progress messages.
timestamp()
{
    date "+%Y-%m-%d %H:%M:%S"
}


# Write a timestamped progress message to the job/terminal log.
log_message()
{
    printf "[%s] %s\n" "$(timestamp)" "$*"
}


# Report a fatal setup or analysis error and stop the entire script.
die()
{
    printf "ERROR: %s\n" "$*" >&2
    exit 1
}


# Count data rows in a generated TSV, excluding its single header row.
family_count()
{
    local input_file="$1"
    local number

    number=$(wc -l < "$input_file")

    if (( number > 0 )); then
        printf "%s\n" "$((number - 1))"
    else
        printf "0\n"
    fi
}


# Check for expected progress markers and absence of initialization failure.
# This is a log-based heuristic, not an independent convergence diagnostic.
# The wrappers below also require CAFE5 to return exit code zero.
cafe_succeeded()
{
    local logfile="$1"

    grep -q "Computing pvalues" "$logfile" &&
        grep -q "Starting reconstruction processes" "$logfile" &&
        ! grep -q "Failed to initialize" "$logfile"
}


# Extract the last reported number of families retained by root filtering.
# Return NA when the expected message is absent from the CAFE5 log.
root_retained()
{
    local logfile="$1"
    local retained

    retained=$(
        grep "Filtering families not present at the root" "$logfile" |
            tail -1 |
            awk '{print $NF}'
    )

    if [[ -z "$retained" ]]; then
        retained="NA"
    fi

    printf "%s\n" "$retained"
}


# Estimate the evolutionary rate and Poisson root-distribution parameter
# from small families. Save both the exact command and CAFE5 output to a log.
run_small_model()
{
    local input_file="$1"
    local tree_file="$2"
    local output_prefix="$3"
    local logfile="$4"

    {
        printf "Command: %q" "$CAFE5_BIN"
        printf " %q" \
            -i "$input_file" \
            -t "$tree_file" \
            -c "$THREADS" \
            -p \
            -o "$output_prefix"
        printf "\n\n"

        "$CAFE5_BIN" \
            -i "$input_file" \
            -t "$tree_file" \
            -c "$THREADS" \
            -p \
            -o "$output_prefix"

    } > "$logfile" 2>&1

    local cafe_exit=$?
    (( cafe_exit == 0 )) && cafe_succeeded "$logfile"
}


# Analyze a large-family subset using the parameters estimated from small
# families in this haplotype analysis; do not re-estimate them for every subset.
run_fixed_model()
{
    local input_file="$1"
    local tree_file="$2"
    local evolutionary_lambda="$3"
    local poisson_lambda="$4"
    local output_prefix="$5"
    local logfile="$6"

    {
        printf "Command: %q" "$CAFE5_BIN"
        printf " %q" \
            -i "$input_file" \
            -t "$tree_file" \
            -c "$THREADS" \
            -l "$evolutionary_lambda" \
            "--poisson=${poisson_lambda}" \
            -o "$output_prefix"
        printf "\n\n"

        "$CAFE5_BIN" \
            -i "$input_file" \
            -t "$tree_file" \
            -c "$THREADS" \
            -l "$evolutionary_lambda" \
            "--poisson=${poisson_lambda}" \
            -o "$output_prefix"

    } > "$logfile" 2>&1

    local cafe_exit=$?
    (( cafe_exit == 0 )) && cafe_succeeded "$logfile"
}


# ============================================================
# 6. CHECK INPUTS
# ============================================================

# STEP 3: Check executable availability, numeric settings, and all input files
# before starting the four sequential analyses. Input files are never modified.
command -v "$CAFE5_BIN" >/dev/null 2>&1 || die "CAFE5 executable not found: $CAFE5_BIN"
command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python executable not found: $PYTHON_BIN"
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || die "THREADS must be a positive integer"
[[ "$FORCE" == 0 || "$FORCE" == 1 ]] || die "FORCE must be 0 or 1"
mkdir -p "$OUTPUT_ROOT" || die "Cannot create output directory: $OUTPUT_ROOT"

for HAP in "${HAPLOTYPES[@]}"
do
    RAW_COUNTS="${ORTHOFINDER_DIR[$HAP]}/${GENECOUNT_RELATIVE}"

    [[ -f "$RAW_COUNTS" ]] ||
        die "Haplotype ${HAP} count file not found: ${RAW_COUNTS}"

    [[ -f "${TREE_FILE[$HAP]}" ]] ||
        die "Haplotype ${HAP} dated tree not found: ${TREE_FILE[$HAP]}"
done


# Refuse existing haplotype outputs before resetting the master summary.
if [[ "$FORCE" == 0 ]]; then
    for HAP in "${HAPLOTYPES[@]}"; do
        [[ ! -e "${OUTPUT_ROOT}/haplotype${HAP}" ]] ||
            die "Haplotype ${HAP} output already exists. Use a new OUTPUT_ROOT or FORCE=1."
    done
fi

# ============================================================
# 7. MASTER SUMMARY
# ============================================================

# STEP 4: Initialize one row per haplotype with model sizes and run statuses.
# This master summary is recreated for each invocation.
MASTER_SUMMARY="${OUTPUT_ROOT}/CAFE5_haplotypes_1_2_3_4_summary.tsv"

printf \
"Haplotype\tSmallInput\tSmallRootRetained\tEvolutionaryLambda\tPoissonLambda\tLargeInput\tLargestSuccessfulCutoff\tFinalLargeInput\tFinalLargeRootRetained\tSmallStatus\tLargeStatus\n" \
> "$MASTER_SUMMARY"


# ============================================================
# 8. RUN EACH HAPLOTYPE
# ============================================================

for HAP in "${HAPLOTYPES[@]}"
do
    log_message "===================================================="
    log_message "Starting CAFE5 analysis for haplotype ${HAP}"
    log_message "===================================================="

    # STEP 5: Establish this haplotype's paths and separate output directories.
    # Each haplotype gets its own parameter estimates and family selections.
    RESULTS_DIR="${ORTHOFINDER_DIR[$HAP]}"
    RAW_COUNTS="${RESULTS_DIR}/${GENECOUNT_RELATIVE}"
    DATED_TREE="${TREE_FILE[$HAP]}"

    HAP_DIR="${OUTPUT_ROOT}/haplotype${HAP}"

    if [[ -e "$HAP_DIR" ]]; then
        if [[ "$FORCE" -eq 1 ]]; then
            ARCHIVE_DIR="${HAP_DIR}_previous_$(date +%Y%m%d_%H%M%S)"

            log_message "Archiving existing directory:"
            log_message "  ${HAP_DIR}"
            log_message "  -> ${ARCHIVE_DIR}"

            mv "$HAP_DIR" "$ARCHIVE_DIR"
        else
            die "${HAP_DIR} already exists. Set FORCE=1 to archive and rerun."
        fi
    fi

    mkdir -p \
        "${HAP_DIR}/input" \
        "${HAP_DIR}/small" \
        "${HAP_DIR}/large/scans" \
        "${HAP_DIR}/large/candidates" \
        "${HAP_DIR}/large/final"

    # Keep the full converted table, size-based subsets, and family statistics
    # for traceability of every subsequent filter.
    PREPARED_COUNTS="${HAP_DIR}/input/Orthogroups.GeneCount.CAFE.tsv"
    SMALL_COUNTS="${HAP_DIR}/input/CAFE_small_lt${LARGE_FAMILY_SIZE}.tsv"
    LARGE_COUNTS="${HAP_DIR}/input/CAFE_large_ge${LARGE_FAMILY_SIZE}.tsv"
    FAMILY_STATS="${HAP_DIR}/input/family_size_statistics.tsv"

    # --------------------------------------------------------
    # 8.1 Prepare CAFE5 input tables
    # --------------------------------------------------------

    # STEP 6: Convert OrthoFinder counts to Desc, Orthogroup, species columns.
    # Remove Total, add Desc=null, and calculate minimum, maximum, and range.
    # Split all families by maximum count, retaining IDs and species order.
    log_message "Preparing count tables for haplotype ${HAP}"

    "$PYTHON_BIN" - \
        "$RAW_COUNTS" \
        "$PREPARED_COUNTS" \
        "$SMALL_COUNTS" \
        "$LARGE_COUNTS" \
        "$FAMILY_STATS" \
        "$LARGE_FAMILY_SIZE" <<'PYTHON'
import csv
import sys
from pathlib import Path

(
    raw_file,
    prepared_file,
    small_file,
    large_file,
    statistics_file,
    large_threshold,
) = sys.argv[1:]

large_threshold = int(large_threshold)

raw_file = Path(raw_file)
prepared_file = Path(prepared_file)
small_file = Path(small_file)
large_file = Path(large_file)
statistics_file = Path(statistics_file)

with raw_file.open(
    "r",
    encoding="utf-8-sig",
    newline="",
) as handle:

    reader = csv.reader(handle, delimiter="\t")

    try:
        original_header = next(reader)
    except StopIteration:
        raise SystemExit(
            f"ERROR: Empty input file: {raw_file}"
        )

    original_header = [
        value.strip().rstrip("\r")
        for value in original_header
    ]

    if "Orthogroup" not in original_header:
        raise SystemExit(
            f"ERROR: Orthogroup column not found in {raw_file}"
        )

    # Remove only the OrthoFinder Total column.
    keep_indices = [
        index
        for index, column in enumerate(original_header)
        if column != "Total"
    ]

    retained_header = [
        original_header[index]
        for index in keep_indices
    ]

    if retained_header[0] != "Orthogroup":
        raise SystemExit(
            "ERROR: Orthogroup is not the first retained column."
        )

    # CAFE5 input begins with Desc followed by Orthogroup.
    cafe_header = ["Desc"] + retained_header

    all_rows = []
    small_rows = []
    large_rows = []
    statistics_rows = []

    for line_number, original_row in enumerate(
        reader,
        start=2,
    ):
        if not original_row:
            continue

        original_row = [
            value.strip().rstrip("\r")
            for value in original_row
        ]

        if len(original_row) < len(original_header):
            original_row.extend(
                [""] * (
                    len(original_header) -
                    len(original_row)
                )
            )

        retained_row = [
            original_row[index]
            for index in keep_indices
        ]

        orthogroup = retained_row[0]

        if not orthogroup:
            continue

        counts = []

        for species, value in zip(
            retained_header[1:],
            retained_row[1:],
        ):
            try:
                count = int(float(value))
            except ValueError:
                raise SystemExit(
                    f"ERROR: Non-numeric count at line "
                    f"{line_number}, species {species}: {value!r}"
                )

            if count < 0:
                raise SystemExit(
                    f"ERROR: Negative family count for "
                    f"{orthogroup}, species {species}"
                )

            counts.append(count)

        minimum = min(counts)
        maximum = max(counts)
        difference = maximum - minimum

        # Add required Desc column with null values.
        cafe_row = [
            "null",
            orthogroup,
            *counts,
        ]

        all_rows.append(cafe_row)

        if maximum < large_threshold:
            dataset = "small"
            small_rows.append(cafe_row)
        else:
            dataset = "large"
            large_rows.append(cafe_row)

        statistics_rows.append(
            [
                orthogroup,
                minimum,
                maximum,
                difference,
                dataset,
            ]
        )


def write_tsv(filename, header, rows):
    with filename.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as output:

        writer = csv.writer(
            output,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(header)
        writer.writerows(rows)


write_tsv(
    prepared_file,
    cafe_header,
    all_rows,
)

write_tsv(
    small_file,
    cafe_header,
    small_rows,
)

write_tsv(
    large_file,
    cafe_header,
    large_rows,
)

statistics_rows.sort(
    key=lambda row: (
        -row[3],
        row[0],
    )
)

write_tsv(
    statistics_file,
    [
        "Orthogroup",
        "Minimum",
        "Maximum",
        "Difference",
        "Dataset",
    ],
    statistics_rows,
)

print(f"All families: {len(all_rows)}")
print(f"Small families: {len(small_rows)}")
print(f"Large families: {len(large_rows)}")
PYTHON

    if [[ $? -ne 0 ]]; then
        die "Count-table preparation failed for haplotype ${HAP}"
    fi

    SMALL_INPUT_N=$(family_count "$SMALL_COUNTS")
    LARGE_INPUT_N=$(family_count "$LARGE_COUNTS")

    log_message "Small-family input: ${SMALL_INPUT_N}"
    log_message "Large-family input: ${LARGE_INPUT_N}"

    # --------------------------------------------------------
    # 8.2 Validate tree tips against count-table species
    # --------------------------------------------------------

    # STEP 7: Require the same unique species labels in the tree and counts.
    # The parser assumes unquoted Newick tip names with branch lengths.
    # This checks names only; it does not establish dating or ultrametricity.
    log_message "Validating tree and count-table species"

    "$PYTHON_BIN" - \
        "$PREPARED_COUNTS" \
        "$DATED_TREE" <<'PYTHON'
import csv
import re
import sys
from pathlib import Path

count_file, tree_file = sys.argv[1:]

with open(
    count_file,
    "r",
    encoding="utf-8-sig",
    newline="",
) as handle:
    reader = csv.reader(handle, delimiter="\t")
    header = next(reader)

count_species = [
    name.strip().rstrip("\r")
    for name in header[2:]
]

tree_text = Path(tree_file).read_text(
    encoding="utf-8"
).strip()

# Terminal tip names occur after "(" or "," and before ":".
tree_tips = re.findall(
    r"(?<=[(,])([^():,;]+)(?=:)",
    tree_text,
)

tree_tips = [
    name.strip()
    for name in tree_tips
]

count_set = set(count_species)
tree_set = set(tree_tips)

missing_from_tree = sorted(
    count_set - tree_set
)

missing_from_counts = sorted(
    tree_set - count_set
)

print("Count-table species:")
for name in count_species:
    print(f"  {name}")

print("\nTree tips:")
for name in tree_tips:
    print(f"  {name}")

if missing_from_tree:
    print("\nIn count table but missing from tree:")
    for name in missing_from_tree:
        print(f"  {name}")

if missing_from_counts:
    print("\nIn tree but missing from count table:")
    for name in missing_from_counts:
        print(f"  {name}")

if missing_from_tree or missing_from_counts:
    raise SystemExit(
        "ERROR: Tree tips and count-table species do not match."
    )

if len(count_species) != len(count_set):
    raise SystemExit(
        "ERROR: Duplicate species names in count table."
    )

if len(tree_tips) != len(tree_set):
    raise SystemExit(
        "ERROR: Duplicate terminal names in tree."
    )

print("\nTree/count validation passed.")
PYTHON

    if [[ $? -ne 0 ]]; then
        die "Tree/count validation failed for haplotype ${HAP}"
    fi

    # --------------------------------------------------------
    # 8.3 Fit small-family Base Poisson model
    # --------------------------------------------------------

    # STEP 8: Fit the small-family model and stop if its fit fails.
    # Later large-family runs depend on these fitted parameters.
    SMALL_PREFIX="${HAP_DIR}/small/cafe5_hap${HAP}_small_poisson"
    SMALL_LOG="${HAP_DIR}/small/cafe5_hap${HAP}_small_poisson.out"

    log_message "Running small-family Base Poisson model"

    if run_small_model \
        "$SMALL_COUNTS" \
        "$DATED_TREE" \
        "$SMALL_PREFIX" \
        "$SMALL_LOG"
    then
        SMALL_STATUS="SUCCESS"
    else
        SMALL_STATUS="FAILED"
    fi

    if [[ "$SMALL_STATUS" != "SUCCESS" ]]; then
        printf \
"Haplotype%s\t%s\tNA\tNA\tNA\t%s\tNA\tNA\tNA\tFAILED\tNOT_RUN\n" \
            "$HAP" \
            "$SMALL_INPUT_N" \
            "$LARGE_INPUT_N" \
            >> "$MASTER_SUMMARY"

        die "Small-family model failed for haplotype ${HAP}: ${SMALL_LOG}"
    fi

    SMALL_ROOT_RETAINED=$(root_retained "$SMALL_LOG")

    # Read the final parameter messages from the small-family log.
    # Extraction relies on the CAFE5 log format used by the supplied workflow.
    EVOLUTIONARY_LAMBDA=$(
        grep "Best match is:" "$SMALL_LOG" |
            tail -1 |
            awk '{print $4}'
    )

    POISSON_LAMBDA=$(
        grep "Poisson lambda:" "$SMALL_LOG" |
            tail -1 |
            awk '{print $3}'
    )

    [[ -n "$EVOLUTIONARY_LAMBDA" ]] ||
        die "Could not extract evolutionary lambda for haplotype ${HAP}"

    [[ -n "$POISSON_LAMBDA" ]] ||
        die "Could not extract root Poisson lambda for haplotype ${HAP}"

    log_message "Evolutionary lambda: ${EVOLUTIONARY_LAMBDA}"
    log_message "Root Poisson lambda: ${POISSON_LAMBDA}"
    log_message "Root-retained small families: ${SMALL_ROOT_RETAINED}"

    {
        printf "Parameter\tValue\n"
        printf "EvolutionaryLambda\t%s\n" "$EVOLUTIONARY_LAMBDA"
        printf "RootPoissonLambda\t%s\n" "$POISSON_LAMBDA"
    } > "${HAP_DIR}/small/model_parameters.tsv"

    # --------------------------------------------------------
    # 8.4 Create large-family statistics
    # --------------------------------------------------------

    # STEP 9: Extract only large-family statistics for cutoff and candidate lists.
    LARGE_STATS="${HAP_DIR}/input/large_family_statistics.tsv"

    awk -F'\t' '
        BEGIN {
            OFS="\t"
        }

        NR==1 {
            print "Orthogroup", "Minimum", "Maximum", "Difference"
            next
        }

        $5=="large" {
            print $1, $2, $3, $4
        }
    ' "$FAMILY_STATS" > "$LARGE_STATS"

    # --------------------------------------------------------
    # 8.5 Generate large-family cutoff inputs
    # --------------------------------------------------------

    # STEP 10/11: Process all configured copy-number difference cutoffs.
    # Input generation retains large families with range <= the given cutoff;
    # the later model scan records success/failure and root-retained counts.
    for CUTOFF in "${CUTOFFS[@]}"
    do
        CUTOFF_INPUT="${HAP_DIR}/large/scans/CAFE_large_diff_le${CUTOFF}.tsv"

        "$PYTHON_BIN" - \
            "$LARGE_COUNTS" \
            "$CUTOFF_INPUT" \
            "$CUTOFF" <<'PYTHON'
import csv
import sys

input_file, output_file, cutoff = sys.argv[1:]
cutoff = int(cutoff)

with open(
    input_file,
    "r",
    encoding="utf-8-sig",
    newline="",
) as source:

    reader = csv.reader(source, delimiter="\t")
    header = next(reader)

    retained = []

    for row in reader:
        if not row:
            continue

        counts = [
            int(float(value))
            for value in row[2:]
        ]

        difference = max(counts) - min(counts)

        if difference <= cutoff:
            retained.append(row)

with open(
    output_file,
    "w",
    encoding="utf-8",
    newline="",
) as destination:

    writer = csv.writer(
        destination,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(header)
    writer.writerows(retained)
PYTHON
    done

    # --------------------------------------------------------
    # 8.6 Run cutoff scan
    # --------------------------------------------------------

    # Run the fixed-parameter model separately on every nonempty cutoff input.
    # BEST_CUTOFF tracks the last successful cutoff in the ascending list.
    SCAN_SUMMARY="${HAP_DIR}/large/large_cutoff_scan_summary.tsv"

    printf \
"Cutoff\tInputFamilies\tRootRetained\tStatus\n" \
        > "$SCAN_SUMMARY"

    BEST_CUTOFF=""
    BEST_CUTOFF_INPUT=""

    log_message "Running large-family cutoff scan"

    # STEP 10/11: Process all configured copy-number difference cutoffs.
    # Input generation retains large families with range <= the given cutoff;
    # the later model scan records success/failure and root-retained counts.
    for CUTOFF in "${CUTOFFS[@]}"
    do
        CUTOFF_INPUT="${HAP_DIR}/large/scans/CAFE_large_diff_le${CUTOFF}.tsv"
        CUTOFF_N=$(family_count "$CUTOFF_INPUT")

        if (( CUTOFF_N == 0 )); then
            printf \
"%s\t0\tNA\tEMPTY\n" \
                "$CUTOFF" \
                >> "$SCAN_SUMMARY"

            continue
        fi

        CUTOFF_PREFIX="${HAP_DIR}/large/scans/cafe5_hap${HAP}_large_le${CUTOFF}"
        CUTOFF_LOG="${HAP_DIR}/large/scans/cafe5_hap${HAP}_large_le${CUTOFF}.out"

        log_message \
            "Haplotype ${HAP}: cutoff ${CUTOFF}, families ${CUTOFF_N}"

        if run_fixed_model \
            "$CUTOFF_INPUT" \
            "$DATED_TREE" \
            "$EVOLUTIONARY_LAMBDA" \
            "$POISSON_LAMBDA" \
            "$CUTOFF_PREFIX" \
            "$CUTOFF_LOG"
        then
            CUTOFF_STATUS="SUCCESS"

            BEST_CUTOFF="$CUTOFF"
            BEST_CUTOFF_INPUT="$CUTOFF_INPUT"
        else
            CUTOFF_STATUS="FAILED"
        fi

        CUTOFF_ROOT=$(root_retained "$CUTOFF_LOG")

        printf \
"%s\t%s\t%s\t%s\n" \
            "$CUTOFF" \
            "$CUTOFF_N" \
            "$CUTOFF_ROOT" \
            "$CUTOFF_STATUS" \
            >> "$SCAN_SUMMARY"
    done

    # --------------------------------------------------------
    # 8.7 No successful large-family subset
    # --------------------------------------------------------

    # STEP 12: If no cutoff succeeds, record the outcome and move to the next
    # haplotype. Its small-family results remain available.
    if [[ -z "$BEST_CUTOFF" ]]; then
        log_message "No successful large-family subset for haplotype ${HAP}"

        printf \
"Haplotype%s\t%s\t%s\t%s\t%s\t%s\tNA\t0\t0\t%s\tNO_SUCCESSFUL_SUBSET\n" \
            "$HAP" \
            "$SMALL_INPUT_N" \
            "$SMALL_ROOT_RETAINED" \
            "$EVOLUTIONARY_LAMBDA" \
            "$POISSON_LAMBDA" \
            "$LARGE_INPUT_N" \
            "$SMALL_STATUS" \
            >> "$MASTER_SUMMARY"

        continue
    fi

    log_message "Largest successful cutoff: ${BEST_CUTOFF}"

    # --------------------------------------------------------
    # 8.8 Test remaining large families incrementally
    # --------------------------------------------------------

    # STEP 13: Seed the accepted set with the largest successful cutoff subset.
    # Sort excluded candidates by increasing range, then family ID.
    # Acceptance is sequential and can depend on the order of candidate testing.
    CURRENT_ACCEPTED="${HAP_DIR}/large/candidates/current_accepted.tsv"
    cp "$BEST_CUTOFF_INPUT" "$CURRENT_ACCEPTED"

    BASE_SUCCESS_LOG="${HAP_DIR}/large/scans/cafe5_hap${HAP}_large_le${BEST_CUTOFF}.out"
    CURRENT_ROOT_RETAINED=$(root_retained "$BASE_SUCCESS_LOG")

    CANDIDATE_LIST="${HAP_DIR}/large/candidates/candidates_above_cutoff.tsv"

    awk -F'\t' \
        -v cutoff="$BEST_CUTOFF" '
        BEGIN {
            OFS="\t"
        }

        NR==1 {
            next
        }

        $4 > cutoff {
            print $1, $4, $2, $3
        }
    ' "$LARGE_STATS" |
        sort -t$'\t' -k2,2n -k1,1 \
        > "$CANDIDATE_LIST"

    CANDIDATE_SUMMARY="${HAP_DIR}/large/candidate_testing_summary.tsv"

    printf \
"Family\tDifference\tBeforeRootRetained\tAfterRootRetained\tStatus\n" \
        > "$CANDIDATE_SUMMARY"

    CANDIDATE_NUMBER=0

    # Add one candidate to the current set and rerun the fixed-parameter model.
    # Keep it only when the trial succeeds; otherwise retain the previous set.
    while IFS=$'\t' read -r FAMILY DIFFERENCE MINIMUM MAXIMUM
    do
        [[ -n "$FAMILY" ]] || continue

        CANDIDATE_NUMBER=$((CANDIDATE_NUMBER + 1))

        TRIAL_INPUT="${HAP_DIR}/large/candidates/trial_${CANDIDATE_NUMBER}_${FAMILY}.tsv"
        TRIAL_PREFIX="${HAP_DIR}/large/candidates/cafe5_trial_${CANDIDATE_NUMBER}_${FAMILY}"
        TRIAL_LOG="${HAP_DIR}/large/candidates/cafe5_trial_${CANDIDATE_NUMBER}_${FAMILY}.out"

        "$PYTHON_BIN" - \
            "$CURRENT_ACCEPTED" \
            "$LARGE_COUNTS" \
            "$FAMILY" \
            "$TRIAL_INPUT" <<'PYTHON'
import csv
import sys

current_file, all_large_file, family, output_file = sys.argv[1:]

with open(
    current_file,
    "r",
    encoding="utf-8-sig",
    newline="",
) as handle:
    current_data = list(
        csv.reader(handle, delimiter="\t")
    )

if not current_data:
    raise SystemExit(
        "ERROR: Current accepted file is empty."
    )

header = current_data[0]
rows = current_data[1:]

existing_families = {
    row[1]
    for row in rows
    if len(row) > 1
}

if family not in existing_families:
    family_found = False

    with open(
        all_large_file,
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as handle:

        reader = csv.reader(handle, delimiter="\t")
        source_header = next(reader)

        if source_header != header:
            raise SystemExit(
                "ERROR: Header mismatch in large-family files."
            )

        for row in reader:
            if len(row) > 1 and row[1] == family:
                rows.append(row)
                family_found = True
                break

    if not family_found:
        raise SystemExit(
            f"ERROR: Candidate family not found: {family}"
        )

with open(
    output_file,
    "w",
    encoding="utf-8",
    newline="",
) as handle:

    writer = csv.writer(
        handle,
        delimiter="\t",
        lineterminator="\n",
    )

    writer.writerow(header)
    writer.writerows(rows)
PYTHON

        # Compare aggregate root-retained counts before and after the trial.
        # An unchanged count is labeled ROOT_FILTERED_ACCEPTED by this workflow;
        # this count comparison does not directly verify individual family status.
        BEFORE_ROOT="$CURRENT_ROOT_RETAINED"

        log_message \
            "Testing ${FAMILY}, difference ${DIFFERENCE}"

        if run_fixed_model \
            "$TRIAL_INPUT" \
            "$DATED_TREE" \
            "$EVOLUTIONARY_LAMBDA" \
            "$POISSON_LAMBDA" \
            "$TRIAL_PREFIX" \
            "$TRIAL_LOG"
        then
            AFTER_ROOT=$(root_retained "$TRIAL_LOG")

            if [[ "$AFTER_ROOT" == "$BEFORE_ROOT" ]]; then
                TRIAL_STATUS="ROOT_FILTERED_ACCEPTED"
            else
                TRIAL_STATUS="TESTABLE_ACCEPTED"
            fi

            cp "$TRIAL_INPUT" "$CURRENT_ACCEPTED"
            CURRENT_ROOT_RETAINED="$AFTER_ROOT"
        else
            AFTER_ROOT=$(root_retained "$TRIAL_LOG")
            TRIAL_STATUS="NUMERICAL_FAILURE_REJECTED"
        fi

        printf \
"%s\t%s\t%s\t%s\t%s\n" \
            "$FAMILY" \
            "$DIFFERENCE" \
            "$BEFORE_ROOT" \
            "$AFTER_ROOT" \
            "$TRIAL_STATUS" \
            >> "$CANDIDATE_SUMMARY"

    done < "$CANDIDATE_LIST"

    # --------------------------------------------------------
    # 8.9 Final large-family run
    # --------------------------------------------------------

    # STEP 14: Save the accepted large-family table and rerun the final model.
    # Report its status explicitly; do not infer success from file creation alone.
    FINAL_LARGE_INPUT="${HAP_DIR}/input/CAFE_large_final.tsv"
    cp "$CURRENT_ACCEPTED" "$FINAL_LARGE_INPUT"

    FINAL_LARGE_N=$(family_count "$FINAL_LARGE_INPUT")

    FINAL_PREFIX="${HAP_DIR}/large/final/cafe5_hap${HAP}_large_final"
    FINAL_LOG="${HAP_DIR}/large/final/cafe5_hap${HAP}_large_final.out"

    log_message \
        "Running final large-family model with ${FINAL_LARGE_N} families"

    if run_fixed_model \
        "$FINAL_LARGE_INPUT" \
        "$DATED_TREE" \
        "$EVOLUTIONARY_LAMBDA" \
        "$POISSON_LAMBDA" \
        "$FINAL_PREFIX" \
        "$FINAL_LOG"
    then
        LARGE_STATUS="SUCCESS"
    else
        LARGE_STATUS="FAILED"
    fi

    FINAL_LARGE_ROOT=$(root_retained "$FINAL_LOG")

    # --------------------------------------------------------
    # 8.10 Haplotype summary
    # --------------------------------------------------------

    # STEP 15: Write a detailed haplotype summary and append the master row.
    # Input counts precede root filtering; root-retained counts are separate.
    HAP_SUMMARY="${HAP_DIR}/CAFE5_haplotype${HAP}_summary.tsv"

    {
        printf "Metric\tValue\n"
        printf "Haplotype\t%s\n" "$HAP"
        printf "OrthoFinderDirectory\t%s\n" "$RESULTS_DIR"
        printf "CountFile\t%s\n" "$RAW_COUNTS"
        printf "DatedTree\t%s\n" "$DATED_TREE"
        printf "SmallInputFamilies\t%s\n" "$SMALL_INPUT_N"
        printf "SmallRootRetained\t%s\n" "$SMALL_ROOT_RETAINED"
        printf "EvolutionaryLambda\t%s\n" "$EVOLUTIONARY_LAMBDA"
        printf "RootPoissonLambda\t%s\n" "$POISSON_LAMBDA"
        printf "LargeInputFamilies\t%s\n" "$LARGE_INPUT_N"
        printf "LargestSuccessfulCutoff\t%s\n" "$BEST_CUTOFF"
        printf "FinalLargeInputFamilies\t%s\n" "$FINAL_LARGE_N"
        printf "FinalLargeRootRetained\t%s\n" "$FINAL_LARGE_ROOT"
        printf "SmallStatus\t%s\n" "$SMALL_STATUS"
        printf "LargeStatus\t%s\n" "$LARGE_STATUS"
    } > "$HAP_SUMMARY"

    printf \
"Haplotype%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$HAP" \
        "$SMALL_INPUT_N" \
        "$SMALL_ROOT_RETAINED" \
        "$EVOLUTIONARY_LAMBDA" \
        "$POISSON_LAMBDA" \
        "$LARGE_INPUT_N" \
        "$BEST_CUTOFF" \
        "$FINAL_LARGE_N" \
        "$FINAL_LARGE_ROOT" \
        "$SMALL_STATUS" \
        "$LARGE_STATUS" \
        >> "$MASTER_SUMMARY"

    log_message "Completed haplotype ${HAP}"
done


# ============================================================
# 9. FINAL REPORT
# ============================================================

log_message "===================================================="
# STEP 16: Display the summary for review. Completion of the workflow does
# not imply every model succeeded; inspect SmallStatus and LargeStatus.
log_message "CAFE5 analyses completed"
log_message "Master summary: ${MASTER_SUMMARY}"
log_message "===================================================="

if command -v column >/dev/null 2>&1; then
    column -t -s $'\t' "$MASTER_SUMMARY"
else
    cat "$MASTER_SUMMARY"
fi
