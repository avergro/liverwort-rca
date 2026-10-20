#!/usr/bin/env bash
# Run BUSCO on each filtered proteome and generate a summary table.
# Lineage: embryophyta_odb10 (for land plants)
# For Chlamydomonas (green alga), also run chlorophyta_odb10 separately.
#
# Requirements: busco on PATH (provided by the liverwort-rca conda env).
# Override the executable with BUSCO=/path/to/busco if needed.
#
# Usage: bash 03_busco_qc.sh [threads]

set -euo pipefail

THREADS=${1:-8}
FILTERED_DIR="data/proteomes/filtered"
BUSCO_DIR="data/busco"
LINEAGE_LAND="embryophyta_odb10"
LINEAGE_ALGAE="chlorophyta_odb10"

mkdir -p "$BUSCO_DIR"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

run_busco() {
    local fasta=$1
    local species=$2
    local lineage=$3
    local outdir="${BUSCO_DIR}/${species}"

    if [[ -d "${outdir}" ]]; then
        log "Skipping ${species} — BUSCO output already exists"
        return
    fi

    log "Running BUSCO for ${species} (lineage: ${lineage})..."
    "${BUSCO:-busco}" \
        --in "$fasta" \
        --out "$outdir" \
        --out_path "$BUSCO_DIR" \
        --mode proteins \
        --lineage_dataset "$lineage" \
        --download_path "$(pwd)/data/busco/lineages" \
        --cpu "$THREADS" \
        --offline \
        --quiet
}

# Land plants — embryophyta lineage
for fasta in "${FILTERED_DIR}"/*.primary.fa; do
    species=$(basename "$fasta" .primary.fa)
    if [[ "$species" == "Creinhardtii" ]]; then
        run_busco "$fasta" "$species" "$LINEAGE_ALGAE"
    else
        run_busco "$fasta" "$species" "$LINEAGE_LAND"
    fi
done

# ── Summary table ─────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════════════════════"
echo " BUSCO Summary — Metabolic Complexity Project"
echo "══════════════════════════════════════════════════════════════════"
printf "%-25s %8s %8s %8s %8s  %s\n" \
    "Species" "Complete" "Duplicat" "Fragment" "Missing" "Decision"
echo "──────────────────────────────────────────────────────────────────"

THRESHOLD=75  # % complete BUSCOs minimum to include species

for species_dir in "${BUSCO_DIR}"/*/; do
    species=$(basename "$species_dir")
    summary="${species_dir}/short_summary"*.txt
    summary=$(ls ${summary} 2>/dev/null | head -1)

    if [[ -z "$summary" ]]; then
        printf "%-25s %s\n" "$species" "No summary found"
        continue
    fi

    complete=$(grep -oE 'C:[0-9.]+%' "$summary" | head -1 | grep -oE '[0-9.]+')
    duplicated=$(grep -oE 'D:[0-9.]+%' "$summary" | head -1 | grep -oE '[0-9.]+')
    fragmented=$(grep -oE 'F:[0-9.]+%' "$summary" | head -1 | grep -oE '[0-9.]+')
    missing=$(grep -oE 'M:[0-9.]+%' "$summary" | head -1 | grep -oE '[0-9.]+')

    complete_int=${complete%.*}
    if (( complete_int >= THRESHOLD )); then
        decision="✅ INCLUDE"
    elif (( complete_int >= 60 )); then
        decision="⚠️  INCLUDE w/ caveat"
    else
        decision="❌ EXCLUDE"
    fi

    printf "%-25s %7s%% %7s%% %7s%% %7s%%  %s\n" \
        "$species" "$complete" "$duplicated" "$fragmented" "$missing" "$decision"
done

echo "──────────────────────────────────────────────────────────────────"
echo "Threshold: ≥${THRESHOLD}% complete BUSCOs required for inclusion"
echo ""
echo "Next step: review table, then run 04_run_orthofinder.sh"
