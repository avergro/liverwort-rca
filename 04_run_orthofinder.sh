#!/usr/bin/env bash
# Run OrthoFinder on filtered proteomes that passed BUSCO QC.
#
# Input:  data/proteomes/filtered/*.primary.fa
# Output: data/orthofinder/
#
# OrthoFinder is on PATH (provided by the liverwort-rca conda env).
# Override the executable with ORTHOFINDER=/path/to/orthofinder if needed.
# Usage: bash 04_run_orthofinder.sh [threads]

set -euo pipefail

THREADS=${1:-8}
FILTERED="data/proteomes/filtered"
ORTHOFINDER_IN="data/orthofinder_v3_input"
ORTHOFINDER_OUT="data/orthofinder_v3"
ORTHOFINDER="${ORTHOFINDER:-orthofinder}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

# ── Build input directory ──────────────────────────────────────────────────
# OrthoFinder takes a directory where every .fa file = one species.
# Species excluded here based on BUSCO QC — edit list as needed.

EXCLUDE=()  # Add species codes to exclude, e.g. ("Pamericana" "Lcruciata")

mkdir -p "$ORTHOFINDER_IN"

for fa in "${FILTERED}"/*.primary.fa; do
    species=$(basename "$fa" .primary.fa)
    skip=false
    for excl in "${EXCLUDE[@]:-}"; do
        [[ "$species" == "$excl" ]] && skip=true && break
    done
    if [[ "$skip" == false ]]; then
        ln -sf "$(pwd)/$fa" "${ORTHOFINDER_IN}/${species}.fa"
        log "Included: $species"
    else
        log "Excluded: $species (manual exclusion)"
    fi
done

echo ""
log "Species included in OrthoFinder run:"
ls "${ORTHOFINDER_IN}/"
echo ""

# ── Run OrthoFinder ────────────────────────────────────────────────────────
# -S diamond_ultra_sens: best sensitivity for distant species (algae → plants)
# -M msa: gene tree inference using MSA (slower but better for deep phylogenies)
# -T fasttree: fast tree inference (use iqtree for publication-quality)

log "Starting OrthoFinder..."
"$ORTHOFINDER" \
    -f "$ORTHOFINDER_IN" \
    -o "$ORTHOFINDER_OUT" \
    -S diamond_ultra_sens \
    -M msa \
    -T fasttree \
    -t "$THREADS" \
    -a "$THREADS" \
    2>&1 | tee data/orthofinder_v3.log

log "OrthoFinder complete. Results in: ${ORTHOFINDER_OUT}/"
log "Key output files:"
log "  Orthogroups/Orthogroups.tsv          — gene counts per species per OG"
log "  Orthogroups/Orthogroups_SingleCopy.txt — single-copy OGs"
log "  Comparative_Genomics_Statistics/     — summary stats"
echo ""
log "Next step: python3 05_build_matrix.py"
