#!/bin/bash
# Run MCi analysis at 4 levels:
#   1. Whole genome
#   2. Primary metabolism (all KEGG metabolic pathways)
#   3. C metabolism (carbohydrate + energy)
#   4. N metabolism (amino acid + nitrogen assimilation)
#
# Prerequisites:
#   - OrthoFinder results in data/orthofinder_v3/
#   - Run 06_kegg_filter.py first to generate data/filters/

set -euo pipefail

PYTHON="${PYTHON:-python3}"
OFD=data/orthofinder_v3

echo "============================================================"
echo " Step 1: Generate KEGG functional filters"
echo "============================================================"
$PYTHON 06_kegg_filter.py --orthofinder-dir $OFD --out-dir data/filters

echo ""
echo "============================================================"
echo " Step 2: MCi — whole genome"
echo "============================================================"
$PYTHON 05_build_matrix.py \
    --orthofinder-dir $OFD \
    --out-dir data/results_v3

echo ""
echo "============================================================"
echo " Step 3: MCi — primary metabolism"
echo "============================================================"
$PYTHON 05_build_matrix.py \
    --orthofinder-dir $OFD \
    --og-filter data/filters/primary_metabolism_ogs.txt \
    --out-dir data/results_primary_metabolism_v3

echo ""
echo "============================================================"
echo " Step 4: MCi — C metabolism"
echo "============================================================"
$PYTHON 05_build_matrix.py \
    --orthofinder-dir $OFD \
    --og-filter data/filters/C_metabolism_ogs.txt \
    --out-dir data/results_C_metabolism_v3

echo ""
echo "============================================================"
echo " Step 5: MCi — N metabolism"
echo "============================================================"
$PYTHON 05_build_matrix.py \
    --orthofinder-dir $OFD \
    --og-filter data/filters/N_metabolism_ogs.txt \
    --out-dir data/results_N_metabolism_v3

echo ""
echo "============================================================"
echo " All analyses complete. Results summary:"
echo "============================================================"
for dir in data/results_v3 data/results_primary_metabolism_v3 data/results_C_metabolism_v3 data/results_N_metabolism_v3; do
    if [ -f "$dir/species_complexity.tsv" ]; then
        echo ""
        echo "  [$dir]"
        $PYTHON -c "
import pandas as pd
df = pd.read_csv('$dir/species_complexity.tsv', sep='\t', index_col=0)
df = df.sort_values('MCi_rca', ascending=False)
print(df[['MCi_rca']].to_string())
"
    fi
done
