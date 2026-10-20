#!/usr/bin/env python3
"""RCA score (MCi-RCA) of every species in the four Gene Ontology contexts.

Reads the GO orthogroup filters of Matias (06_go_filter_multispecies.py) and the
same OrthoFinder output as the KEGG contexts, and writes
data/results_go/MCi_rca_GO.tsv (species x context) for Fig. 2c and Table 2.
Usage: python3 31_go_rankings.py
"""
import importlib.util
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
GO_DIR = ROOT / "Matias/alexander_metodologia/GO_files"
OUT = ROOT / "data/results_go"
CONTEXTS = {"primary": "primary_metabolism", "carbon": "C_metabolism",
            "nitrogen": "N_metabolism", "regulatory": "regulatory_networks"}

spec = importlib.util.spec_from_file_location("bm", ROOT / "05_build_matrix.py")
bm = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bm)

counts = bm.parse_orthogroups(bm.find_orthogroups_tsv(ROOT / "data/orthofinder_v3"))
scores, sizes = {}, {}
for ctx, f in CONTEXTS.items():
    ids = set((GO_DIR / f"{f}_ogs.txt").read_text().split())
    M = bm.filter_matrix(bm.build_rca_matrix(counts[counts.index.isin(ids)]), 2)
    scores[ctx] = bm.compute_complexity(M)[0]
    sizes[ctx] = (len(ids), M.shape[0])
OUT.mkdir(parents=True, exist_ok=True)
df = pd.DataFrame(scores)
df.to_csv(OUT / "MCi_rca_GO.tsv", sep="\t")
print(df.round(3).sort_values("primary", ascending=False))
for ctx, (n, m) in sizes.items():
    print(f"{ctx}: {n} orthogroups in the filter, {m} after RCA filtering")
