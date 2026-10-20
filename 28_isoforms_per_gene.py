#!/usr/bin/env python3
"""Isoforms per protein-coding gene in the eleven raw proteomes.

Counts protein sequences and genes in data/proteomes/raw (before the
longest-isoform filter) and reports isoforms per gene. Gene identity is
taken from the 'gene:' / 'gene=' header tag when present; otherwise the
header ID is collapsed by the species-specific isoform suffix rule below.
Output: data/results_isoforms/isoforms_per_gene.tsv
"""
import gzip, re
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RAW = ROOT / "data/proteomes/raw"
FILT = ROOT / "data/proteomes/filtered"
OUT = ROOT / "data/results_isoforms"

# when no gene tag exists: strip an isoform suffix from the sequence ID
SUFFIX = {
    "Afiliculoides": r"\.\d+$",      # Azfi_s0001.g000001(.1)
    "Scucullata": r"\.\d+$",
    "Lcruciata": None,               # one sequence per gene ID
    "Pamericana": None,              # NCBI protein accessions, gene not in header
}

def opener(f):
    return gzip.open(f, "rt") if f.suffix == ".gz" else open(f)

rows = []
for f in sorted(RAW.iterdir()):
    sp = f.name.split(".")[0]
    seqs, genes = 0, Counter()
    with opener(f) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            seqs += 1
            m = re.search(r"gene[:=](\S+)", line)
            if m:
                g = m.group(1)
            else:
                g = line[1:].split()[0]
                if SUFFIX.get(sp):
                    g = re.sub(SUFFIX[sp], "", g)
            genes[g] += 1
    n_filt = sum(1 for l in open(FILT / f"{sp}.primary.fa") if l.startswith(">"))
    multi = sum(1 for c in genes.values() if c > 1)
    rows.append((sp, seqs, len(genes), n_filt, seqs / len(genes), multi, 100 * multi / len(genes), max(genes.values())))

OUT.mkdir(exist_ok=True)
with open(OUT / "isoforms_per_gene.tsv", "w") as out:
    out.write("species\tprotein_sequences\tgenes\tgenes_after_filter\tisoforms_per_gene\tgenes_with_>1_isoform\tpct_genes_multi_isoform\tmax_isoforms\n")
    for r in sorted(rows, key=lambda r: -r[4]):
        out.write("%s\t%d\t%d\t%d\t%.2f\t%d\t%.1f\t%d\n" % r)
print(open(OUT / "isoforms_per_gene.tsv").read())
