#!/usr/bin/env python3
"""
Filter proteome FASTAs to one representative protein per gene.
Strategy: keep the longest isoform per gene locus.

Handles header formats from:
  - Ensembl Plants:  >PROTEINID ... gene:GENEID ...
  - FernBase:        >Azfi_s001.g000001.1 (gene inferred from dropping last .N)
  - NCBI RefSeq:     >XP_... gene=GENENAME ...
  - CoGe (Lunularia): >gene_id.isoform or similar

Output: data/proteomes/filtered/<species>.primary.fa
"""

import gzip
import re
import sys
from pathlib import Path
from collections import defaultdict


def open_fasta(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_gene_id(header: str, source_hint: str) -> str:
    """Extract gene locus ID from a FASTA header line."""

    # Ensembl: gene:GENEID in header
    m = re.search(r'gene:(\S+)', header)
    if m:
        return m.group(1)

    # NCBI RefSeq/Genbank annotated: gene=GENENAME
    m = re.search(r'gene=(\S+)', header)
    if m:
        return m.group(1)

    # NCBI GenBank WGS locus_tag pattern: e.g. "hypothetical protein ABC123_042440"
    # Locus tags follow format PREFIX_DIGITS (e.g. MRB53_042440, Mapoly0001s0001)
    m = re.search(r'\b([A-Za-z0-9]+_\d{6,})\b', header)
    if m:
        return m.group(1)

    # FernBase / CoGe: strip last .N from protein ID
    protein_id = header.split()[0].lstrip('>')
    gene_id = re.sub(r'\.\d+$', '', protein_id)
    return gene_id


def read_fasta(path: Path) -> dict[str, str]:
    """Return {header: sequence} preserving full headers."""
    records = {}
    current_header = None
    current_seq = []
    with open_fasta(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if current_header:
                    records[current_header] = ''.join(current_seq)
                current_header = line
                current_seq = []
            elif line:
                current_seq.append(line)
        if current_header:
            records[current_header] = ''.join(current_seq)
    return records


def filter_longest_isoform(records: dict[str, str]) -> dict[str, tuple[str, str]]:
    """
    Return {gene_id: (header, sequence)} keeping longest isoform per gene.
    """
    best: dict[str, tuple[str, str, int]] = {}  # gene_id -> (header, seq, length)
    for header, seq in records.items():
        gene_id = parse_gene_id(header, "")
        seq_len = len(seq.replace('*', ''))  # ignore stop codon char
        if gene_id not in best or seq_len > best[gene_id][2]:
            best[gene_id] = (header, seq, seq_len)
    return {gid: (h, s) for gid, (h, s, _) in best.items()}


def write_fasta(filtered: dict[str, tuple[str, str]], out_path: Path) -> int:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    count = 0
    with open(out_path, 'w') as fh:
        for gene_id, (header, seq) in filtered.items():
            fh.write(f"{header}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
            count += 1
    return count


def main():
    raw_dir = Path("data/proteomes/raw")
    out_dir = Path("data/proteomes/filtered")

    extensions = [".pep.all.fa.gz", ".pep.fa", ".pep.fa.gz", ".faa.gz", ".fa"]

    inputs = []
    for ext in extensions:
        inputs.extend(raw_dir.glob(f"*{ext}"))
    inputs = list(set(inputs))  # deduplicate

    if not inputs:
        print(f"No proteome files found in {raw_dir}", file=sys.stderr)
        sys.exit(1)

    print(f"{'Species':<25} {'Raw proteins':>14} {'Genes kept':>12} {'Reduction':>10}")
    print("-" * 65)

    for fasta_path in sorted(inputs):
        species = fasta_path.name.split('.')[0]
        records = read_fasta(fasta_path)
        filtered = filter_longest_isoform(records)
        out_path = out_dir / f"{species}.primary.fa"
        n_written = write_fasta(filtered, out_path)
        reduction = (1 - n_written / len(records)) * 100 if records else 0
        print(f"{species:<25} {len(records):>14,} {n_written:>12,} {reduction:>9.1f}%")

    print(f"\nFiltered FASTAs written to: {out_dir}/")
    print("Next step: run 03_busco_qc.sh on each filtered FASTA")


if __name__ == "__main__":
    main()
