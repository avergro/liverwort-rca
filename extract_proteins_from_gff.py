#!/usr/bin/env python3
"""
Extract protein sequences from genome FASTA + GFF3 annotation (CoGe format).

Usage:
    python3 extract_proteins_from_gff.py \
        --genome  data/proteomes/raw/Lcruciata.genome.fa \
        --gff     data/proteomes/raw/Lcruciata.gff \
        --out     data/proteomes/raw/Lcruciata.pep.fa
"""

import argparse
import sys
from collections import defaultdict

CODON_TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

RC = str.maketrans('ACGTacgt', 'TGCAtgca')

def revcomp(seq):
    return seq.translate(RC)[::-1]

def translate(nuc):
    prot = []
    for i in range(0, len(nuc) - 2, 3):
        codon = nuc[i:i+3].upper()
        aa = CODON_TABLE.get(codon, 'X')
        if aa == '*':
            break
        prot.append(aa)
    return ''.join(prot)

def load_genome(fa_path):
    print(f"Loading genome: {fa_path}", flush=True)
    genome = {}
    seq_id = None
    parts = []
    with open(fa_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if seq_id:
                    genome[seq_id] = ''.join(parts)
                seq_id = line[1:].split()[0].removeprefix('lcl|')
                parts = []
            else:
                parts.append(line)
    if seq_id:
        genome[seq_id] = ''.join(parts)
    print(f"  {len(genome):,} scaffolds loaded", flush=True)
    return genome

def parse_gff(gff_path):
    """Returns dict: mrna_id -> {'chrom','strand','cds':[(start,end), ...],'name':gene_name}"""
    mrna_info = {}   # mrna_id → {chrom, strand, name}
    cds_map   = defaultdict(list)  # mrna_id → [(start1,end1), ...]

    with open(gff_path) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            parts = line.rstrip().split('\t')
            if len(parts) < 9:
                continue
            feat_type = parts[2]
            if feat_type == 'mRNA':
                chrom  = parts[0]
                strand = parts[6]
                attrs  = dict(kv.split('=', 1) for kv in parts[8].split(';') if '=' in kv)
                mrna_id = attrs.get('ID', '')
                name    = attrs.get('Name', mrna_id)
                mrna_info[mrna_id] = {'chrom': chrom, 'strand': strand, 'name': name}
            elif feat_type == 'CDS':
                start = int(parts[3])
                end   = int(parts[4])
                attrs = dict(kv.split('=', 1) for kv in parts[8].split(';') if '=' in kv)
                parent = attrs.get('Parent', '')
                cds_map[parent].append((start, end))

    return mrna_info, cds_map

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--genome', required=True)
    ap.add_argument('--gff',    required=True)
    ap.add_argument('--out',    required=True)
    args = ap.parse_args()

    genome = load_genome(args.genome)
    print(f"Parsing GFF: {args.gff}", flush=True)
    mrna_info, cds_map = parse_gff(args.gff)
    print(f"  {len(mrna_info):,} mRNAs, {len(cds_map):,} mRNAs with CDS", flush=True)

    written = skipped_missing = skipped_short = 0

    with open(args.out, 'w') as out:
        for mrna_id, info in mrna_info.items():
            cdss = cds_map.get(mrna_id)
            if not cdss:
                skipped_missing += 1
                continue

            chrom  = info['chrom']
            strand = info['strand']
            name   = info['name']

            if chrom not in genome:
                skipped_missing += 1
                continue

            seq = genome[chrom]
            cdss_sorted = sorted(cdss, key=lambda x: x[0])

            nuc = ''.join(seq[s-1:e] for s, e in cdss_sorted)  # GFF is 1-based inclusive
            if strand == '-':
                nuc = revcomp(nuc)

            prot = translate(nuc)
            if len(prot) < 30:
                skipped_short += 1
                continue

            out.write(f'>{name}\n')
            for i in range(0, len(prot), 60):
                out.write(prot[i:i+60] + '\n')
            written += 1

    print(f"\nDone: {written:,} proteins written")
    print(f"      {skipped_missing:,} skipped (missing scaffold or no CDS)")
    print(f"      {skipped_short:,}  skipped (<30 aa)")
    print(f"Output: {args.out}")

if __name__ == '__main__':
    main()
