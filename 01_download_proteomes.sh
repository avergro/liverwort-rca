#!/usr/bin/env bash
# Download proteomes for metabolic complexity analysis
# Verified filenames as of 2026-05-03 (Ensembl Plants release 62)
#
# Run from project root. Creates data/proteomes/raw/

set -euo pipefail

OUTDIR="data/proteomes/raw"
mkdir -p "$OUTDIR"

ENSEMBL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-62/fasta"
FERNBASE="https://fernbase.org/ftp/ftp"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

dl() {
    local url=$1
    local out=$2
    if [[ -f "$out" ]]; then
        log "Already exists, skipping: $(basename "$out")"
        return
    fi
    log "Downloading: $(basename "$out")"
    wget -q --show-progress -O "$out" "$url"
}

# ── Ensembl Plants release 62 ──────────────────────────────────────────────

dl "${ENSEMBL}/marchantia_polymorpha/pep/Marchantia_polymorpha.MpTak1_v7.1.pep.all.fa.gz" \
   "${OUTDIR}/Mpolymorpha.pep.all.fa.gz"

dl "${ENSEMBL}/physcomitrium_patens/pep/Physcomitrium_patens.Phypa_V3.pep.all.fa.gz" \
   "${OUTDIR}/Ppatens.pep.all.fa.gz"

dl "${ENSEMBL}/arabidopsis_thaliana/pep/Arabidopsis_thaliana.TAIR10.pep.all.fa.gz" \
   "${OUTDIR}/Athaliana.pep.all.fa.gz"

dl "${ENSEMBL}/chlamydomonas_reinhardtii/pep/Chlamydomonas_reinhardtii.Chlamydomonas_reinhardtii_v5.5.pep.all.fa.gz" \
   "${OUTDIR}/Creinhardtii.pep.all.fa.gz"

dl "${ENSEMBL}/solanum_lycopersicum/pep/Solanum_lycopersicum.SL3.0.pep.all.fa.gz" \
   "${OUTDIR}/Slycopersicum.pep.all.fa.gz"

# ── FernBase ───────────────────────────────────────────────────────────────
# High-confidence only — low-confidence excluded by design decision

dl "${FERNBASE}/Azolla_filiculoides/Azolla_asm_v1.0/Azolla_filiculoides.protein.highconfidence_v1.fasta" \
   "${OUTDIR}/Afiliculoides.pep.fa"

dl "${FERNBASE}/Salvinia_cucullata/Salvinia_asm_v1.2/Salvinia_cucullata.protein.highconfidence_v1.2.fasta" \
   "${OUTDIR}/Scucullata.pep.fa"

# ── Persea americana (avocado) — NCBI RefSeq GCF_002303995.2 ──────────────
# Quality uncertain — subject to BUSCO QC before inclusion in OrthoFinder

dl "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/029/852/735/GCA_029852735.1_ASM2985273v1/GCA_029852735.1_ASM2985273v1_protein.faa.gz" \
   "${OUTDIR}/Pamericana.pep.fa.gz"

# ── Picea abies V2 — local files (ConGenIE / Nat. Plants 2024) ────────────
# Source: <Picea abies V2 files>/
# File: Picab02_230926_at01_all_aa.fa.txt (181,046 isoforms, 54,510 genes)
# Pre-filtered longest isoform: Picab02_230926_at01_longest_aa.fa.txt
if [[ ! -L "${OUTDIR}/Pabies.pep.fa" ]]; then
    ln -sf <Picea abies V2 files>/Picab02_230926_at01_all_aa.fa.txt \
       "${OUTDIR}/Pabies.pep.fa"
    log "Pabies: symlink created"
else
    log "Pabies: already present"
fi

# ── Selaginella moellendorffii — Ensembl Plants r62 ──────────────────────
# Lycophyte outgroup — phylogenetic bridge between ferns and angiosperms

dl "${ENSEMBL}/selaginella_moellendorffii/pep/Selaginella_moellendorffii.v1.0.pep.all.fa.gz" \
   "${OUTDIR}/Smoellendorffii.pep.all.fa.gz"

# ── Lunularia cruciata — manual download required ─────────────────────────
if [[ ! -f "${OUTDIR}/Lcruciata.pep.fa" ]]; then
    cat <<'MSG'

──────────────────────────────────────────────────────────────────────────────
Lunularia cruciata — MANUAL DOWNLOAD REQUIRED
──────────────────────────────────────────────────────────────────────────────
Source: Linde et al. 2023, Genome Biology and Evolution (doi:10.1093/gbe/evad014)

Option A — CoGe web interface:
  1. Go to https://genomevolution.org/coge/GenomeInfo.pl?gid=64629
  2. Export → Protein sequences (CDS translation)
  3. Save as: data/proteomes/raw/Lcruciata.pep.fa

Option B — Check NCBI BioProject PRJNA927668 (same paper)
  wget -O data/proteomes/raw/Lcruciata.pep.fa <NCBI_URL>
──────────────────────────────────────────────────────────────────────────────
MSG
else
    log "Lunularia cruciata already present — OK"
fi

# ── Summary ────────────────────────────────────────────────────────────────
echo ""
log "Download complete. Files in ${OUTDIR}/:"
ls -lh "${OUTDIR}/"
echo ""
log "Next step: python3 02_filter_primary_transcripts.py"
