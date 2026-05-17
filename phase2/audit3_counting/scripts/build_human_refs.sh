#!/usr/bin/env bash
# Build human references for all 3 binaries (4 configs).
#
# Inputs (already downloaded to REF_DIR):
#   GRCh38.primary_assembly.genome.fa.gz       — GENCODE v45 primary FASTA
#   gencode.v45.primary_assembly.annotation.gtf.gz  — GENCODE v45 GTF
#
# Outputs in REF_DIR/{star_index, salmon_splici_index, kallisto_index}:
#   STAR genome index            (slots 1 + 2: STARsolo default + CR-mimic)
#   salmon splici index          (slot 3: alevin-fry)
#   kallisto cDNA index + t2g    (slot 4: kb count)
#
# Resource expectations:
#   - STAR index: ~1-2 hr, ~30 GB peak RAM, ~30 GB on-disk index
#   - salmon splici index: ~20 min, ~10 GB peak RAM, ~10 GB on-disk
#   - kallisto cDNA index: ~10 min, ~5 GB peak RAM, ~3 GB on-disk
#
# Total time: ~2 hr serial, ~1-1.5 hr if STAR and salmon run concurrently.
# RAM budget on 64 GB machine: STAR alone hits ~30 GB. Run STAR FIRST and
# SEQUENTIALLY relative to the others to avoid OOM.
set -euo pipefail

THREADS="${THREADS:-12}"
REF_DIR="/mnt/nvme2/refs/GRCh38_gencode45"
ENV="/home/ross/miniforge3/envs/audit3_counting"
SJDB_OVERHANG="${SJDB_OVERHANG:-90}"  # 10x R2 = 90 cycles typical for v3.1

# Put env binaries on PATH so kb-python can resolve kallisto / bustools
export PATH="$ENV/bin:$PATH"

cd "$REF_DIR"

# Unpack FASTA + GTF (kept .gz too)
if [[ ! -f GRCh38.primary_assembly.genome.fa ]]; then
  echo "[$(date -Iseconds)] gunzip FASTA"
  gunzip -k GRCh38.primary_assembly.genome.fa.gz
fi
if [[ ! -f gencode.v45.primary_assembly.annotation.gtf ]]; then
  echo "[$(date -Iseconds)] gunzip GTF"
  gunzip -k gencode.v45.primary_assembly.annotation.gtf.gz
fi

# ---- STAR genome index (slots 1 + 2) ----
STAR_DIR="$REF_DIR/star_index_sjdb${SJDB_OVERHANG}"
if [[ ! -f "$STAR_DIR/SAindex" ]]; then
  echo "[$(date -Iseconds)] Building STAR index (sjdbOverhang=$SJDB_OVERHANG) — expect ~1-2 hr"
  mkdir -p "$STAR_DIR"
  "$ENV/bin/STAR" \
    --runMode genomeGenerate \
    --genomeDir "$STAR_DIR" \
    --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
    --sjdbGTFfile gencode.v45.primary_assembly.annotation.gtf \
    --sjdbOverhang "$SJDB_OVERHANG" \
    --runThreadN "$THREADS" \
    --genomeSAsparseD 3
  echo "[$(date -Iseconds)] STAR index DONE"
else
  echo "[$(date -Iseconds)] STAR index already exists — skipping"
fi

# ---- alevin-fry splici reference + salmon index (slot 3) ----
SPLICI_DIR="$REF_DIR/splici_fl91"
SALMON_DIR="$REF_DIR/salmon_splici_index"
if [[ ! -f "$SALMON_DIR/info.json" ]]; then
  echo "[$(date -Iseconds)] Building splici reference via pyroe — expect ~5 min"
  mkdir -p "$SPLICI_DIR"
  "$ENV/bin/pyroe" make-splici \
    GRCh38.primary_assembly.genome.fa \
    gencode.v45.primary_assembly.annotation.gtf \
    91 \
    "$SPLICI_DIR" \
    --flank-trim-length 5 \
    --filename-prefix splici
  echo "[$(date -Iseconds)] Building salmon index from splici — expect ~15 min"
  mkdir -p "$SALMON_DIR"
  "$ENV/bin/salmon" index \
    -t "$SPLICI_DIR/splici_fl86.fa" \
    -i "$SALMON_DIR" \
    -k 31 \
    --threads "$THREADS"
  cp "$SPLICI_DIR/splici_fl86_t2g_3col.tsv" "$SALMON_DIR/t2g_3col.tsv"
  echo "[$(date -Iseconds)] Salmon splici index DONE"
else
  echo "[$(date -Iseconds)] Salmon splici index already exists — skipping"
fi

# ---- kallisto-bustools index (slot 4) ----
KB_DIR="$REF_DIR/kb_index"
if [[ ! -f "$KB_DIR/index.idx" ]]; then
  echo "[$(date -Iseconds)] Building kallisto index via kb ref — expect ~10 min"
  mkdir -p "$KB_DIR"
  "$ENV/bin/kb" ref \
    -i "$KB_DIR/index.idx" \
    -g "$KB_DIR/t2g.txt" \
    -f1 "$KB_DIR/transcripts.fa" \
    --workflow standard \
    GRCh38.primary_assembly.genome.fa \
    gencode.v45.primary_assembly.annotation.gtf
  echo "[$(date -Iseconds)] kallisto index DONE"
else
  echo "[$(date -Iseconds)] kallisto index already exists — skipping"
fi

echo "[$(date -Iseconds)] All human references built. Sizes:"
du -sh "$STAR_DIR" "$SALMON_DIR" "$KB_DIR" 2>/dev/null || true
