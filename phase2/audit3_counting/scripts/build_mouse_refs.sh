#!/usr/bin/env bash
# Build the two mouse references not already cached on /mnt/nvme2/refs/mm39:
#   - alevin-fry splici reference + salmon index (cached salmon_index is plain cDNA, not splici)
#   - kallisto cDNA index (no cache)
#
# Reuses cached:
#   GRCm39.primary_assembly.genome.fa.gz
#   gencode.vM34.primary_assembly.annotation.gtf
#   star_index_sjdb150 (used as-is for STAR slots 1 + 2 — sjdb=150 vs 90 is
#                        sub-optimal for 10x R2=90 reads but held constant
#                        across all tool slots so audit isolates tool effect)
set -euo pipefail

THREADS="${THREADS:-12}"
REF_DIR="/mnt/nvme2/refs/mm39"
ENV="/home/ross/miniforge3/envs/audit3_counting"
export PATH="$ENV/bin:$PATH"

cd "$REF_DIR"

# Unpack FASTA if needed (cached GTF is already uncompressed)
if [[ ! -f GRCm39.primary_assembly.genome.fa ]]; then
  echo "[$(date -Iseconds)] gunzip mouse FASTA"
  gunzip -k GRCm39.primary_assembly.genome.fa.gz
fi

# ---- alevin-fry splici reference + salmon splici index ----
SPLICI_DIR="$REF_DIR/splici_fl91"
SALMON_SPLICI_DIR="$REF_DIR/salmon_splici_index"
if [[ ! -f "$SALMON_SPLICI_DIR/info.json" ]]; then
  echo "[$(date -Iseconds)] Building mouse splici reference via pyroe"
  mkdir -p "$SPLICI_DIR"
  pyroe make-splici \
    GRCm39.primary_assembly.genome.fa \
    gencode.vM34.primary_assembly.annotation.gtf \
    91 \
    "$SPLICI_DIR" \
    --flank-trim-length 5 \
    --filename-prefix splici
  echo "[$(date -Iseconds)] Building mouse salmon splici index"
  mkdir -p "$SALMON_SPLICI_DIR"
  salmon index \
    -t "$SPLICI_DIR/splici_fl86.fa" \
    -i "$SALMON_SPLICI_DIR" \
    -k 31 \
    --threads "$THREADS"
  cp "$SPLICI_DIR/splici_fl86_t2g_3col.tsv" "$SALMON_SPLICI_DIR/t2g_3col.tsv"
  echo "[$(date -Iseconds)] Mouse salmon splici index DONE"
else
  echo "[$(date -Iseconds)] Mouse salmon splici index already exists — skipping"
fi

# ---- kallisto-bustools index ----
KB_DIR="$REF_DIR/kb_index"
if [[ ! -f "$KB_DIR/index.idx" ]]; then
  echo "[$(date -Iseconds)] Building mouse kallisto index via kb ref"
  mkdir -p "$KB_DIR"
  kb ref \
    -i "$KB_DIR/index.idx" \
    -g "$KB_DIR/t2g.txt" \
    -f1 "$KB_DIR/transcripts.fa" \
    --workflow standard \
    GRCm39.primary_assembly.genome.fa \
    gencode.vM34.primary_assembly.annotation.gtf
  echo "[$(date -Iseconds)] Mouse kallisto index DONE"
else
  echo "[$(date -Iseconds)] Mouse kallisto index already exists — skipping"
fi

echo "[$(date -Iseconds)] Mouse references complete. Sizes:"
du -sh "$REF_DIR/star_index_sjdb150" "$SALMON_SPLICI_DIR" "$KB_DIR" 2>/dev/null || true
