#!/bin/bash
# CP1 Deliverable A — run all 3 ambient tools on ONE dataset.
# Stop-on-first-failure (set -e) so per-dataset failures surface, not swallowed.
# Usage: cp1_run_dataset.sh <dataset_name>
set -euo pipefail
DS="$1"
BASE=/mnt/nvme1/omics-audit/phase2/audit_ambient_correction
GENE=/mnt/nvme1/omics-audit/phase2/audit3_counting/processed/$DS/star_default/Solo.out/Gene
RAW=$GENE/raw; FILT=$GENE/filtered
OUT=$BASE/cp1/per_dataset/$DS; mkdir -p "$OUT"
STAGE=$BASE/cp1/staging/${DS}_raw_v3; mkdir -p "$STAGE"
CB=$HOME/miniforge3/envs/ambient_cb/bin
LOG=$BASE/cp1/logs
echo "=== $DS START $(date +%H:%M:%S) ==="

# stage gzipped v3 layout for CellBender (CP0 finding: plain v3 misread as v2)
[ -s "$STAGE/matrix.mtx.gz" ]   || gzip -c "$RAW/matrix.mtx"   > "$STAGE/matrix.mtx.gz"
[ -s "$STAGE/features.tsv.gz" ] || gzip -c "$RAW/features.tsv" > "$STAGE/features.tsv.gz"
[ -s "$STAGE/barcodes.tsv.gz" ] || gzip -c "$RAW/barcodes.tsv" > "$STAGE/barcodes.tsv.gz"

# 1. CellBender (GPU, default parameters)
echo "[$DS] CellBender..."; cd "$OUT"
$CB/cellbender remove-background --input "$STAGE" --output cb.h5 --cuda \
    > "$LOG/${DS}_cellbender.log" 2>&1
$CB/python "$BASE/cp1/scripts/cp1_extract_cellbender.py" "$OUT/cb_filtered.h5" "$RAW" "$OUT/CellBender" \
    > "$LOG/${DS}_cb_extract.log" 2>&1

# 2. SoupX (Rscript subprocess)
echo "[$DS] SoupX..."
Rscript "$BASE/cp1/scripts/cp1_soupx.R" "$RAW" "$FILT" "$OUT/SoupX" \
    > "$LOG/${DS}_soupx.log" 2>&1

# 3. DecontX (Rscript subprocess)
echo "[$DS] DecontX..."
Rscript "$BASE/cp1/scripts/cp1_decontx.R" "$FILT" "$OUT/DecontX" \
    > "$LOG/${DS}_decontx.log" 2>&1

echo "=== $DS DONE $(date +%H:%M:%S) ==="
