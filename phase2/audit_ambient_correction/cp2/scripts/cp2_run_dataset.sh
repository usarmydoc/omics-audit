#!/bin/bash
# CP2 Deliverable B — run all 3 tools x 2 orderings for ONE dataset.
# Optimized per system: SoupX/DecontX both orderings run CONCURRENTLY (CPU,
# orthogonal); CellBender reuses CP1 corrected output (no GPU re-run), QC-only
# rerun per ordering. Resumable: skips a (tool,ordering) whose summary exists.
set -uo pipefail
DS="$1"
BASE=/mnt/nvme1/omics-audit/phase2/audit_ambient_correction
GENE=/mnt/nvme1/omics-audit/phase2/audit3_counting/processed/$DS/star_default/Solo.out/Gene
RAW=$GENE/raw; FILT=$GENE/filtered
CB_H5=$BASE/cp1/per_dataset/$DS/cb_filtered.h5
OUT=$BASE/cp2/per_dataset/$DS; mkdir -p "$OUT"
LOG=$BASE/cp2/logs; SC=$BASE/cp2/scripts
CBPY=$HOME/miniforge3/envs/ambient_cb/bin/python
echo "=== $DS START $(date +%H:%M:%S) ==="
rc=0
need(){ [ ! -s "$OUT/${1}_summary.json" ]; }   # resume: rerun only if summary missing

# --- 4 R runs concurrently (SoupX/DecontX x O1/O2) ---
declare -A PID
for tool in SoupX DecontX; do for ord in O1 O2; do
  tag="${tool}_${ord}"
  if need "$tag"; then
    if [ "$tool" = SoupX ]; then
      ( Rscript "$SC/cp2_soupx.R" "$RAW" "$FILT" "$OUT/$tag" "$ord" > "$LOG/${DS}_${tag}.log" 2>&1 ) &
    else
      ( Rscript "$SC/cp2_decontx.R" "$FILT" "$OUT/$tag" "$ord" > "$LOG/${DS}_${tag}.log" 2>&1 ) &
    fi
    PID[$tag]=$!
  else echo "[$DS] skip $tag (exists)"; fi
done; done
for tag in "${!PID[@]}"; do
  if ! wait "${PID[$tag]}"; then echo "[$DS] FAIL $tag — see logs/${DS}_${tag}.log"; rc=1; fi
done

# --- CellBender both orderings (reuse CP1; QC-only; sequential, light) ---
for ord in O1 O2; do
  tag="CellBender_${ord}"
  if need "$tag"; then
    if ! $CBPY "$SC/cp2_cellbender.py" "$CB_H5" "$RAW" "$FILT/features.tsv" "$OUT/$tag" "$ord" \
         > "$LOG/${DS}_${tag}.log" 2>&1; then echo "[$DS] FAIL $tag"; rc=1; fi
  else echo "[$DS] skip $tag (exists)"; fi
done

echo "=== $DS END $(date +%H:%M:%S) rc=$rc ==="
exit $rc
