#!/bin/bash
# CP1 — QC methods + 12 downstream pipeline runs + compare.
# Assumes matrices/<short>_raw.h5ad already pulled (cp1_pull.py, base env).
# Per-condition PASS/FAIL recorded; continue on single failure; halt on 2+ same error.
set -uo pipefail
CP1=/mnt/nvme1/omics-audit/phase2/audit_qc_mad_propagation/cp1
PY=$HOME/miniforge3/envs/audit3_counting/bin/python
LOG=$CP1/logs; ST=$CP1/RUN_STATUS.tsv
echo -e "dataset\tmethod\tstatus\ttime" > "$ST"
DATASETS=(blood liver small_intestine)   # clean control first, then edge-case
METHODS=(C1 C2 MAD3 MAD5)

# Step 2: QC methods per dataset (verify cell counts vs QC-MAD)
for ds in "${DATASETS[@]}"; do
  echo "=== QC methods: $ds $(date +%H:%M:%S) ==="
  if ! $PY "$CP1/scripts/cp1_qc.py" "$ds" > "$LOG/${ds}_qc.log" 2>&1; then
    echo "[FATAL] $ds QC method counts mismatch QC-MAD — see ${ds}_qc.log; HALT"; exit 1
  fi
  grep -E "vs|reproduce" "$LOG/${ds}_qc.log" | tail -5
done

# Step 3: 12 pipeline runs (sequential, RAM-bound)
declare -A errsig; halt=0
for ds in "${DATASETS[@]}"; do
  for m in "${METHODS[@]}"; do
    [ "$halt" -eq 1 ] && break 2
    if [ -s "$CP1/per_condition/$ds/${m}.summary.json" ]; then echo "skip $ds/$m"; echo -e "$ds\t$m\tSKIP\t$(date +%H:%M:%S)" >>"$ST"; continue; fi
    if $PY "$CP1/scripts/cp1_pipeline.py" "$ds" "$m" > "$LOG/${ds}_${m}.log" 2>&1; then
      echo -e "$ds\t$m\tDONE\t$(date +%H:%M:%S)" >>"$ST"; echo "  [$ds/$m] DONE"
    else
      sig=$(grep -oE "Error[:A-Za-z ]*|RuntimeError|KeyError|ValueError" "$LOG/${ds}_${m}.log" | tail -1)
      echo -e "$ds\t$m\tFAIL\t$(date +%H:%M:%S)" >>"$ST"; echo "  [$ds/$m] FAIL ($sig) — see ${ds}_${m}.log"
      errsig[$sig]=$(( ${errsig[$sig]:-0} + 1 ))
      if [ "${errsig[$sig]}" -ge 2 ]; then echo "HALT: 2+ conditions failed with same error: $sig"; halt=1; fi
    fi
  done
done
[ "$halt" -eq 1 ] && { echo "CP1_HALTED_REPEATED_ERROR"; exit 1; }

# Step 4-5: compare
$PY "$CP1/scripts/cp1_compare.py" "${DATASETS[@]}" > "$LOG/compare.log" 2>&1
echo "CP1_COMPLETE"; tail -20 "$LOG/compare.log"
