#!/bin/bash
# CP1 (parallel) — QC methods + 12 pipeline runs CONCURRENTLY (cap 3, RAM-bounded
# ~9.5GB each) + serial repro.lock registration + compare. Resumable (skip existing
# summaries). Per-condition PASS/FAIL recorded. Uses the 24-thread / 60GB machine.
set -uo pipefail
CP1=/mnt/nvme1/omics-audit/phase2/audit_qc_mad_propagation/cp1
PY=$HOME/miniforge3/envs/audit3_counting/bin/python
LOG=$CP1/logs; ST=$CP1/RUN_STATUS.tsv; MAXJOBS=3
# CAP per-run threads so MAXJOBS concurrent x N threads <= 24 cores (avoid the
# load-56 oversubscription from uncapped scanpy/BLAS/numba grabbing all cores each).
# 3 jobs x 7 threads = 21. Exported -> inherited by python AND its Rscript subprocs.
export OMP_NUM_THREADS=7 OPENBLAS_NUM_THREADS=7 MKL_NUM_THREADS=7 NUMBA_NUM_THREADS=7 VECLIB_MAXIMUM_THREADS=7
echo -e "dataset\tmethod\tstatus\ttime" > "$ST"
DATASETS=(blood liver small_intestine); METHODS=(C1 C2 MAD3 MAD5)

# Step 2: QC methods (skip if survivors already written)
for ds in "${DATASETS[@]}"; do
  if [ -s "$CP1/survivors/${ds}_MAD5.txt" ]; then echo "QC methods $ds: cached"; continue; fi
  echo "=== QC methods: $ds ==="
  $PY "$CP1/scripts/cp1_qc.py" "$ds" > "$LOG/${ds}_qc.log" 2>&1 || { echo "[FATAL] $ds QC mismatch"; exit 1; }
done

# Step 3: 12 pipeline runs, up to MAXJOBS concurrent
run_one(){ local ds=$1 m=$2
  if [ -s "$CP1/per_condition/$ds/${m}.summary.json" ]; then echo -e "$ds\t$m\tSKIP\t$(date +%H:%M:%S)" >>"$ST"; echo "skip $ds/$m"; return; fi
  if $PY "$CP1/scripts/cp1_pipeline.py" "$ds" "$m" > "$LOG/${ds}_${m}.log" 2>&1; then
    echo -e "$ds\t$m\tDONE\t$(date +%H:%M:%S)" >>"$ST"; echo "  [$ds/$m] DONE"
  else echo -e "$ds\t$m\tFAIL\t$(date +%H:%M:%S)" >>"$ST"; echo "  [$ds/$m] FAIL — see ${ds}_${m}.log"; fi
}
for ds in "${DATASETS[@]}"; do for m in "${METHODS[@]}"; do
  run_one "$ds" "$m" &
  while [ "$(jobs -rp | wc -l)" -ge "$MAXJOBS" ]; do wait -n; done
done; done
wait
nfail=$(grep -c $'\tFAIL\t' "$ST" || true)
echo "pipeline batch done; failures=$nfail"

# Step 4: SERIAL repro.lock registration (avoids parallel write race)
$PY - <<'PYREG'
import sys, glob; sys.path.insert(0,"/mnt/nvme1/omics-audit/phase2/scripts")
import dge_native as dn
n=0
for p in sorted(glob.glob("/mnt/nvme1/omics-audit/phase2/audit_qc_mad_propagation/cp1/per_condition/*/*")):
    if p.endswith((".obs.tsv",".markers.tsv",".summary.json")):
        ds=p.split("/per_condition/")[1].split("/")[0]
        dn.register_output(__import__("pathlib").Path(p), kind="qcprop_cp1", dataset=ds); n+=1
print(f"registered {n} CP1 condition outputs (serial)")
PYREG

# Step 5: compare
$PY "$CP1/scripts/cp1_compare.py" "${DATASETS[@]}" > "$LOG/compare.log" 2>&1
echo "CP1_COMPLETE failures=$nfail"; tail -18 "$LOG/compare.log"
