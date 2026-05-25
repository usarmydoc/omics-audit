#!/bin/bash
# CP3 Deliverable C — biological propagation. pbmc_5k (control) then intestine.
# no-correction baseline reused from Audit 3 C3 (star_default). 5 corrected conditions
# generated + run through the identical cp6 downstream pipeline. Resumable (skip if
# condition summary exists). Matrix-gen SoupX∥DecontX concurrent; pipelines sequential.
set -uo pipefail
BASE=/mnt/nvme1/omics-audit/phase2/audit_ambient_correction/cp3
A3=/mnt/nvme1/omics-audit/phase2/audit3_counting
PROC=$A3/processed
CB=$HOME/miniforge3/envs/ambient_cb/bin/python
A3PY=$HOME/miniforge3/envs/audit3_counting/bin/python
LOG=$BASE/logs; ST=$BASE/RUN_STATUS.tsv
[ -f "$ST" ] || echo -e "dataset\tcondition\tstatus\ttime" > "$ST"

run_pipe(){ local ds=$1 cond=$2 sub=$3
  if [ -s "$BASE/$sub/${cond}.summary.json" ]; then echo "skip $sub/$cond"; return 0; fi
  if $A3PY "$BASE/scripts/cp3_pipeline.py" "$ds" "$cond" "$BASE/matrices/$ds/$cond" > "$LOG/${sub}_${cond}.log" 2>&1; then
    echo -e "$ds\t$cond\tDONE\t$(date +%H:%M:%S)" >> "$ST"; echo "  [$sub/$cond] DONE"
  else echo -e "$ds\t$cond\tFAIL\t$(date +%H:%M:%S)" >> "$ST"; echo "  [$sub/$cond] FAIL — see $LOG/${sub}_${cond}.log"; return 1; fi
}

declare -A DSMAP=( [pbmc]=10x_pbmc_5k_v3.1 [intestine]=gse288156_mouse_intestine_scrna )
declare -A C3SUB=( [pbmc]=control/10x_pbmc_5k_v3.1 [intestine]=intestine )
for sub in pbmc intestine; do
  ds=${DSMAP[$sub]}; G=$PROC/$ds/star_default/Solo.out/Gene
  echo "=== $sub ($ds) $(date +%H:%M:%S) ==="
  mkdir -p "$BASE/$sub" "$BASE/matrices/$ds"

  # no-correction baseline: reuse C3 star_default outputs
  c3="$A3/c3/${C3SUB[$sub]}/per_tool_pipeline_outputs"
  for ext in obs.tsv markers.tsv summary.json; do cp -n "$c3/star_default.$ext" "$BASE/$sub/nocorr.$ext" 2>/dev/null; done

  # generate condition matrices
  [ -s "$BASE/matrices/$ds/CellBender/matrix.mtx" ] || \
    $CB "$BASE/scripts/cp3_cellbender_matrix.py" \
      "/mnt/nvme1/omics-audit/phase2/audit_ambient_correction/cp1/per_dataset/$ds/cb_filtered.h5" \
      "$BASE/matrices/$ds/CellBender" > "$LOG/${sub}_genCB.log" 2>&1
  for ord in O1 O2; do
    declare -A P=()
    for tool in SoupX DecontX; do
      cond="${tool}_${ord}"
      if [ ! -s "$BASE/matrices/$ds/$cond/matrix.mtx" ]; then
        ( Rscript "$BASE/scripts/cp3_genmatrix.R" "$tool" "$G/raw" "$G/filtered" "$BASE/matrices/$ds/$cond" "$ord" \
            > "$LOG/${sub}_gen_${cond}.log" 2>&1 ) & P[$cond]=$!
      fi
    done
    for c in "${!P[@]}"; do wait "${P[$c]}" || echo "  gen $c FAILED"; done
  done

  # run pipelines sequentially (RAM-bound)
  for cond in CellBender SoupX_O1 SoupX_O2 DecontX_O1 DecontX_O2; do run_pipe "$ds" "$cond" "$sub" || true; done
done
echo "CP3_ORCHESTRATE_COMPLETE"
