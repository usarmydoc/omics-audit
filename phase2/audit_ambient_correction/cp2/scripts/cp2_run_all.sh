#!/bin/bash
# CP2 Deliverable B — all datasets, sequential smallest->largest.
# Resumable (per-dataset runner skips existing (tool,ordering) summaries).
# Per-dataset PASS/FAIL recorded. Halts if 2 CONSECUTIVE datasets fail.
BASE=/mnt/nvme1/omics-audit/phase2/audit_ambient_correction/cp2
STATUS=$BASE/RUN_STATUS.tsv
[ -f "$STATUS" ] || echo -e "dataset\texit_code\tstart\tend" > "$STATUS"
DATASETS=(
  10x_pbmc_1k_v3 10x_neuron_1k_v3 10x_t_3k_v2 10x_pbmc_4k_v2 10x_pbmc_5k_v3.1
  10x_pbmc_10k_v3.1 gse325955_mouse_kidney_E18_5 gse287209_human_lung_organoid
  gse288156_mouse_intestine_scrna
)
consec=0; fails=0
for ds in "${DATASETS[@]}"; do
  st=$(date +%H:%M:%S)
  bash "$BASE/scripts/cp2_run_dataset.sh" "$ds" >> "$BASE/logs/run_all.log" 2>&1
  ec=$?
  en=$(date +%H:%M:%S)
  echo -e "${ds}\t${ec}\t${st}\t${en}" >> "$STATUS"
  if [ "$ec" -ne 0 ]; then
    echo "[FAIL] $ds (exit $ec)" >> "$BASE/logs/run_all.log"; fails=$((fails+1)); consec=$((consec+1))
    if [ "$consec" -ge 2 ]; then echo "HALT: 2 consecutive dataset failures" >> "$BASE/logs/run_all.log"; echo "CP2_HALTED_2_CONSECUTIVE_FAILS"; exit 1; fi
  else consec=0; fi
done
echo "ALL_DATASETS_DONE fails=$fails" >> "$BASE/logs/run_all.log"
echo "CP2_RUN_ALL_COMPLETE fails=$fails"
