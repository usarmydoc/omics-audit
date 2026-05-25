#!/bin/bash
# CP1 Deliverable A — run remaining datasets sequentially, smallest->largest.
# Per-dataset PASS/FAIL recorded to RUN_STATUS.tsv (not swallowed); continues
# so one failure doesn't lose the others, but every failure is surfaced at end.
BASE=/mnt/nvme1/omics-audit/phase2/audit_ambient_correction/cp1
STATUS=$BASE/RUN_STATUS.tsv
echo -e "dataset\texit_code\tstart\tend" > "$STATUS"
# smallest -> largest by filtered cell count
DATASETS=(
  10x_neuron_1k_v3
  10x_t_3k_v2
  10x_pbmc_4k_v2
  10x_pbmc_5k_v3.1
  10x_pbmc_10k_v3.1
  gse325955_mouse_kidney_E18_5
  gse287209_human_lung_organoid
  gse288156_mouse_intestine_scrna
)
fails=0
for ds in "${DATASETS[@]}"; do
  st=$(date +%H:%M:%S)
  bash "$BASE/scripts/cp1_run_dataset.sh" "$ds" >> "$BASE/logs/run_all.log" 2>&1
  ec=$?
  en=$(date +%H:%M:%S)
  echo -e "${ds}\t${ec}\t${st}\t${en}" >> "$STATUS"
  if [ "$ec" -ne 0 ]; then echo "[FAIL] $ds (exit $ec) — see logs/${ds}_*.log" >> "$BASE/logs/run_all.log"; fails=$((fails+1)); fi
done
echo "ALL_DATASETS_DONE fails=$fails" >> "$BASE/logs/run_all.log"
echo "CP1_RUN_ALL_COMPLETE fails=$fails"
