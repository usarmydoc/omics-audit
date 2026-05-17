#!/usr/bin/env bash
# Process one dataset across all 4 counting tool configurations.
#
# Usage: process_one_dataset.sh <dataset_id>
#
# Reads dataset metadata from audit3_inputs.tsv to determine chemistry,
# species, FASTQ location. Runs:
#   slot 1: STARsolo default
#   slot 2: STARsolo CR-mimic
#   slot 3: alevin-fry
#   slot 4: kb count
#
# Each slot writes to processed/<dataset_id>/<slot>/.
# Timing + RSS captured to processed/<dataset_id>/<slot>/timing.json.
# Skips a slot if its sentinel output file already exists.
set -euo pipefail

DATASET_ID="${1:?Usage: process_one_dataset.sh <dataset_id>}"

AUDIT_ROOT="/mnt/nvme1/omics-audit/phase2/audit3_counting"
ENV="/home/ross/miniforge3/envs/audit3_counting"
export PATH="$ENV/bin:$PATH"
THREADS="${THREADS:-12}"

INVENTORY="$AUDIT_ROOT/inventory/audit3_inputs.tsv"
PROCESSED_ROOT="$AUDIT_ROOT/processed"
FASTQ_ROOT="/mnt/nvme2/audit3_fastqs"
WHITELIST_DIR="/mnt/nvme2/refs/whitelists"

# --- Look up dataset metadata from inventory (safe: no eval) ---
get_col() {
  # get_col <column_name>  →  prints the value of that column for the current $DATASET_ID
  local col="$1"
  awk -F'\t' -v id="$DATASET_ID" -v col="$col" '
    NR==1 {for(i=1;i<=NF;i++) if ($i==col) ci=i; next}
    $1==id {print $ci}
  ' "$INVENTORY"
}
chemistry=$(get_col chemistry)
chemistry_exact=$(get_col chemistry_exact)
species=$(get_col species)
tissue=$(get_col tissue)

if [[ -z "$chemistry" ]]; then
  echo "ERROR: dataset_id '$DATASET_ID' not found in $INVENTORY" >&2
  exit 1
fi

echo "[$(date -Iseconds)] Processing $DATASET_ID"
echo "  chemistry=$chemistry  exact=$chemistry_exact  species=$species  tissue=$tissue"

# --- Chemistry-dependent parameters ---
case "$chemistry" in
  3p_v2)
    CB_LEN=16; UMI_LEN=10
    WL="$WHITELIST_DIR/737K-august-2016.txt"
    SALMON_FLAG="--chromium"
    KB_X="10xv2"
    SOLO_FLAG=""
    ;;
  3p_v3)
    CB_LEN=16; UMI_LEN=12
    WL="$WHITELIST_DIR/3M-february-2018.txt"
    SALMON_FLAG="--chromiumV3"
    KB_X="10xv3"
    SOLO_FLAG=""
    ;;
  5p_v2)
    # 5' v2 uses same barcode + whitelist structure as 3' v2; only orientation differs
    CB_LEN=16; UMI_LEN=10
    WL="$WHITELIST_DIR/737K-august-2016.txt"
    SALMON_FLAG="--chromium"
    KB_X="10xv2"  # kb's 5p_v2 syntax: 10xv2 with read structure spec; we keep simple
    SOLO_FLAG=""
    ;;
  *)
    echo "ERROR: unsupported chemistry '$chemistry'" >&2
    exit 1
    ;;
esac

# --- Species-dependent references ---
if [[ "$species" == "mouse" ]]; then
  STAR_IDX="/mnt/nvme2/refs/mm39/star_index_sjdb150"
  SALMON_IDX="/mnt/nvme2/refs/mm39/salmon_splici_index"
  KB_IDX="/mnt/nvme2/refs/mm39/kb_index/index.idx"
  KB_T2G="/mnt/nvme2/refs/mm39/kb_index/t2g.txt"
  T2G_3COL="/mnt/nvme2/refs/mm39/salmon_splici_index/t2g_3col.tsv"
else
  STAR_IDX="/mnt/nvme2/refs/GRCh38_gencode45/star_index_sjdb90"
  SALMON_IDX="/mnt/nvme2/refs/GRCh38_gencode45/salmon_splici_index"
  KB_IDX="/mnt/nvme2/refs/GRCh38_gencode45/kb_index/index.idx"
  KB_T2G="/mnt/nvme2/refs/GRCh38_gencode45/kb_index/t2g.txt"
  T2G_3COL="/mnt/nvme2/refs/GRCh38_gencode45/salmon_splici_index/t2g_3col.tsv"
fi

# --- Locate FASTQ files ---
# Convention: $FASTQ_ROOT/<dataset_id_short>/<fastqs_dir>/
# Find R1 + R2 files. Concatenate multi-lane via comma-separated paths for STAR.
case "$DATASET_ID" in
  10x_pbmc_1k_v3)             FQ_DIR="$FASTQ_ROOT/pbmc_1k_v3/pbmc_1k_v3_fastqs" ;;
  10x_t_3k_v2)                FQ_DIR="$FASTQ_ROOT/t_3k_v2/t_3k_fastqs" ;;
  10x_pbmc_4k_v2)             FQ_DIR="$FASTQ_ROOT/pbmc_4k/pbmc4k_fastqs" ;;
  10x_pbmc_5k_v3.1)           FQ_DIR="$FASTQ_ROOT/pbmc_5k_v3_1/Chromium_3p_GEX_Human_PBMC_fastqs" ;;
  10x_pbmc_10k_v3.1)          FQ_DIR="$FASTQ_ROOT/pbmc_10k_v3_1/10k_PBMC_3p_nextgem_Chromium_Controller_fastqs" ;;
  10x_neuron_1k_v3)           FQ_DIR="$FASTQ_ROOT/neuron_1k_v3/neuron_1k_v3_fastqs" ;;
  tabula_muris_liver_droplet) FQ_DIR="$FASTQ_ROOT/tabula_muris_liver" ;;
  tabula_muris_heart_droplet) FQ_DIR="$FASTQ_ROOT/tabula_muris_heart" ;;
  *) echo "ERROR: no FQ_DIR mapping for $DATASET_ID" >&2; exit 1 ;;
esac

# Tabula Muris SRRs land at the root of $FASTQ_ROOT/tabula_muris_<tissue>/
# rather than in a fastqs/ subdir. Handle both.
if [[ ! -d "$FQ_DIR" ]]; then
  parent="$(dirname "$FQ_DIR")"
  if [[ -d "$parent" ]]; then FQ_DIR="$parent"; fi
fi

# Use nullglob to make non-matching globs expand to empty (not error under pipefail)
shopt -s nullglob
R1_FILES=( "$FQ_DIR"/*_R1_*.fastq.gz "$FQ_DIR"/*_1.fastq.gz )
R2_FILES=( "$FQ_DIR"/*_R2_*.fastq.gz "$FQ_DIR"/*_2.fastq.gz )
shopt -u nullglob
# Dedupe + comma-join
R1_LIST=$(printf "%s\n" "${R1_FILES[@]}" | sort -u | grep -v '^$' | tr '\n' ',' | sed 's/,$//')
R2_LIST=$(printf "%s\n" "${R2_FILES[@]}" | sort -u | grep -v '^$' | tr '\n' ',' | sed 's/,$//')

if [[ -z "$R1_LIST" || -z "$R2_LIST" ]]; then
  echo "ERROR: no FASTQs found in $FQ_DIR" >&2; exit 1
fi
echo "  R1 files: $(echo $R1_LIST | tr ',' '\n' | wc -l)"
echo "  R2 files: $(echo $R2_LIST | tr ',' '\n' | wc -l)"

OUTROOT="$PROCESSED_ROOT/$DATASET_ID"
mkdir -p "$OUTROOT"

# --- Helper: timed run + RSS peak via /usr/bin/time ---
run_timed() {
  local slot="$1"; shift
  local out_dir="$OUTROOT/$slot"
  local time_log="$out_dir/timing_raw.txt"
  mkdir -p "$out_dir"
  /usr/bin/time -v -o "$time_log" "$@" 2>&1 | tee "$out_dir/run.log" | tail -3
  # Parse timing.json
  python3 - "$time_log" "$out_dir/timing.json" "$slot" "$DATASET_ID" <<'PY'
import sys, json, re, os
inp, out, slot, ds = sys.argv[1:5]
text = open(inp).read()
def find(pat, default=None):
    m = re.search(pat, text)
    return m.group(1) if m else default
res = {
    "dataset_id": ds, "slot": slot,
    "user_seconds": float(find(r'User time \(seconds\): ([0-9.]+)', 0)),
    "system_seconds": float(find(r'System time \(seconds\): ([0-9.]+)', 0)),
    "elapsed_str": find(r'Elapsed \(wall clock\) time .*: (.+)', ''),
    "max_rss_kb": int(find(r'Maximum resident set size .*: ([0-9]+)', 0)),
    "pct_cpu": find(r'Percent of CPU this job got: (.+)', ''),
}
res["max_rss_mb"] = res["max_rss_kb"] // 1024 if res["max_rss_kb"] else 0
res["max_rss_gb"] = round(res["max_rss_kb"] / 1024 / 1024, 2) if res["max_rss_kb"] else 0
json.dump(res, open(out, "w"), indent=2)
print(json.dumps({k:v for k,v in res.items() if k != "max_rss_kb"}))
PY
}

# ============================================================
# Slot 1: STARsolo default
# ============================================================
SLOT1="$OUTROOT/star_default"
if [[ ! -f "$SLOT1/Solo.out/Gene/raw/matrix.mtx" ]]; then
  echo "[$(date -Iseconds)] [SLOT 1] STARsolo default — $DATASET_ID"
  mkdir -p "$SLOT1"
  run_timed star_default STAR \
    --runMode alignReads \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WL" \
    --soloFeatures Gene \
    --genomeDir "$STAR_IDX" \
    --readFilesIn "$R2_LIST" "$R1_LIST" \
    --readFilesCommand zcat \
    --soloCBstart 1 --soloCBlen $CB_LEN \
    --soloUMIstart $((CB_LEN+1)) --soloUMIlen $UMI_LEN \
    --outSAMtype None \
    --runThreadN "$THREADS" \
    --outFileNamePrefix "$SLOT1/"
else
  echo "[$(date -Iseconds)] [SLOT 1] $DATASET_ID — already done, skipping"
fi

# ============================================================
# Slot 2: STARsolo CR-mimic
# ============================================================
SLOT2="$OUTROOT/star_cr_mimic"
if [[ ! -f "$SLOT2/Solo.out/Gene/raw/matrix.mtx" ]]; then
  echo "[$(date -Iseconds)] [SLOT 2] STARsolo CR-mimic — $DATASET_ID"
  mkdir -p "$SLOT2"
  run_timed star_cr_mimic STAR \
    --runMode alignReads \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "$WL" \
    --soloFeatures Gene \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --clipAdapterType CellRanger4 \
    --outFilterScoreMin 30 \
    --genomeDir "$STAR_IDX" \
    --readFilesIn "$R2_LIST" "$R1_LIST" \
    --readFilesCommand zcat \
    --soloCBstart 1 --soloCBlen $CB_LEN \
    --soloUMIstart $((CB_LEN+1)) --soloUMIlen $UMI_LEN \
    --outSAMtype None \
    --runThreadN "$THREADS" \
    --outFileNamePrefix "$SLOT2/"
else
  echo "[$(date -Iseconds)] [SLOT 2] $DATASET_ID — already done, skipping"
fi

# ============================================================
# Slot 3: alevin-fry
# ============================================================
SLOT3="$OUTROOT/alevin_fry"
if [[ ! -f "$SLOT3/quant/alevin/quants_mat.mtx" ]]; then
  echo "[$(date -Iseconds)] [SLOT 3] alevin-fry — $DATASET_ID"
  mkdir -p "$SLOT3"
  # First: salmon alevin to produce RAD file
  run_timed alevin_fry bash -c "
    salmon alevin -i $SALMON_IDX -l ISR $SALMON_FLAG \
      -1 $(echo "$R1_LIST" | tr ',' ' ') \
      -2 $(echo "$R2_LIST" | tr ',' ' ') \
      -p $THREADS --sketch -o $SLOT3 && \
    alevin-fry generate-permit-list --input $SLOT3 --expected-ori fw \
      --unfiltered-pl $WL --output-dir $SLOT3/permit && \
    alevin-fry collate -t $THREADS --input-dir $SLOT3/permit --rad-dir $SLOT3 && \
    alevin-fry quant -t $THREADS --input-dir $SLOT3/permit \
      -m $T2G_3COL -r cr-like --use-mtx -o $SLOT3/quant
  "
else
  echo "[$(date -Iseconds)] [SLOT 3] $DATASET_ID — already done, skipping"
fi

# ============================================================
# Slot 4: kb count
# ============================================================
SLOT4="$OUTROOT/kb_count"
if [[ ! -f "$SLOT4/counts_unfiltered/cells_x_genes.mtx" ]]; then
  echo "[$(date -Iseconds)] [SLOT 4] kb count — $DATASET_ID"
  mkdir -p "$SLOT4"
  R1_SPACE=$(echo "$R1_LIST" | tr ',' ' ')
  R2_SPACE=$(echo "$R2_LIST" | tr ',' ' ')
  PAIRS=""
  R1_ARR=($R1_SPACE); R2_ARR=($R2_SPACE)
  for i in "${!R1_ARR[@]}"; do
    PAIRS="$PAIRS ${R1_ARR[$i]} ${R2_ARR[$i]}"
  done
  run_timed kb_count bash -c "
    kb count -i $KB_IDX -g $KB_T2G -x $KB_X \
      -o $SLOT4 --filter bustools -t $THREADS $PAIRS
  "
else
  echo "[$(date -Iseconds)] [SLOT 4] $DATASET_ID — already done, skipping"
fi

# --- Summary ---
echo "[$(date -Iseconds)] $DATASET_ID DONE"
for slot in star_default star_cr_mimic alevin_fry kb_count; do
  case $slot in
    star_default|star_cr_mimic) f="$OUTROOT/$slot/Solo.out/Gene/raw/matrix.mtx" ;;
    alevin_fry) f="$OUTROOT/$slot/quant/alevin/quants_mat.mtx" ;;
    kb_count) f="$OUTROOT/$slot/counts_unfiltered/cells_x_genes.mtx" ;;
  esac
  if [[ -f "$f" ]]; then
    sz=$(du -h "$f" 2>/dev/null | cut -f1)
    t=$(jq -r '.elapsed_str // "?"' "$OUTROOT/$slot/timing.json" 2>/dev/null || echo "?")
    rss=$(jq -r '.max_rss_gb // "?"' "$OUTROOT/$slot/timing.json" 2>/dev/null || echo "?")
    printf "  %-15s OK  matrix=%s  elapsed=%s  rss_peak=%sG\n" "$slot" "$sz" "$t" "$rss"
  else
    printf "  %-15s MISSING (%s)\n" "$slot" "$f"
  fi
done
