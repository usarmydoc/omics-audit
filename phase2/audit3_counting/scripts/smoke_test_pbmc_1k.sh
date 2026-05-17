#!/usr/bin/env bash
# CP2 smoke test — run all 4 tool configurations on a HEAD subset of
# pbmc_1k_v3 to confirm each tool reads its index and produces output
# end-to-end. Not a benchmark; just sanity that the pipelines fire.
#
# Strategy:
#   1. Use head subset of the FASTQs (first 200k reads) for speed
#   2. Run each tool config end-to-end
#   3. Check that a count matrix file exists at the expected output path
#   4. Report per-tool timing + output sanity
set -euo pipefail

ENV="/home/ross/miniforge3/envs/audit3_counting"
export PATH="$ENV/bin:$PATH"
THREADS="${THREADS:-12}"

FASTQ_DIR="/mnt/nvme2/audit3_fastqs/pbmc_1k_v3/pbmc_1k_v3_fastqs"
REF_DIR="/mnt/nvme2/refs/GRCh38_gencode45"
WORK="/mnt/nvme2/audit3_fastqs/pbmc_1k_v3/_smoke_test"
mkdir -p "$WORK"
cd "$WORK"

# Build head subset (200k reads per R1/R2/I1) — fast enough for smoke test
SUB_R1="$WORK/sub_R1.fastq.gz"
SUB_R2="$WORK/sub_R2.fastq.gz"
N_READS=200000
if [[ ! -f "$SUB_R2" ]]; then
  echo "[$(date -Iseconds)] Building $N_READS-read head subset"
  # head closes the pipe early → zcat gets SIGPIPE; locally relax pipefail for the subset build
  set +o pipefail
  zcat "$FASTQ_DIR"/pbmc_1k_v3_S1_L001_R1_001.fastq.gz "$FASTQ_DIR"/pbmc_1k_v3_S1_L002_R1_001.fastq.gz | head -n $((N_READS * 4)) | gzip > "$SUB_R1"
  zcat "$FASTQ_DIR"/pbmc_1k_v3_S1_L001_R2_001.fastq.gz "$FASTQ_DIR"/pbmc_1k_v3_S1_L002_R2_001.fastq.gz | head -n $((N_READS * 4)) | gzip > "$SUB_R2"
  set -o pipefail
fi
echo "Subset sizes:"
ls -lh "$SUB_R1" "$SUB_R2"

# 10x v3 whitelist (3M-february-2018.txt.gz packaged with CellRanger;
# also distributed by STARsolo / pyroe / kb-python). Pull from pyroe data.
WL="/mnt/nvme2/refs/whitelists/3M-february-2018.txt.gz"
if [[ ! -f "$WL" ]]; then
  echo "[$(date -Iseconds)] Downloading 10x v3 whitelist"
  mkdir -p "$(dirname "$WL")"
  # Working source: teichlab/scg_lib_structs mirror (10x's GitHub repo is gated)
  curl -sSL -o "$WL" https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/3M-february-2018.txt.gz
  gunzip -t "$WL"
fi
# Tools want an uncompressed plain-text whitelist; un-gz it once
if [[ ! -f "${WL%.gz}" ]]; then
  gunzip -k "$WL"
fi

# ============================================================
# Slot 1: STARsolo default
# ============================================================
SLOT1="$WORK/star_default"
if [[ ! -d "$SLOT1/Solo.out" ]]; then
  echo "[$(date -Iseconds)] [SLOT 1] STARsolo default"
  mkdir -p "$SLOT1"
  cd "$SLOT1"
  STAR \
    --runMode alignReads \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WL%.gz}" \
    --soloFeatures Gene \
    --genomeDir "$REF_DIR/star_index_sjdb90" \
    --readFilesIn "$SUB_R2" "$SUB_R1" \
    --readFilesCommand zcat \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --outSAMtype None \
    --runThreadN "$THREADS" \
    --outFileNamePrefix "$SLOT1/" \
    2>&1 | tail -5
  cd "$WORK"
fi

# ============================================================
# Slot 2: STARsolo CR-mimic (CellRanger-compatible flags)
# ============================================================
SLOT2="$WORK/star_cr_mimic"
if [[ ! -d "$SLOT2/Solo.out" ]]; then
  echo "[$(date -Iseconds)] [SLOT 2] STARsolo CR-mimic"
  mkdir -p "$SLOT2"
  cd "$SLOT2"
  STAR \
    --runMode alignReads \
    --soloType CB_UMI_Simple \
    --soloCBwhitelist "${WL%.gz}" \
    --soloFeatures Gene \
    --soloUMIfiltering MultiGeneUMI_CR \
    --soloUMIdedup 1MM_CR \
    --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
    --clipAdapterType CellRanger4 \
    --outFilterScoreMin 30 \
    --genomeDir "$REF_DIR/star_index_sjdb90" \
    --readFilesIn "$SUB_R2" "$SUB_R1" \
    --readFilesCommand zcat \
    --soloCBstart 1 --soloCBlen 16 \
    --soloUMIstart 17 --soloUMIlen 12 \
    --outSAMtype None \
    --runThreadN "$THREADS" \
    --outFileNamePrefix "$SLOT2/" \
    2>&1 | tail -5
  cd "$WORK"
fi

# ============================================================
# Slot 3: alevin-fry
# ============================================================
SLOT3="$WORK/alevin_fry"
if [[ ! -d "$SLOT3/quant/alevin" ]]; then
  echo "[$(date -Iseconds)] [SLOT 3] alevin-fry (salmon alevin + alevin-fry)"
  mkdir -p "$SLOT3"
  salmon alevin \
    -i "$REF_DIR/salmon_splici_index" \
    -l ISR \
    --chromiumV3 \
    -1 "$SUB_R1" \
    -2 "$SUB_R2" \
    -p "$THREADS" --sketch \
    -o "$SLOT3" \
    2>&1 | tail -3
  alevin-fry generate-permit-list \
    --input "$SLOT3" \
    --expected-ori fw \
    --unfiltered-pl "${WL%.gz}" \
    --output-dir "$SLOT3/permit" \
    2>&1 | tail -3
  alevin-fry collate -t "$THREADS" \
    --input-dir "$SLOT3/permit" --rad-dir "$SLOT3" \
    2>&1 | tail -3
  alevin-fry quant -t "$THREADS" \
    --input-dir "$SLOT3/permit" \
    -m "$REF_DIR/salmon_splici_index/t2g_3col.tsv" \
    -r cr-like \
    --use-mtx \
    -o "$SLOT3/quant" \
    2>&1 | tail -3
fi

# ============================================================
# Slot 4: kb count (kallisto + bustools)
# ============================================================
SLOT4="$WORK/kb_count"
if [[ ! -d "$SLOT4/counts_unfiltered" ]]; then
  echo "[$(date -Iseconds)] [SLOT 4] kb count"
  mkdir -p "$SLOT4"
  kb count \
    -i "$REF_DIR/kb_index/index.idx" \
    -g "$REF_DIR/kb_index/t2g.txt" \
    -x 10xv3 \
    -o "$SLOT4" \
    --filter bustools \
    -t "$THREADS" \
    "$SUB_R1" "$SUB_R2" \
    2>&1 | tail -3
fi

# ============================================================
# Report
# ============================================================
echo ""
echo "[$(date -Iseconds)] === smoke test summary ==="
for slot in star_default star_cr_mimic alevin_fry kb_count; do
  printf "%-15s " "$slot"
  case $slot in
    star_default|star_cr_mimic)
      f="$WORK/$slot/Solo.out/Gene/raw/matrix.mtx"
      ;;
    alevin_fry)
      f="$WORK/$slot/quant/alevin/quants_mat.gz"
      ;;
    kb_count)
      f="$WORK/$slot/counts_unfiltered/cells_x_genes.mtx"
      ;;
  esac
  if [[ -f "$f" ]]; then
    sz=$(du -h "$f" 2>/dev/null | cut -f1)
    echo "OK ($f $sz)"
  else
    echo "MISSING ($f)"
  fi
done
echo ""
echo "[$(date -Iseconds)] smoke test DONE"
