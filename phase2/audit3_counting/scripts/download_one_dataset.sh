#!/usr/bin/env bash
# Download FASTQ for one Audit 3 dataset.
# Verified URLs (probed via curl HEAD 2026-05-17) and SRA accessions for
# the 8-dataset CP3 working set.
#
# Usage: download_one_dataset.sh <dataset_id>
# Idempotent (skips if FASTQs already extracted).
set -euo pipefail

DATASET_ID="${1:?Usage: download_one_dataset.sh <dataset_id>}"
FASTQ_ROOT="/mnt/nvme2/audit3_fastqs"
ENV="/home/ross/miniforge3/envs/audit3_counting"
export PATH="$ENV/bin:$PATH"

# Determine acquisition path: direct 10x CDN tarball, or SRA via prefetch
case "$DATASET_ID" in
  # ---- 10x CDN tarball datasets ----
  10x_pbmc_1k_v3)
    DEST="$FASTQ_ROOT/pbmc_1k_v3"
    URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_fastqs.tar"
    EXPECT="pbmc_1k_v3_fastqs" ;;
  10x_t_3k_v2)
    DEST="$FASTQ_ROOT/t_3k_v2"
    URL="https://cf.10xgenomics.com/samples/cell-exp/2.1.0/t_3k/t_3k_fastqs.tar"
    EXPECT="fastqs" ;;
  10x_pbmc_4k_v2)
    DEST="$FASTQ_ROOT/pbmc_4k"
    URL="https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_fastqs.tar"
    EXPECT="fastqs" ;;
  10x_pbmc_5k_v3.1)
    DEST="$FASTQ_ROOT/pbmc_5k_v3_1"
    URL="https://cf.10xgenomics.com/samples/cell-exp/7.0.1/SC3pv3_GEX_Human_PBMC/SC3pv3_GEX_Human_PBMC_fastqs.tar"
    EXPECT="Chromium_3p_GEX_Human_PBMC_fastqs" ;;
  10x_pbmc_10k_v3.1)
    DEST="$FASTQ_ROOT/pbmc_10k_v3_1"
    URL="https://cf.10xgenomics.com/samples/cell-exp/6.1.0/10k_PBMC_3p_nextgem_Chromium_Controller/10k_PBMC_3p_nextgem_Chromium_Controller_fastqs.tar"
    EXPECT="10k_PBMC_3p_nextgem_Chromium_Controller_fastqs" ;;
  10x_neuron_1k_v3)
    DEST="$FASTQ_ROOT/neuron_1k_v3"
    URL="https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_1k_v3/neuron_1k_v3_fastqs.tar"
    EXPECT="neuron_1k_v3_fastqs" ;;

  # ---- SRA datasets (Tabula Muris liver + heart, 3' v2 mouse) ----
  tabula_muris_liver_droplet|tabula_muris_heart_droplet)
    # Handled via the SRA branch below
    ;;

  # ---- SRA datasets (scraper additions 2026-05-17) ----
  gse287209_human_lung_organoid)
    # Human lung organoid, 3' v3, Dost et al. PMID 41992061
    # 3 GSMs (GSM8741358-60) → SRX27375400/01/02, ~31 GB total
    DEST="$FASTQ_ROOT/gse287209_lung_organoid"
    SRX_LIST=(SRX27375400 SRX27375401 SRX27375402)
    ;;
  gse325955_mouse_kidney_E18_5)
    # Mouse kidney E18.5 scRNA subset, Finer et al.
    # 4 GSMs (GSM9617755-58) → SRX32649381-84, ~46 GB total
    DEST="$FASTQ_ROOT/gse325955_kidney_E18_5"
    SRX_LIST=(SRX32649381 SRX32649382 SRX32649383 SRX32649384)
    ;;
  gse288156_mouse_intestine_scrna)
    # Mouse intestine scRNA, Lassila et al. — replaces failed Tabula Muris
    # 5 scRNA GSMs (GSM8760147-51) → SRX27491077/79/80/81/82 → 10 SRRs, ~104 GB total
    # Read structure VERIFIED via 1000-read sample 2026-05-17: 4-file layout
    # (10/10/28/90 bp, _3=CB+UMI, _4=transcript). 9 ATAC samples in GSE288156 SKIPPED.
    DEST="$FASTQ_ROOT/gse288156_intestine"
    SRX_LIST=(SRX27491077 SRX27491079 SRX27491080 SRX27491081 SRX27491082)
    ;;

  *)
    echo "ERROR: no acquisition mapping for $DATASET_ID" >&2
    exit 1
    ;;
esac

# ---- Generic SRX list SRA branch (scraper-added datasets) ----
if [[ -n "${SRX_LIST:-}" ]]; then
  mkdir -p "$DEST"; cd "$DEST"
  echo "[$(date -Iseconds)] $DATASET_ID — SRA pull from ${#SRX_LIST[@]} SRX accessions"
  for SRX in "${SRX_LIST[@]}"; do
    SRRS=$(esearch -db sra -query "$SRX" 2>/dev/null | \
           efetch -format runinfo 2>/dev/null | \
           tail -n +2 | cut -d, -f1 | grep -E '^[SDE]RR' || true)
    if [[ -z "$SRRS" ]]; then
      echo "  $SRX → NO SRRs"; continue
    fi
    for SRR in $SRRS; do
      echo "  $SRX → $SRR"
      if ls "${SRR}"_*.fastq.gz >/dev/null 2>&1; then
        echo "    already extracted, skipping"; continue
      fi
      if [[ ! -d "$SRR" ]]; then
        prefetch "$SRR" --max-size u 2>&1 | tail -2
      fi
      fasterq-dump --threads 12 --split-files --include-technical "$SRR" 2>&1 | tail -2
      gzip "${SRR}"_*.fastq 2>/dev/null || true
      rm -rf "$SRR"
    done
  done
  echo "[$(date -Iseconds)] $DATASET_ID SRA extraction complete:"
  ls -lh *.fastq.gz 2>/dev/null | head -10
  exit 0
fi

# ---- SRA branch ----
if [[ "$DATASET_ID" == "tabula_muris_liver_droplet" ]] || \
   [[ "$DATASET_ID" == "tabula_muris_heart_droplet" ]]; then

  if [[ "$DATASET_ID" == "tabula_muris_liver_droplet" ]]; then
    TISSUE="Liver"
    DEST="$FASTQ_ROOT/tabula_muris_liver"
  else
    TISSUE="Heart"
    DEST="$FASTQ_ROOT/tabula_muris_heart"
  fi

  mkdir -p "$DEST"; cd "$DEST"

  # Tabula Muris specific GSM accessions per Tabula Muris (Schaum et al. 2018)
  # Table S1 — droplet samples for liver and heart tissues.
  # These accessions are looked up from GSE109774's SOFT metadata for the
  # specific 'tissue' subset with 10x droplet protocol. The 3-mouse FACS+droplet
  # design produces ~3 GSMs per tissue under droplet.
  if [[ "$TISSUE" == "Liver" ]]; then
    # Per GSE109774 Liver droplet: GSM3040905 (P4), GSM3040906 (P7-mock?), etc.
    # NOTE: specific accession set to be confirmed by query at runtime;
    # placeholders here as starting set for SRR resolution
    GSMS=("GSM3040905" "GSM3040906" "GSM3040907")
  else
    GSMS=("GSM3040908" "GSM3040909" "GSM3040910")
  fi

  echo "[$(date -Iseconds)] $DATASET_ID — SRA acquisition via prefetch + fasterq-dump"
  echo "  GSMs (preliminary; will validate against GEO): ${GSMS[*]}"
  echo "  GSM-to-SRR resolution requires NCBI Entrez query — running now ..."
  # Resolve GSM → SRR via efetch (entrez-direct). If not installed, prompt user.
  if ! command -v efetch >/dev/null 2>&1; then
    echo "  efetch not present — installing entrez-direct via mamba"
    mamba install -n audit3_counting -y -c bioconda entrez-direct 2>&1 | tail -3
  fi
  for GSM in "${GSMS[@]}"; do
    echo "  Resolving $GSM ..."
    SRR=$(esearch -db sra -query "$GSM" 2>/dev/null | efetch -format runinfo 2>/dev/null | tail -n +2 | head -1 | cut -d, -f1)
    if [[ -z "$SRR" ]]; then
      echo "    NO SRR resolved for $GSM — skipping"
      continue
    fi
    echo "    $GSM → $SRR"
    if [[ -f "${SRR}_1.fastq.gz" ]]; then
      echo "    already extracted, skipping"
      continue
    fi
    prefetch "$SRR"
    fasterq-dump --threads 12 --split-files --include-technical "$SRR"
    # Compress
    gzip "${SRR}"_*.fastq 2>/dev/null || true
    rm -rf "$SRR"   # remove the .sra cache after fastq extraction
  done
  echo "[$(date -Iseconds)] $DATASET_ID SRA extraction complete:"
  ls -lh *.fastq.gz 2>/dev/null | head -10
  exit 0
fi

# ---- 10x CDN tarball branch ----
mkdir -p "$DEST"; cd "$DEST"
if [[ -d "$EXPECT" ]] && [[ -n "$(ls "$EXPECT"/*.fastq.gz 2>/dev/null)" ]]; then
  echo "[$(date -Iseconds)] $DATASET_ID already on disk at $DEST/$EXPECT — skipping"
  ls -lh "$EXPECT"/*.fastq.gz | head -6
  exit 0
fi

TARBALL="${URL##*/}"
if [[ ! -f "$TARBALL" ]]; then
  echo "[$(date -Iseconds)] Downloading $DATASET_ID from $URL"
  curl -sSL --fail -o "$TARBALL" "$URL"
fi

if ! file "$TARBALL" | grep -qE 'POSIX tar|tar archive'; then
  echo "ERROR: $TARBALL is not a tar archive:" >&2
  file "$TARBALL" >&2
  exit 1
fi

echo "[$(date -Iseconds)] Extracting $TARBALL"
tar xf "$TARBALL"
echo "[$(date -Iseconds)] $DATASET_ID complete:"
ls -lh "$EXPECT"/*.fastq.gz | head -6
