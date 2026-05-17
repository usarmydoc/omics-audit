#!/usr/bin/env bash
# Resolve Tabula Muris (GSE109774) liver + heart droplet GSMs to SRR accessions,
# then prefetch + fasterq-dump them. Replaces the placeholder GSMs in the
# original download script with real ones via NCBI Entrez.
#
# Usage: resolve_tabula_muris.sh <liver|heart>
set -euo pipefail

TISSUE_QUERY="${1:?Usage: resolve_tabula_muris.sh <liver|heart>}"
case "$TISSUE_QUERY" in
  liver) TISSUE_PRETTY="Liver"; DEST="/mnt/nvme2/audit3_fastqs/tabula_muris_liver" ;;
  heart) TISSUE_PRETTY="Heart"; DEST="/mnt/nvme2/audit3_fastqs/tabula_muris_heart" ;;
  *) echo "ERROR: tissue must be 'liver' or 'heart'" >&2; exit 1 ;;
esac

ENV="/home/ross/miniforge3/envs/audit3_counting"
export PATH="$ENV/bin:$PATH"
mkdir -p "$DEST"; cd "$DEST"

echo "[$(date -Iseconds)] Resolving GSE109774 → SRA runinfo (Tabula Muris)"
RUNINFO="$DEST/sra_runinfo.csv"
if [[ ! -s "$RUNINFO" ]]; then
  # Fetch SRA runinfo for the entire GSE; will filter to droplet+tissue downstream
  esearch -db sra -query "GSE109774[ACCN]" 2>/dev/null | \
    efetch -format runinfo 2>/dev/null > "$RUNINFO"
  echo "  rows: $(wc -l < "$RUNINFO")"
fi

# RA runinfo columns include LibraryName, SampleName, source_name etc.
# Tabula Muris convention: LibraryName like "10X_P4_3" where P4 = mouse/plate ID.
# Tissue + method come from the GEO sample title.
# Pull each unique SRA Sample (BioSample) and its title via efetch -format docsum.
echo "[$(date -Iseconds)] Resolving GSM titles to identify ${TISSUE_PRETTY} droplet samples"

# Get title per BioSample via a second Entrez query
# (BioSamples link to GSMs which have descriptive titles like "Liver: 10X Droplet 3M-7")
TITLES_FILE="$DEST/gsm_titles.tsv"
if [[ ! -s "$TITLES_FILE" ]]; then
  # Pull the GSE-level sample table from GEO via the gds database
  esearch -db gds -query "GSE109774[ACCN]" 2>/dev/null | \
    efetch -format docsum 2>/dev/null > "$DEST/gse_docsum.xml"

  # Per-sample IDs in the GSE soft format
  esearch -db gds -query "GSE109774[ACCN]" 2>/dev/null | \
    elink -target gds -name gds_gds 2>/dev/null | \
    efetch -format docsum 2>/dev/null | \
    xtract -pattern DocumentSummary -element Accession,title 2>/dev/null > "$TITLES_FILE" || true

  if [[ ! -s "$TITLES_FILE" ]]; then
    echo "  Falling back to NCBI E-utilities REST API for sample lookup"
    # Direct call: list GSM accessions + titles from GSE matrix metadata
    GSM_LINK=$(esearch -db gds -query "GSE109774[ACCN]" 2>/dev/null | \
               elink -target gds -name gds_gds 2>/dev/null | \
               efetch -format uid 2>/dev/null)
    : # leaving empty fallback; we'll handle below if titles missing
  fi
  echo "  resolved $(wc -l < "$TITLES_FILE") sample titles"
fi

# If titles missing, use SAMPLE_TITLES via direct GEO matrix download
if [[ ! -s "$TITLES_FILE" ]]; then
  echo "  Pulling GSE matrix series_matrix.txt.gz for title parsing"
  curl -sSL -o "$DEST/series_matrix.txt.gz" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE109nnn/GSE109774/matrix/GSE109774_series_matrix.txt.gz"
  zcat "$DEST/series_matrix.txt.gz" | \
    awk -F'\t' '
      /^!Sample_geo_accession/ { for(i=2;i<=NF;i++){ gsub(/"/,"",$i); acc[i]=$i } }
      /^!Sample_title/         { for(i=2;i<=NF;i++){ gsub(/"/,"",$i); print acc[i]"\t"$i } }
    ' > "$TITLES_FILE"
  echo "  resolved $(wc -l < "$TITLES_FILE") sample titles from matrix"
fi

# Filter: tissue match + droplet protocol
MATCHED_GSMS=$(awk -F'\t' -v tissue="$TISSUE_PRETTY" '
  /Droplet/ && tolower($2) ~ tolower(tissue) { print $1 }
' "$TITLES_FILE")

if [[ -z "$MATCHED_GSMS" ]]; then
  echo "ERROR: no $TISSUE_PRETTY droplet GSMs found in TITLES_FILE" >&2
  echo "Sample of titles file:" >&2
  head -10 "$TITLES_FILE" >&2
  exit 1
fi

echo "[$(date -Iseconds)] ${TISSUE_PRETTY} droplet GSMs found:"
echo "$MATCHED_GSMS"

# Resolve each GSM → SRR via runinfo
for GSM in $MATCHED_GSMS; do
  # GSMs map to BioSamples; BioSamples map to SRA runs
  SRRS=$(esearch -db sra -query "$GSM" 2>/dev/null | \
         efetch -format runinfo 2>/dev/null | \
         tail -n +2 | cut -d, -f1 | grep -E '^[SDE]RR' || true)
  if [[ -z "$SRRS" ]]; then
    echo "  $GSM → NO SRR (skipping)"
    continue
  fi
  for SRR in $SRRS; do
    echo "  $GSM → $SRR"
    if ls "${SRR}"_*.fastq.gz >/dev/null 2>&1; then
      echo "    already extracted, skipping"
      continue
    fi
    if [[ ! -d "$SRR" ]]; then
      prefetch "$SRR" --max-size u 2>&1 | tail -2
    fi
    fasterq-dump --threads 12 --split-files --include-technical "$SRR" 2>&1 | tail -2
    gzip "${SRR}"_*.fastq 2>/dev/null || true
    rm -rf "$SRR"
  done
done

echo "[$(date -Iseconds)] $TISSUE_PRETTY droplet acquisition complete:"
ls -lh *.fastq.gz 2>/dev/null | head -10
echo "Total bytes: $(du -sh "$DEST" | cut -f1)"
