"""
B2 data access: download GTEx V10 raw counts and TPM from Google Cloud bucket.

Downloads per-tissue GCT files (gene_reads + gene_tpm) for 20 tissues.
"""

import gzip
import io
from pathlib import Path

import pandas as pd
import requests
from tqdm import tqdm

BUCKET_BASE = "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq"
COUNTS_DIR = f"{BUCKET_BASE}/counts-by-tissue"
TPM_DIR = f"{BUCKET_BASE}/tpms-by-tissue"

# 20 tissues prioritizing diversity across organ systems
TISSUES = {
    # Brain
    "brain_cortex": "Brain - Cortex",
    # Heart
    "heart_left_ventricle": "Heart - Left Ventricle",
    # Liver
    "liver": "Liver",
    # Kidney
    "kidney_cortex": "Kidney - Cortex",
    # Lung
    "lung": "Lung",
    # Muscle
    "muscle_skeletal": "Muscle - Skeletal",
    # Skin
    "skin_sun_exposed_lower_leg": "Skin - Sun Exposed (Lower leg)",
    # Gut
    "colon_transverse": "Colon - Transverse",
    "small_intestine_terminal_ileum": "Small Intestine - Terminal Ileum",
    "stomach": "Stomach",
    # Blood
    "whole_blood": "Whole Blood",
    # Reproductive
    "testis": "Testis",
    "ovary": "Ovary",
    "uterus": "Uterus",
    "prostate": "Prostate",
    # Endocrine
    "thyroid": "Thyroid",
    "pituitary": "Pituitary",
    "adrenal_gland": "Adrenal Gland",
    # Vascular
    "artery_aorta": "Artery - Aorta",
    # Nerve
    "nerve_tibial": "Nerve - Tibial",
}


def _parse_gct(content: str) -> pd.DataFrame:
    """Parse GCT v1.2 format into a DataFrame (genes × samples)."""
    lines = content.split("\n")
    # Line 0: #1.2, Line 1: nrows ncols, Line 2+: header + data
    header = lines[2].split("\t")
    sample_ids = header[2:]

    data_lines = []
    gene_ids = []
    gene_descs = []
    for line in lines[3:]:
        if not line.strip():
            continue
        parts = line.split("\t")
        gene_ids.append(parts[0])
        gene_descs.append(parts[1])
        data_lines.append(parts[2:])

    df = pd.DataFrame(data_lines, index=gene_ids, columns=sample_ids)
    df.index.name = "gene_id"
    return df, gene_descs


def download_tissue(tissue_key: str, cache_dir: Path,
                    max_samples: int = 100) -> dict:
    """Download counts and TPM for one tissue.

    Returns dict with 'counts' (int DataFrame), 'tpm' (float DataFrame),
    'tissue_key', 'tissue_name', 'n_samples', 'gene_descriptions'.
    """
    counts_cache = cache_dir / f"{tissue_key}_counts.parquet"
    tpm_cache = cache_dir / f"{tissue_key}_tpm.parquet"

    if counts_cache.exists() and tpm_cache.exists():
        counts = pd.read_parquet(counts_cache)
        tpm = pd.read_parquet(tpm_cache)
        return {
            "counts": counts,
            "tpm": tpm,
            "tissue_key": tissue_key,
            "tissue_name": TISSUES[tissue_key],
            "n_samples": counts.shape[1],
        }

    cache_dir.mkdir(parents=True, exist_ok=True)

    # Download counts
    counts_url = f"{COUNTS_DIR}/gene_reads_v10_{tissue_key}.gct.gz"
    print(f"    Downloading counts: {tissue_key}...")
    r = requests.get(counts_url, timeout=120)
    r.raise_for_status()
    content = gzip.decompress(r.content).decode("utf-8")
    counts_df, gene_descs = _parse_gct(content)

    # Download TPM
    tpm_url = f"{TPM_DIR}/gene_tpm_v10_{tissue_key}.gct.gz"
    print(f"    Downloading TPM: {tissue_key}...")
    r2 = requests.get(tpm_url, timeout=120)
    r2.raise_for_status()
    content2 = gzip.decompress(r2.content).decode("utf-8")
    tpm_df, _ = _parse_gct(content2)

    # Ensure same samples in both
    common_samples = counts_df.columns.intersection(tpm_df.columns)
    if max_samples and len(common_samples) > max_samples:
        common_samples = common_samples[:max_samples]

    counts_df = counts_df[common_samples].apply(pd.to_numeric, errors="coerce").fillna(0).astype(int)
    tpm_df = tpm_df[common_samples].apply(pd.to_numeric, errors="coerce").fillna(0.0)

    # Cache
    counts_df.to_parquet(counts_cache)
    tpm_df.to_parquet(tpm_cache)

    return {
        "counts": counts_df,
        "tpm": tpm_df,
        "tissue_key": tissue_key,
        "tissue_name": TISSUES[tissue_key],
        "n_samples": counts_df.shape[1],
    }
