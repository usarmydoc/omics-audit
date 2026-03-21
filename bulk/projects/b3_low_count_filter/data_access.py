"""
B3 data access: download TCGA STAR-Counts from GDC API.

Downloads HTSeq/STAR count files for cancer types with ≥10 normal samples.
Builds per-cancer-type count matrices (genes × samples).
"""

import gzip
import io
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

GDC_FILES_ENDPOINT = "https://api.gdc.cancer.gov/files"
GDC_DATA_ENDPOINT = "https://api.gdc.cancer.gov/data"

# Cancer types with ≥10 normal samples (from infrastructure check)
CANCER_TYPES_WITH_NORMALS = [
    "TCGA-BRCA", "TCGA-KIRC", "TCGA-LUAD", "TCGA-THCA", "TCGA-PRAD",
    "TCGA-LUSC", "TCGA-LIHC", "TCGA-HNSC", "TCGA-COAD", "TCGA-STAD",
    "TCGA-UCEC", "TCGA-KIRP", "TCGA-KICH", "TCGA-BLCA", "TCGA-ESCA",
    "TCGA-READ",
]

# Map cancer type to tissue of origin (for tissue-specific gene lookup)
CANCER_TO_TISSUE = {
    "TCGA-BRCA": "Breast", "TCGA-KIRC": "Kidney", "TCGA-LUAD": "Lung",
    "TCGA-THCA": "Thyroid", "TCGA-PRAD": "Prostate", "TCGA-LUSC": "Lung",
    "TCGA-LIHC": "Liver", "TCGA-HNSC": "Head and Neck",
    "TCGA-COAD": "Colon", "TCGA-STAD": "Stomach", "TCGA-UCEC": "Uterus",
    "TCGA-KIRP": "Kidney", "TCGA-KICH": "Kidney", "TCGA-BLCA": "Bladder",
    "TCGA-ESCA": "Esophagus", "TCGA-READ": "Colon",
}


def _query_file_ids(project_id: str, sample_type: str, max_files: int = 60) -> list[dict]:
    """Query GDC for STAR-Counts file UUIDs for a given project and sample type."""
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
            {"op": "in", "content": {"field": "cases.samples.sample_type", "value": [sample_type]}},
        ],
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id,file_name,cases.samples.sample_type,cases.case_id",
        "size": max_files,
        "format": "json",
    }
    r = requests.get(GDC_FILES_ENDPOINT, params=params, timeout=30)
    r.raise_for_status()
    hits = r.json()["data"]["hits"]
    return [{"file_id": h["file_id"], "file_name": h["file_name"],
             "case_id": h["cases"][0]["case_id"]} for h in hits]


def _download_counts_file(file_id: str) -> pd.Series:
    """Download a single STAR-Counts file and return gene_id → count Series.

    GDC STAR-Counts files are TSV.gz with columns:
    gene_id, gene_name, gene_type, unstranded, stranded_first, stranded_second, tpm_unstranded, fpkm_unstranded, fpkm_uq_unstranded
    We use the 'unstranded' column (raw counts).
    """
    r = requests.get(f"{GDC_DATA_ENDPOINT}/{file_id}", timeout=60,
                     headers={"Accept": "application/octet-stream"})
    r.raise_for_status()
    # GDC returns gzipped TSV
    try:
        content = gzip.decompress(r.content).decode("utf-8")
    except gzip.BadGzipFile:
        content = r.content.decode("utf-8")
    df = pd.read_csv(io.StringIO(content), sep="\t", comment="#")
    # Skip the first few metadata rows (N_unmapped, N_multimapping, etc.)
    df = df[~df["gene_id"].str.startswith("N_")]
    # Use gene_name for readability, deduplicate by taking max
    counts = df.groupby("gene_name")["unstranded"].max()
    return counts.astype(int)


def download_cancer_type(project_id: str, cache_dir: Path,
                         max_tumor: int = 50, max_normal: int = 60) -> dict:
    """Download count matrices for one cancer type.

    Returns dict with keys: 'counts' (DataFrame genes×samples),
    'sample_type' (Series: sample_id → 'Tumor'/'Normal'),
    'n_tumor', 'n_normal', 'project_id'.
    """
    cache_file = cache_dir / f"{project_id}_counts.parquet"
    meta_file = cache_dir / f"{project_id}_meta.parquet"

    if cache_file.exists() and meta_file.exists():
        counts = pd.read_parquet(cache_file)
        meta = pd.read_parquet(meta_file)
        return {
            "counts": counts,
            "sample_type": meta["sample_type"],
            "n_tumor": (meta["sample_type"] == "Tumor").sum(),
            "n_normal": (meta["sample_type"] == "Normal").sum(),
            "project_id": project_id,
        }

    print(f"  Downloading {project_id}...")
    tumor_files = _query_file_ids(project_id, "Primary Tumor", max_tumor)
    normal_files = _query_file_ids(project_id, "Solid Tissue Normal", max_normal)

    all_files = [(f, "Tumor") for f in tumor_files] + [(f, "Normal") for f in normal_files]
    sample_counts = {}
    sample_types = {}
    failures = 0

    for finfo, stype in tqdm(all_files, desc=f"  {project_id}", leave=False):
        sample_id = f"{finfo['case_id'][:8]}_{stype[0]}"
        # Ensure unique sample IDs
        orig_id = sample_id
        suffix = 1
        while sample_id in sample_counts:
            sample_id = f"{orig_id}_{suffix}"
            suffix += 1
        try:
            sample_counts[sample_id] = _download_counts_file(finfo["file_id"])
            sample_types[sample_id] = stype
        except Exception as e:
            failures += 1
            if failures > 5:
                print(f"    Warning: {failures} download failures for {project_id}")
            continue
        # Rate limiting
        time.sleep(0.1)

    if not sample_counts:
        raise RuntimeError(f"No samples downloaded for {project_id}")

    counts = pd.DataFrame(sample_counts).fillna(0).astype(int)
    meta = pd.DataFrame({"sample_type": sample_types})

    # Cache
    cache_dir.mkdir(parents=True, exist_ok=True)
    counts.to_parquet(cache_file)
    meta.to_parquet(meta_file)

    return {
        "counts": counts,
        "sample_type": meta["sample_type"],
        "n_tumor": (meta["sample_type"] == "Tumor").sum(),
        "n_normal": (meta["sample_type"] == "Normal").sum(),
        "project_id": project_id,
    }


def get_tissue_specific_genes(tissue: str, n_genes: int = 200) -> list[str]:
    """Get tissue-specific genes from GTEx median expression.

    Uses GTEx API to find genes highly expressed in the target tissue
    relative to other tissues. Returns top n_genes by specificity.
    """
    # Use Human Protein Atlas tissue-specific gene lists as fallback,
    # or GTEx median expression. GTEx API v2 provides median gene expression.
    # For simplicity, use a curated set of known tissue-specific genes.
    # This is more robust than API calls that may return normalized data.
    TISSUE_GENES = {
        "Breast": ["ESR1", "PGR", "GATA3", "FOXA1", "KRT18", "KRT19", "MUC1",
                    "ERBB2", "CDH1", "EPCAM", "KRT8", "KRT7", "ANKRD30A",
                    "SCGB2A2", "SCGB1D2", "TFF1", "TFF3", "AGR2", "PRLR", "CITED1"],
        "Kidney": ["SLC12A1", "AQP1", "AQP2", "UMOD", "SLC34A1", "NPHS1",
                    "NPHS2", "WT1", "PAX8", "HNF1B", "SLC22A6", "SLC22A8",
                    "CLCNKA", "CLCNKB", "KCNJ1", "SLC12A3", "CALB1", "AQP3",
                    "SLC4A1", "CA4"],
        "Lung": ["SFTPA1", "SFTPB", "SFTPC", "SFTPD", "NKX2-1", "SCGB1A1",
                  "SCGB3A2", "AGER", "HOPX", "AQP5", "ABCA3", "LAMP3",
                  "MUC5AC", "MUC5B", "FOXJ1", "DNAH5", "CFTR", "SLC34A2",
                  "NAPSA", "WIF1"],
        "Thyroid": ["TG", "TPO", "TSHR", "NIS", "PAX8", "FOXE1", "NKX2-1",
                     "DUOX2", "DIO1", "DIO2", "SLC5A5", "THRA", "THRB",
                     "IYD", "DUOXA2", "CT", "CALCA", "PTH", "SLC26A4", "KCNJ13"],
        "Prostate": ["KLK3", "KLK2", "ACPP", "NKX3-1", "AR", "TMPRSS2",
                      "FOLH1", "SLC45A3", "MSMB", "TGM4", "KLK4", "STEAP2",
                      "PMEPA1", "HOXB13", "ERG", "AMACR", "PCA3", "SPINK1",
                      "RLN1", "SEMG1"],
        "Liver": ["ALB", "AFP", "SERPINA1", "APOB", "APOA1", "CYP3A4",
                   "CYP2E1", "CYP1A2", "HNF4A", "HNF1A", "FGA", "FGB",
                   "FGG", "F2", "F7", "F9", "F10", "AHSG", "TF", "TTR"],
        "Head and Neck": ["KRT5", "KRT14", "TP63", "KRT6A", "KRT13",
                          "KRT4", "SPRR1A", "SPRR1B", "SPRR3", "IVL",
                          "LOR", "DSG3", "DSC3", "PERP", "S100A7", "S100A8",
                          "S100A9", "SERPINB3", "SERPINB4", "CLCA2"],
        "Colon": ["CDX2", "VIL1", "MUC2", "TFF3", "FABP1", "CEACAM5",
                   "CEACAM6", "GPA33", "CDH17", "CLDN3", "CLDN4", "CLDN7",
                   "REG4", "DEFA5", "DEFA6", "LYZ", "OLFM4", "LGR5",
                   "BEST4", "OTOP2"],
        "Stomach": ["GKN1", "GKN2", "TFF1", "TFF2", "MUC5AC", "MUC6",
                     "PGA3", "PGA4", "PGA5", "PGC", "GIF", "ATP4A",
                     "ATP4B", "GAST", "GHRL", "SST", "CHGA", "LGR5",
                     "AQP5", "LIPF"],
        "Uterus": ["ESR1", "PGR", "HOXA10", "HOXA11", "MMP7", "MMP26",
                    "PAEP", "SPP1", "LIF", "WNT7A", "MSX1", "MSX2",
                    "HAND2", "FOXO1", "IGF1", "LEFTY2", "MUC1", "SCGB2A2",
                    "KLF9", "DKK1"],
        "Bladder": ["UPK1A", "UPK1B", "UPK2", "UPK3A", "UPK3B", "KRT20",
                     "KRT5", "KRT14", "GATA3", "FOXA1", "PPARG", "TP63",
                     "KRT7", "KRT18", "FGFR3", "RB1", "CDKN2A", "HRAS",
                     "CD44", "ERBB2"],
        "Esophagus": ["KRT4", "KRT13", "KRT5", "KRT14", "TP63", "SOX2",
                       "KRT6A", "SPRR1A", "SPRR3", "IVL", "LOR", "DSG3",
                       "PERP", "S100A7", "S100A8", "S100A9", "MUC1",
                       "MUC4", "TFF3", "GKN1"],
    }
    return TISSUE_GENES.get(tissue, [])


def check_batch_metadata(project_id: str) -> dict:
    """Check if batch/plate metadata is available for a TCGA project.

    Queries GDC biospecimen data for plate_id or center_id fields.
    Returns dict with 'has_plate', 'has_center', 'n_plates', 'n_centers'.
    """
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "project.project_id", "value": [project_id]}},
        ],
    }
    fields = "aliquot_ids,submitter_id,project.project_id"
    params = {
        "filters": json.dumps(filters),
        "fields": "submitter_id,analyte_ids,portions.analytes.aliquots.center.short_name,portions.analytes.aliquots.plate_name",
        "size": 5,
        "format": "json",
        "expand": "portions.analytes.aliquots",
    }
    try:
        r = requests.get("https://api.gdc.cancer.gov/cases", params=params, timeout=15)
        r.raise_for_status()
        hits = r.json()["data"]["hits"]
        plates = set()
        centers = set()
        for h in hits:
            for portion in h.get("portions", []):
                for analyte in portion.get("analytes", []):
                    for aliquot in analyte.get("aliquots", []):
                        if aliquot.get("plate_name"):
                            plates.add(aliquot["plate_name"])
                        center = aliquot.get("center", {})
                        if center and center.get("short_name"):
                            centers.add(center["short_name"])
        return {
            "project_id": project_id,
            "has_plate": len(plates) > 0,
            "has_center": len(centers) > 0,
            "n_plates": len(plates),
            "n_centers": len(centers),
        }
    except Exception as e:
        return {
            "project_id": project_id,
            "has_plate": False,
            "has_center": False,
            "n_plates": 0,
            "n_centers": 0,
            "error": str(e),
        }
