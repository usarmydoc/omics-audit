"""
B1 data access: reuse TCGA count matrices from B3 cache.
"""

from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
B3_CACHE = REPO_ROOT / "output" / "b3_cache"

CANCER_TYPES = [
    "TCGA-BRCA", "TCGA-KIRC", "TCGA-LUAD", "TCGA-THCA", "TCGA-PRAD",
    "TCGA-LUSC", "TCGA-LIHC", "TCGA-HNSC", "TCGA-COAD", "TCGA-STAD",
    "TCGA-UCEC", "TCGA-KIRP", "TCGA-KICH", "TCGA-BLCA", "TCGA-ESCA",
    "TCGA-READ",
]

CANCER_TO_TISSUE = {
    "TCGA-BRCA": "Breast", "TCGA-KIRC": "Kidney", "TCGA-LUAD": "Lung",
    "TCGA-THCA": "Thyroid", "TCGA-PRAD": "Prostate", "TCGA-LUSC": "Lung",
    "TCGA-LIHC": "Liver", "TCGA-HNSC": "Head and Neck",
    "TCGA-COAD": "Colon", "TCGA-STAD": "Stomach", "TCGA-UCEC": "Uterus",
    "TCGA-KIRP": "Kidney", "TCGA-KICH": "Kidney", "TCGA-BLCA": "Bladder",
    "TCGA-ESCA": "Esophagus", "TCGA-READ": "Colon",
}


def load_cancer_type(project_id: str) -> dict:
    """Load cached count matrix and sample metadata from B3."""
    counts_file = B3_CACHE / f"{project_id}_counts.parquet"
    meta_file = B3_CACHE / f"{project_id}_meta.parquet"

    if not counts_file.exists():
        raise FileNotFoundError(f"B3 cache not found for {project_id}")

    counts = pd.read_parquet(counts_file)
    meta = pd.read_parquet(meta_file)

    return {
        "counts": counts,
        "sample_type": meta["sample_type"],
        "n_tumor": (meta["sample_type"] == "Tumor").sum(),
        "n_normal": (meta["sample_type"] == "Normal").sum(),
        "project_id": project_id,
    }
