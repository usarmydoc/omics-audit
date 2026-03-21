"""
B4 data access: load TCGA data with plate-level batch info from barcodes.

Extracts plate IDs from TCGA aliquot barcodes in the GDC file metadata.
Reuses count matrices from B3 cache.
"""

import json
from pathlib import Path

import pandas as pd
import requests

REPO_ROOT = Path(__file__).resolve().parent.parent.parent
B3_CACHE = REPO_ROOT / "output" / "b3_cache"

CANCER_TYPES = [
    "TCGA-BRCA", "TCGA-KIRC", "TCGA-LUAD", "TCGA-THCA", "TCGA-PRAD",
    "TCGA-LUSC", "TCGA-LIHC", "TCGA-HNSC", "TCGA-COAD", "TCGA-STAD",
    "TCGA-UCEC", "TCGA-KIRP", "TCGA-KICH", "TCGA-BLCA", "TCGA-ESCA",
    "TCGA-READ",
]


def _fetch_barcodes(project_id: str, max_files: int = 200) -> dict:
    """Fetch aliquot barcodes from GDC for STAR-Counts files.

    Returns dict mapping case_id_prefix → {barcode, plate, center, tss}.
    """
    filters = {
        "op": "and",
        "content": [
            {"op": "in", "content": {"field": "cases.project.project_id", "value": [project_id]}},
            {"op": "in", "content": {"field": "data_type", "value": ["Gene Expression Quantification"]}},
            {"op": "in", "content": {"field": "analysis.workflow_type", "value": ["STAR - Counts"]}},
        ],
    }
    params = {
        "filters": json.dumps(filters),
        "fields": "cases.case_id,cases.samples.sample_type,cases.samples.portions.analytes.aliquots.submitter_id",
        "size": max_files,
        "format": "json",
    }
    r = requests.get("https://api.gdc.cancer.gov/files", params=params, timeout=30)
    r.raise_for_status()

    barcode_map = {}
    for hit in r.json()["data"]["hits"]:
        for case in hit.get("cases", []):
            case_id = case["case_id"]
            case_prefix = case_id[:8]
            for sample in case.get("samples", []):
                sample_type = sample.get("sample_type", "")
                stype_code = "T" if "Tumor" in sample_type else "N"
                for portion in sample.get("portions", []):
                    for analyte in portion.get("analytes", []):
                        for aliquot in analyte.get("aliquots", []):
                            bc = aliquot.get("submitter_id", "")
                            if not bc.startswith("TCGA-"):
                                continue
                            parts = bc.split("-")
                            if len(parts) >= 7:
                                plate = parts[5]
                                center = parts[6]
                                tss = parts[1]
                            else:
                                plate = "UNKNOWN"
                                center = "UNKNOWN"
                                tss = parts[1] if len(parts) > 1 else "UNKNOWN"
                            key = f"{case_prefix}_{stype_code}"
                            barcode_map[key] = {
                                "barcode": bc,
                                "plate": plate,
                                "center": center,
                                "tss": tss,
                                "sample_type": sample_type,
                            }
    return barcode_map


def load_with_batch(project_id: str, cache_dir: Path = None) -> dict:
    """Load count matrix + batch metadata for one cancer type.

    Returns dict with 'counts', 'sample_type', 'plate', 'tss', 'n_plates'.
    Matches GDC barcodes to B3 cache sample IDs using case_id prefix.
    """
    counts_file = B3_CACHE / f"{project_id}_counts.parquet"
    meta_file = B3_CACHE / f"{project_id}_meta.parquet"

    if not counts_file.exists():
        raise FileNotFoundError(f"B3 cache not found for {project_id}")

    counts = pd.read_parquet(counts_file)
    meta = pd.read_parquet(meta_file)

    # Fetch barcodes from GDC
    barcode_map = _fetch_barcodes(project_id)

    # Match B3 sample IDs (format: case_prefix_T or case_prefix_N) to barcodes
    plates = {}
    tss_codes = {}
    for sample_id in counts.columns:
        if sample_id in barcode_map:
            plates[sample_id] = barcode_map[sample_id]["plate"]
            tss_codes[sample_id] = barcode_map[sample_id]["tss"]
        else:
            # Try matching by prefix
            matched = False
            for key, info in barcode_map.items():
                if sample_id.startswith(key[:8]):
                    plates[sample_id] = info["plate"]
                    tss_codes[sample_id] = info["tss"]
                    matched = True
                    break
            if not matched:
                plates[sample_id] = "UNKNOWN"
                tss_codes[sample_id] = "UNKNOWN"

    plate_series = pd.Series(plates)
    tss_series = pd.Series(tss_codes)

    # Filter out UNKNOWN plates
    known = plate_series[plate_series != "UNKNOWN"].index
    if len(known) < len(counts.columns) * 0.5:
        print(f"    WARNING: Only {len(known)}/{len(counts.columns)} samples matched to plates")

    n_plates = plate_series[plate_series != "UNKNOWN"].nunique()

    return {
        "counts": counts,
        "sample_type": meta["sample_type"],
        "plate": plate_series,
        "tss": tss_series,
        "n_plates": n_plates,
        "n_tss": tss_series[tss_series != "UNKNOWN"].nunique(),
        "project_id": project_id,
    }
