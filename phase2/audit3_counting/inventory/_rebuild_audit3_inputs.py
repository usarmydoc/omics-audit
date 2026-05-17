#!/usr/bin/env python3
"""Rebuild audit3_inputs.tsv per user direction 2026-05-17:
- Remove Mereu HCA benchmark (out of scope per user)
- Add 2 mouse datasets in 3' v2 chemistry from distinct tissues
- Add chemistry_exact column for v3 vs v3.1 stratification
- Keep all original rows otherwise unchanged
"""
from __future__ import annotations
import csv
from pathlib import Path

OUT = Path("/mnt/nvme1/omics-audit/phase2/audit3_counting/inventory/audit3_inputs.tsv")

# Column order matches CP1 spec + new chemistry_exact column
COLS = [
    "dataset_id", "source", "accession", "chemistry", "chemistry_exact",
    "tissue", "species", "n_samples", "n_cells_estimated",
    "fastq_size_gb_estimated", "size_confidence", "reference_genome_expected",
    "original_counter_used_by_submitter", "license_constraints",
    "raw_fastq_url_or_accession", "has_publication", "citation",
    "access_status", "acquisition_status", "chemistry_in_scope", "notes",
]

ROWS = [
    # --- 3' v3 (broad category) ---
    {
        "dataset_id": "10x_pbmc_1k_v3", "source": "10x_public",
        "accession": "pbmc_1k_v3", "chemistry": "3p_v3", "chemistry_exact": "v3",
        "tissue": "PBMC", "species": "human", "n_samples": 1, "n_cells_estimated": 1222,
        "fastq_size_gb_estimated": 5, "size_confidence": "estimated",
        "reference_genome_expected": "GRCh38",
        "original_counter_used_by_submitter": "cellranger",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/",
        "has_publication": "no", "citation": "10x Genomics demo data, Cell Ranger 3.0.0 release",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "Small entry-cost dataset; good smoke test for all four tools",
    },
    {
        "dataset_id": "10x_pbmc_5k_v3.1", "source": "10x_public",
        "accession": "5k_pbmc_v3_nextgem", "chemistry": "3p_v3", "chemistry_exact": "v3.1",
        "tissue": "PBMC", "species": "human", "n_samples": 1, "n_cells_estimated": 5140,
        "fastq_size_gb_estimated": 30, "size_confidence": "estimated",
        "reference_genome_expected": "GRCh38",
        "original_counter_used_by_submitter": "cellranger_7.0.1",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://www.10xgenomics.com/datasets/5k-human-pbmcs-3-v3-1-chromium-controller-3-1-standard",
        "has_publication": "yes (10x app note)",
        "citation": "10x Genomics 5k PBMCs 3' v3.1 Chromium Controller dataset",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "VERIFIED n_cells from catalog page; size estimated; R1=28cy R2=90cy ~35k reads/cell",
    },
    {
        "dataset_id": "10x_pbmc_10k_v3.1", "source": "10x_public",
        "accession": "10k_pbmc_v3_nextgem", "chemistry": "3p_v3", "chemistry_exact": "v3.1",
        "tissue": "PBMC", "species": "human", "n_samples": 1, "n_cells_estimated": 11485,
        "fastq_size_gb_estimated": 75, "size_confidence": "estimated",
        "reference_genome_expected": "GRCh38",
        "original_counter_used_by_submitter": "cellranger",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://www.10xgenomics.com/datasets/10k-human-pbmcs-3-v3-1-chromium-controller-3-1-high",
        "has_publication": "yes (10x app note)",
        "citation": "10x Genomics 10k PBMCs 3' v3.1 Chromium Controller dataset",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "VERIFIED n_cells from catalog page; size estimated",
    },
    {
        "dataset_id": "10x_neurons_900_v3", "source": "10x_public",
        "accession": "neurons_900_v3", "chemistry": "3p_v3", "chemistry_exact": "v3",
        "tissue": "brain", "species": "mouse", "n_samples": 1, "n_cells_estimated": 931,
        "fastq_size_gb_estimated": 5, "size_confidence": "estimated",
        "reference_genome_expected": "mm10",
        "original_counter_used_by_submitter": "cellranger",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://cf.10xgenomics.com/samples/cell-exp/3.0.0/neuron_900_v3/",
        "has_publication": "no", "citation": "10x Genomics demo data, Cell Ranger 3.0.0 release",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "MOUSE — can reuse /mnt/nvme2/refs/mm39 ref set per CP0 addendum",
    },

    # --- 3' v2 ---
    {
        "dataset_id": "10x_pbmc_3k_v2", "source": "10x_public",
        "accession": "pbmc_3k", "chemistry": "3p_v2", "chemistry_exact": "v2",
        "tissue": "PBMC", "species": "human", "n_samples": 1, "n_cells_estimated": 2700,
        "fastq_size_gb_estimated": 18, "size_confidence": "estimated",
        "reference_genome_expected": "GRCh38",
        "original_counter_used_by_submitter": "cellranger",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://cf.10xgenomics.com/samples/cell-exp/1.1.0/pbmc3k/",
        "has_publication": "no", "citation": "Classic 10x dataset, widely-used Seurat tutorial set",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "Canonical benchmark for 3' v2",
    },
    {
        "dataset_id": "10x_pbmc_4k_v2", "source": "10x_public",
        "accession": "pbmc_4k", "chemistry": "3p_v2", "chemistry_exact": "v2",
        "tissue": "PBMC", "species": "human", "n_samples": 1, "n_cells_estimated": 4340,
        "fastq_size_gb_estimated": 25, "size_confidence": "estimated",
        "reference_genome_expected": "GRCh38",
        "original_counter_used_by_submitter": "cellranger",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/",
        "has_publication": "no", "citation": "10x Genomics demo, Cell Ranger 2.1.0 release",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "Pairs with pbmc_3k for v2 chemistry replication",
    },
    {
        "dataset_id": "10x_heart_1k_v2_mouse", "source": "10x_public",
        "accession": "heart_1k_v2", "chemistry": "3p_v2", "chemistry_exact": "v2",
        "tissue": "heart", "species": "mouse", "n_samples": 1,
        "n_cells_estimated": 1011, "fastq_size_gb_estimated": 6,
        "size_confidence": "estimated",
        "reference_genome_expected": "mm10",
        "original_counter_used_by_submitter": "cellranger",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://cf.10xgenomics.com/samples/cell-exp/2.1.0/heart_1k_v2/",
        "has_publication": "no",
        "citation": "10x Genomics E18 mouse heart 1k cells demo, Cell Ranger 2.1.0",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "NEW per user direction 2026-05-17 — mouse heart, 3' v2; expands mouse coverage into v2 chemistry and tissue beyond brain",
    },
    {
        "dataset_id": "tabula_muris_liver_droplet", "source": "GEO",
        "accession": "GSE109774 (liver subset)",
        "chemistry": "3p_v2", "chemistry_exact": "v2",
        "tissue": "liver", "species": "mouse", "n_samples": 3,
        "n_cells_estimated": 1500, "fastq_size_gb_estimated": 25,
        "size_confidence": "estimated",
        "reference_genome_expected": "mm10",
        "original_counter_used_by_submitter": "cellranger_1.3 (per Tabula Muris paper)",
        "license_constraints": "open (CZ Biohub release)",
        "raw_fastq_url_or_accession": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109774 (filter to Liver tissue 10x droplet samples)",
        "has_publication": "yes",
        "citation": "Tabula Muris Consortium 2018, Nature 562:367–372 (DOI:10.1038/s41586-018-0590-4)",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "NEW per user direction 2026-05-17 — Tabula Muris is the standard mouse-tissue-diversity benchmark; liver subset chosen for 3' v2 chemistry replication with heart_1k_v2 and additional tissue coverage. Cell count and size estimated from Tabula Muris paper Table 1; need to confirm specific GSM accessions at download time (likely GSM2967037 + GSM2967038 + GSM2967039 or similar)",
    },

    # --- 5' v2 ---
    {
        "dataset_id": "10x_sc5p_v2_hs_PBMC_5k", "source": "10x_public",
        "accession": "sc5p_v2_hs_PBMC_5k", "chemistry": "5p_v2", "chemistry_exact": "5p_v2",
        "tissue": "PBMC", "species": "human", "n_samples": 1, "n_cells_estimated": 5000,
        "fastq_size_gb_estimated": 35, "size_confidence": "estimated",
        "reference_genome_expected": "GRCh38",
        "original_counter_used_by_submitter": "cellranger",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://www.10xgenomics.com/datasets/5k-human-pbmcs-5-v2-0",
        "has_publication": "yes", "citation": "10x Genomics 5' v2 5k PBMC dataset",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "Pairs with 10k_5p_v2 for chemistry replication",
    },
    {
        "dataset_id": "10x_sc5p_v2_hs_PBMC_10k", "source": "10x_public",
        "accession": "sc5p_v2_hs_PBMC_10k", "chemistry": "5p_v2", "chemistry_exact": "5p_v2",
        "tissue": "PBMC", "species": "human", "n_samples": 1, "n_cells_estimated": 10000,
        "fastq_size_gb_estimated": 70, "size_confidence": "estimated",
        "reference_genome_expected": "GRCh38",
        "original_counter_used_by_submitter": "cellranger",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://www.10xgenomics.com/datasets/10k-human-pbmcs-stained-with-totalseqc-human-universal-cocktail-v1-0-5-v2-0",
        "has_publication": "yes", "citation": "10x Genomics 5' v2 10k PBMC dataset",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "Size estimated; cell count rounded",
    },

    # --- Multiome RNA ---
    {
        "dataset_id": "10x_pbmc_10k_multiome", "source": "10x_public",
        "accession": "pbmc_granulocyte_sorted_10k_arc",
        "chemistry": "multiome_rna", "chemistry_exact": "multiome_v1.0",
        "tissue": "PBMC", "species": "human", "n_samples": 1, "n_cells_estimated": 10970,
        "fastq_size_gb_estimated": 50, "size_confidence": "estimated",
        "reference_genome_expected": "GRCh38",
        "original_counter_used_by_submitter": "cellranger_arc_2.0.0",
        "license_constraints": "CC-BY-4.0 (10x public)",
        "raw_fastq_url_or_accession": "https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-x-1-standard-2-0-0",
        "has_publication": "yes",
        "citation": "10x Genomics 10k PBMCs Multiome v1.0, Chromium X dataset",
        "access_status": "open", "acquisition_status": "needs_download",
        "chemistry_in_scope": "yes",
        "notes": "VERIFIED n_cells (10,970 nuclei); ~50k reads/cell GEX. ATAC libraries NOT downloaded (Audit 3 scopes only RNA modality of multiome).",
    },
]


def main():
    assert all(set(r.keys()) == set(COLS) for r in ROWS), "column mismatch"
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
        writer.writeheader()
        writer.writerows(ROWS)
    # Summary
    by_chem: dict[str, list] = {}
    by_species: dict[str, list] = {}
    by_tissue: dict[str, list] = {}
    total_gb = 0.0
    for r in ROWS:
        by_chem.setdefault(r["chemistry"], []).append(r)
        by_species.setdefault(r["species"], []).append(r)
        by_tissue.setdefault(r["tissue"], []).append(r)
        total_gb += r["fastq_size_gb_estimated"]
    print(f"Wrote {OUT}: {len(ROWS)} datasets, {len(COLS)} columns")
    print(f"\nBy chemistry:")
    for c, rs in sorted(by_chem.items()):
        print(f"  {c:15s}  n={len(rs)}")
    print(f"\nBy species:")
    for s, rs in sorted(by_species.items()):
        print(f"  {s:10s}  n={len(rs)}")
    print(f"\nBy tissue:")
    for t, rs in sorted(by_tissue.items()):
        print(f"  {t:10s}  n={len(rs)}")
    print(f"\nEstimated total FASTQ download: {total_gb:.0f} GB")


if __name__ == "__main__":
    main()
