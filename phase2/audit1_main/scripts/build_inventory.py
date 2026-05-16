#!/usr/bin/env python3
"""Audit 1 main — Checkpoint 1: build inventory of DEG inputs.

For each per-gene DEG TSV in the three input categories, capture:
  input_id, source, comparison_id, tool, file_path, sha256, row_count,
  gene_id_format, has_padj, has_log2fc, generation_date, phase_lock_status,
  input_category, n_significant_padj_05

Also re-registers any input file whose hash drifted from the lock entry, and
reports verification failures + sample-size adequacy per sub-audit.
"""

from __future__ import annotations

import hashlib
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

LOCK = Path("/mnt/nvme1/omics-audit/phase2/repro.lock")
OUT = Path("/mnt/nvme1/omics-audit/phase2/audit1_main/datasets/audit1_main_inputs.tsv")
REPORT = Path("/mnt/nvme1/omics-audit/phase2/audit1_main/datasets/CP1_report.md")

B1_DIR = Path("/mnt/nvme1/omics-audit/phase2/per_gene_degs/b1_tcga")
P4_DIR = Path("/mnt/nvme1/omics-audit/phase2/per_gene_degs/p4_census")
GTEX_DIR = Path("/mnt/nvme1/omics-audit/phase2/per_gene_degs/gtex_pairs")

# Min sig genes for E1/E2 sanity (we want enough signal for enrichment to be meaningful)
SANITY_MIN_SIG = 50


def sha256_of(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def detect_gene_id_format(genes: pd.Series) -> str:
    sample = genes.dropna().astype(str).head(50).tolist()
    if any(g.startswith("ENSG") for g in sample): return "ensembl_human"
    if any(g.startswith("ENSMUSG") for g in sample): return "ensembl_mouse"
    if any(g.startswith("ENS") for g in sample): return "ensembl_other"
    if all(re.match(r"^[A-Z][A-Z0-9-]*$", g) for g in sample[:10]): return "symbol_human"
    if all(re.match(r"^[A-Z][a-z0-9-]*$", g) for g in sample[:10]): return "symbol_mouse"
    return "mixed_or_unknown"


def parse_filename(name: str, category: str) -> tuple[str, str]:
    """Return (comparison_id, tool) from filename."""
    stem = name.removesuffix(".tsv")
    if category == "tcga_cancer":
        # BLCA_DESeq2 → cancer=BLCA, tool=DESeq2
        parts = stem.split("_", 1)
        return parts[0], parts[1] if len(parts) > 1 else "unknown"
    if category == "census_scrna":
        # 16023185_pseudobulk → ds=16023185, tool=pseudobulk
        parts = stem.rsplit("_", 1)
        return parts[0], parts[1] if len(parts) > 1 else "unknown"
    if category == "gtex_tissue_pair":
        # GTEx_liver_vs_kidney_cortex_DESeq2 → pair=GTEx_liver_vs_kidney_cortex, tool=DESeq2
        # last underscore-segment is the tool
        idx = stem.rfind("_")
        return stem[:idx], stem[idx + 1:]
    return stem, "unknown"


def main() -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    lock = json.loads(LOCK.read_text())
    vo = lock["verified_outputs"]

    rows: list[dict] = []
    failures: list[str] = []

    categories = [
        ("tcga_cancer", B1_DIR),
        ("census_scrna", P4_DIR),
        ("gtex_tissue_pair", GTEX_DIR),
    ]

    for category, d in categories:
        for f in sorted(d.glob("*.tsv")):
            rel = str(f.resolve())
            try:
                actual_sha = sha256_of(f)
                lock_entry = vo.get(rel)
                if lock_entry is None:
                    lock_status = "NOT_IN_LOCK"
                elif lock_entry["sha256"] == actual_sha:
                    lock_status = "verified"
                else:
                    lock_status = f"HASH_MISMATCH(lock={lock_entry['sha256'][:8]}, actual={actual_sha[:8]})"
                    failures.append(f"{f}: hash mismatch")

                # Parse the TSV
                df = pd.read_csv(f, sep="\t")
                row_count = len(df)
                has_padj = "padj" in df.columns
                has_log2fc = "log2FC" in df.columns
                gene_id_format = detect_gene_id_format(df.get("gene_id", pd.Series([], dtype=str)))
                n_sig = int((df["padj"].fillna(1) < 0.05).sum()) if has_padj else 0
                comparison_id, tool = parse_filename(f.name, category)
                input_id = f"{category}__{f.stem}"
                gen_date = (lock_entry or {}).get("registered_at", datetime.fromtimestamp(f.stat().st_mtime).isoformat(timespec="seconds"))
                source = f"phase2 regen ({category})"
                rows.append({
                    "input_id": input_id,
                    "source": source,
                    "comparison_id": comparison_id,
                    "tool": tool,
                    "file_path": rel,
                    "sha256": actual_sha,
                    "row_count": row_count,
                    "gene_id_format": gene_id_format,
                    "has_padj": has_padj,
                    "has_log2fc": has_log2fc,
                    "generation_date": gen_date,
                    "phase_lock_status": lock_status,
                    "input_category": category,
                    "n_significant_padj_05": n_sig,
                })
            except Exception as e:
                failures.append(f"{f}: parse error {e}")
                continue

    df_out = pd.DataFrame(rows)
    df_out.to_csv(OUT, sep="\t", index=False)
    print(f"Wrote {len(df_out)} rows to {OUT}")

    # Build report
    lines = ["# Audit 1 main — Checkpoint 1 report\n"]
    lines.append(f"_Generated: {datetime.now().isoformat(timespec='seconds')}_\n")
    lines.append("## Counts by input category\n")
    lines.append("| Category | Total files | Verified | Hash mismatch | Not in lock | n with ≥50 sig | Median sig genes |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|")
    for cat in ["tcga_cancer", "census_scrna", "gtex_tissue_pair"]:
        sub = df_out[df_out["input_category"] == cat]
        n_verified = int((sub["phase_lock_status"] == "verified").sum())
        n_mismatch = int(sub["phase_lock_status"].str.startswith("HASH_MISMATCH").sum())
        n_missing = int((sub["phase_lock_status"] == "NOT_IN_LOCK").sum())
        n_sane = int((sub["n_significant_padj_05"] >= SANITY_MIN_SIG).sum())
        med_sig = int(sub["n_significant_padj_05"].median()) if len(sub) else 0
        lines.append(f"| {cat} | {len(sub)} | {n_verified} | {n_mismatch} | {n_missing} | {n_sane} | {med_sig} |")

    lines.append("\n## Sample size adequacy per sub-audit (per AUDIT_STANDARDS.md §5.3)")
    lines.append("\nConditional confidence tier requires ≥10 per category. Hard default requires ≥15 across ≥3 tissues/categories.\n")
    lines.append("| Category | n inputs | n comparisons (unique) | n tools | OK for conditional (≥10)? |")
    lines.append("|---|---:|---:|---:|:---:|")
    for cat in ["tcga_cancer", "census_scrna", "gtex_tissue_pair"]:
        sub = df_out[df_out["input_category"] == cat]
        n_comp = int(sub["comparison_id"].nunique())
        n_tool = int(sub["tool"].nunique())
        ok = "✓" if n_comp >= 10 else "✗"
        lines.append(f"| {cat} | {len(sub)} | {n_comp} | {n_tool} | {ok} |")

    if failures:
        lines.append("\n## Verification failures\n")
        for f in failures:
            lines.append(f"- {f}")
    else:
        lines.append("\n## Verification: PASS — all files parse cleanly, no hash mismatches.")

    # Gene ID format breakdown
    lines.append("\n## Gene ID format coverage")
    lines.append("\n| Category | symbol_human | ensembl_human | ensembl_mouse | symbol_mouse | other |")
    lines.append("|---|---:|---:|---:|---:|---:|")
    for cat in ["tcga_cancer", "census_scrna", "gtex_tissue_pair"]:
        sub = df_out[df_out["input_category"] == cat]
        cells = []
        for fmt in ["symbol_human", "ensembl_human", "ensembl_mouse", "symbol_mouse"]:
            cells.append(str(int((sub["gene_id_format"] == fmt).sum())))
        other = int((~sub["gene_id_format"].isin(["symbol_human", "ensembl_human", "ensembl_mouse", "symbol_mouse"])).sum())
        cells.append(str(other))
        lines.append(f"| {cat} | " + " | ".join(cells) + " |")

    REPORT.write_text("\n".join(lines) + "\n")
    print(f"Wrote report to {REPORT}")

    # Register both outputs in the lock
    sys.path.insert(0, "/mnt/nvme1/omics-audit/phase2/scripts")
    from dge_native import register_output
    register_output(OUT, kind="audit1_main_inventory", n_inputs=len(df_out))
    register_output(REPORT, kind="audit1_main_cp1_report")
    print("Registered inventory + report in lock.")

    if failures:
        print(f"\nFAILURES: {len(failures)} — see report")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
