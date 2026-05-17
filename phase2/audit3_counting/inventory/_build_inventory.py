#!/usr/bin/env python3
"""Build local_fastq_inventory.tsv from the CP0 scan output.

Reads /tmp/audit3_cp0_all.txt (one path per line) and emits a TSV with
columns required by the CP0 spec plus a category/disposition column.

User-excluded folders (2026-05-17): any file under one of the following
paths is categorized 'user_excluded' with disposition 'skip':
  /mnt/nvme1/12282025Jihoon
  /mnt/nvme1/EGSEA_01122026Jiuen
  /mnt/nvme1/miRNA_01052026Ang
  /mnt/nvme2/Agnieszka 02-20-2026
  /mnt/nvme2/Katherine
(Equivalent via Desktop symlinks ~/Desktop/nvme1, ~/Desktop/nvme2.)
"""
from __future__ import annotations
import csv
from datetime import datetime
from pathlib import Path

SCAN = Path("/tmp/audit3_cp0_all.txt")
OUT  = Path("/mnt/nvme1/omics-audit/phase2/audit3_counting/inventory/local_fastq_inventory.tsv")

EXCLUDE_PREFIXES = [
    "/mnt/nvme1/12282025Jihoon",
    "/mnt/nvme1/EGSEA_01122026Jiuen",
    "/mnt/nvme1/miRNA_01052026Ang",
    "/mnt/nvme2/Agnieszka 02-20-2026",
    "/mnt/nvme2/Katherine",
]

# Test-fixture path patterns — match the R-package and SPAdes test data
TEST_FIXTURE_SUBSTRINGS = [
    "/R/x86_64-pc-linux-gnu-library/",   # any R package extdata/testdata/unitTests
    "/miniforge3/pkgs/spades-",
    "/miniforge3/envs/bio-linux/share/spades/",
]


def categorize(path: str, size_bytes: int) -> tuple[str, str, str]:
    """Return (category, disposition, reason)."""
    # User exclusion takes precedence
    for prefix in EXCLUDE_PREFIXES:
        if path.startswith(prefix):
            return ("user_excluded", "skip",
                    "Folder on user's 2026-05-17 exclude list — do not use.")
    # Test fixtures
    for s in TEST_FIXTURE_SUBSTRINGS:
        if s in path:
            return ("test_fixture", "skip",
                    "Library/tool test data, not real biological data.")
    # Anything else surviving — flag as unknown for manual review
    return ("unclassified", "skip",
            "Did not match any known category; defaulting to skip pending review.")


def main():
    rows = []
    with SCAN.open() as fh:
        paths = [line.strip() for line in fh if line.strip()]
    for p_str in sorted(paths):
        p = Path(p_str)
        try:
            stat = p.stat()
            size_bytes = stat.st_size
            mtime = datetime.fromtimestamp(stat.st_mtime).strftime("%Y-%m-%d")
        except FileNotFoundError:
            size_bytes = 0
            mtime = ""
        size_gb = round(size_bytes / 1e9, 4)
        parent = str(p.parent)
        # No reliable likely_dataset_id, chemistry_hint, species_hint for these
        # skipped files; leave blank rather than guess. md5 sibling check:
        md5_sibling = p.parent / f"{p.name}.md5"
        category, disposition, reason = categorize(p_str, size_bytes)
        # appears_complete: assume yes (filename has no .partial or .tmp suffix
        # and we won't process these anyway)
        appears_complete = "yes" if not (p_str.endswith(".partial")
                                          or p_str.endswith(".tmp")) else "no"
        rows.append({
            "file_path": p_str,
            "file_size_gb": size_gb,
            "mtime": mtime,
            "parent_dir": parent,
            "likely_dataset_id": "",
            "chemistry_hint": "",
            "species_hint": "",
            "appears_complete": appears_complete,
            "md5_available": "yes" if md5_sibling.exists() else "no",
            "category": category,
            "disposition": disposition,
            "reason": reason,
        })

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with OUT.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # Summary stats
    by_cat: dict[str, list[dict]] = {}
    for r in rows:
        by_cat.setdefault(r["category"], []).append(r)
    print(f"Wrote {OUT} with {len(rows)} rows")
    print(f"Category breakdown:")
    for cat, rs in sorted(by_cat.items()):
        total_gb = sum(r["file_size_gb"] for r in rs)
        print(f"  {cat:20s}  n={len(rs):3d}  total={total_gb:7.2f} GB")


if __name__ == "__main__":
    main()
