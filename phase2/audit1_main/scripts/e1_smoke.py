#!/usr/bin/env python3
"""E1 smoke test — 1 input from each category × 3 tools = 9 runs.

Validates the full e1_run.py pipeline before launching 375-run production.
Picks one representative input per category by alphabetical first entry.
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from e1_run import (
    OUT_DIR, RSCRATCH, write_fgsea_helper, write_cp_ora_helper,
    load_hallmark, detect_organism, detect_id_system,
    run_gseapy, run_r_tool, FGSEA_SCRIPT, CP_ORA_SCRIPT,
)
from pathlib import Path as P

INVENTORY = P("/mnt/nvme1/omics-audit/phase2/audit1_main/datasets/audit1_main_inputs.tsv")
DB_DIR = P("/mnt/nvme1/omics-audit/phase2/audit1_main/databases")
SMOKE_OUT = P("/mnt/nvme1/omics-audit/phase2/audit1_main/e1/smoke")

def main():
    SMOKE_OUT.mkdir(parents=True, exist_ok=True)
    write_fgsea_helper()
    write_cp_ora_helper()
    inv = pd.read_csv(INVENTORY, sep="\t")

    # Pick 1 representative per category
    picks = []
    for cat in ["tcga_cancer", "census_scrna", "gtex_tissue_pair"]:
        sub = inv[inv["input_category"] == cat].iloc[0]
        picks.append(sub)

    hm_cache = {
        "human": (load_hallmark("human"), DB_DIR / "msigdb_hallmark_human.tsv"),
        "mouse": (load_hallmark("mouse"), DB_DIR / "msigdb_hallmark_mouse.tsv"),
    }

    print(f"{'input':40s} {'tool':22s} {'n_pwy':>6s} {'sec':>6s}")
    print("-" * 80)
    for row in picks:
        deg = pd.read_csv(row["file_path"], sep="\t")
        organism = detect_organism(row)
        id_system = detect_id_system(row, deg)
        hm_df, hm_path = hm_cache[organism]
        for tool in ["fgsea", "gseapy_enrichr", "clusterProfiler_ORA"]:
            out = SMOKE_OUT / f"{row['input_id']}__{tool}.tsv"
            import time
            t0 = time.time()
            if tool == "fgsea":
                n = run_r_tool(FGSEA_SCRIPT, row["file_path"], hm_path, id_system, out, row["input_id"])
            elif tool == "clusterProfiler_ORA":
                n = run_r_tool(CP_ORA_SCRIPT, row["file_path"], hm_path, id_system, out, row["input_id"])
            else:
                n = run_gseapy(deg, hm_df, id_system, row["input_id"], out)
            print(f"{row['input_id'][:39]:40s} {tool:22s} {n:6d} {time.time()-t0:6.1f}")

if __name__ == "__main__":
    main()
