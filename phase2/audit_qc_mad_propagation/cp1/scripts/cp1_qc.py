#!/usr/bin/env python
"""CP1 Step 2 — compute QC metrics + apply the 4 QC methods to a pulled dataset.

Reuses QC-MAD's exact code paths: C1 (5th-pct floors + 95th mito), C2 (200/500 +
95th mito), MAD k=3/5 via QC-MAD's mad_filter.R (scuttle::isOutlier log=TRUE).
Writes per-method survivor soma_joinid lists and verifies counts vs QC-MAD saved.
Run in audit3_counting env. Usage: cp1_qc.py <short>
"""
import sys, subprocess
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc

short = sys.argv[1]
CP1 = Path("/mnt/nvme1/omics-audit/phase2/audit_qc_mad_propagation/cp1")
QCMAD = Path("/mnt/nvme1/omics-audit/phase2/audit_qc_mad/cp1")
(CP1 / "survivors").mkdir(exist_ok=True)

a = sc.read_h5ad(CP1 / "matrices" / f"{short}_raw.h5ad")
a.var["mt"] = a.var["feature_name"].astype(str).str.upper().str.startswith(("MT-", "MT."))
sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None)
df = pd.DataFrame({"cell": a.obs_names.astype(int),
                   "n_genes_by_counts": a.obs.n_genes_by_counts.values,
                   "total_counts": a.obs.total_counts.values,
                   "pct_counts_mt": a.obs.pct_counts_mt.values})
metrics_tsv = CP1 / "survivors" / f"{short}_metrics.tsv"
df.to_csv(metrics_tsv, sep="\t", index=False)

# MAD via QC-MAD's mad_filter.R
mad_out = CP1 / "survivors" / f"{short}_mad.tsv"
subprocess.run(["Rscript", str(QCMAD / "mad_filter.R"), str(metrics_tsv), str(mad_out)], check=True)
mad = pd.read_csv(mad_out, sep="\t").set_index("cell")

g, c, m = df.n_genes_by_counts.values, df.total_counts.values, df.pct_counts_mt.values
cells = df.cell.values
g5, c5, m95 = np.percentile(g, 5), np.percentile(c, 5), np.percentile(m, 95)
surv = {
    "C1": cells[(g >= g5) & (c >= c5) & (m <= m95)],
    "C2": cells[(g >= 200) & (c >= 500) & (m <= m95)],
    "MAD3": cells[~mad.loc[cells, "mad3_outlier"].values],
    "MAD5": cells[~mad.loc[cells, "mad5_outlier"].values],
}
# verify counts vs QC-MAD saved
tgt = pd.read_csv(QCMAD / "per_dataset_metrics.tsv", sep="\t")
tgt = tgt[tgt["short"] == short].iloc[0]
print(f"[{short}] cell counts (got vs QC-MAD saved):")
ok = True
for meth in ["C1", "C2", "MAD3", "MAD5"]:
    got = len(surv[meth]); exp = int(tgt[f"n_cells_{meth}"])
    flag = "OK" if got == exp else "MISMATCH!!"
    if got != exp: ok = False
    print(f"  {meth}: {got} vs {exp}  {flag}")
    np.savetxt(CP1 / "survivors" / f"{short}_{meth}.txt", surv[meth], fmt="%d")
print(f"[{short}] QC methods {'reproduce QC-MAD exactly' if ok else 'MISMATCH — investigate'}")
sys.exit(0 if ok else 1)
