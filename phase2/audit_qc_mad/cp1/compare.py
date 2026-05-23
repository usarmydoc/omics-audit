#!/usr/bin/env python3
"""CP1 step 3: apply 4 filtering methods to each dataset's QC metrics, compute
pairwise agreement + disagreement concentration + dataset-level bootstrap.
Methods: C1 (quantile-data 5th/5th/95th), C2 (fixed floors 200/500 + 95th mito),
MAD3, MAD5 (from mad_filter.R). Run in audit3_counting (scipy/pandas) or base."""
import sys, subprocess
from pathlib import Path
from itertools import combinations
import numpy as np, pandas as pd

BASE = Path("/mnt/nvme1/omics-audit/phase2/audit_qc_mad/cp1")
QC = BASE / "qc_metrics"
MADDIR = BASE / "mad_flags"; MADDIR.mkdir(exist_ok=True)
sys.path.insert(0, "/mnt/nvme1/omics-audit/phase2/scripts")
from dge_native import register_output

PAIRS = [("C1", "MAD3"), ("C1", "MAD5"), ("C2", "MAD3"), ("C2", "MAD5"),
         ("C1", "C2"), ("MAD3", "MAD5")]

def retained_sets(df, mad):
    cells = df["cell"].values
    g, c, m = df["n_genes_by_counts"].values, df["total_counts"].values, df["pct_counts_mt"].values
    g5, c5 = np.percentile(g, 5), np.percentile(c, 5)
    m95 = np.percentile(m, 95)
    C1 = set(cells[(g >= g5) & (c >= c5) & (m <= m95)])
    C2 = set(cells[(g >= 200) & (c >= 500) & (m <= m95)])
    md = mad.set_index("cell")
    MAD3 = set(cells[~md.loc[cells, "mad3_outlier"].values])
    MAD5 = set(cells[~md.loc[cells, "mad5_outlier"].values])
    return {"C1": C1, "C2": C2, "MAD3": MAD3, "MAD5": MAD5}

def concentration(df, only_cells):
    """classify disagreeing cells by which QC extreme they hit."""
    sub = df[df["cell"].isin(only_cells)]
    g5 = np.percentile(df["n_genes_by_counts"], 5); c5 = np.percentile(df["total_counts"], 5)
    m95 = np.percentile(df["pct_counts_mt"], 95)
    low_g = sub["n_genes_by_counts"] < g5
    low_c = sub["total_counts"] < c5
    hi_m = sub["pct_counts_mt"] > m95
    n_multi = ((low_g.astype(int) + low_c.astype(int) + hi_m.astype(int)) >= 2).sum()
    return {"low_genes": int(low_g.sum()), "low_counts": int(low_c.sum()),
            "high_mito": int(hi_m.sum()), "multiple": int(n_multi),
            "none_of_3": int(len(sub) - (low_g | low_c | hi_m).sum())}

def main():
    ws = pd.read_csv(BASE / "working_set.tsv", sep="\t")
    rows, jac_rows, conc_rows = [], [], []
    failures = []
    for _, r in ws.iterrows():
        short, tissue = r["short"], r["tissue"]
        f = QC / f"{short}.tsv"
        if not f.exists():
            failures.append(f"{short}: QC metrics missing"); continue
        df = pd.read_csv(f, sep="\t")
        # run MAD via scuttle (audit3_counting env has scuttle)
        madf = MADDIR / f"{short}.mad.tsv"
        if not madf.exists():
            cmd = ["conda", "run", "-n", "audit3_counting", "Rscript",
                   str(BASE / "mad_filter.R"), str(f), str(madf)]
            out = subprocess.run(cmd, capture_output=True, text=True)
            if out.returncode != 0:
                failures.append(f"{short}: mad_filter.R failed: {out.stderr[-300:]}"); continue
            print("  " + out.stdout.strip(), flush=True)
        mad = pd.read_csv(madf, sep="\t")
        S = retained_sets(df, mad)
        rows.append(dict(dataset_id=r["dataset_id"], short=short, tissue=tissue,
                         n_cells_raw=len(df), n_cells_C1=len(S["C1"]), n_cells_C2=len(S["C2"]),
                         n_cells_MAD3=len(S["MAD3"]), n_cells_MAD5=len(S["MAD5"])))
        for a, b in PAIRS:
            A, B = S[a], S[b]
            inter, union = A & B, A | B
            jac = len(inter) / len(union) if union else float("nan")
            jac_rows.append(dict(dataset_id=r["dataset_id"], short=short, tissue=tissue,
                                 pair=f"{a}_vs_{b}", jaccard=round(jac, 4),
                                 n_intersection=len(inter), n_union=len(union),
                                 n_method_a_only=len(A - B), n_method_b_only=len(B - A),
                                 more_permissive=a if len(A) > len(B) else (b if len(B) > len(A) else "tie")))
            conc = concentration(df, (A - B) | (B - A))
            for cat, n in conc.items():
                conc_rows.append(dict(short=short, tissue=tissue, pair=f"{a}_vs_{b}",
                                      concentration_category=cat, n_cells=n))
        print(f"  [{short}] C1={len(S['C1'])} C2={len(S['C2'])} MAD3={len(S['MAD3'])} MAD5={len(S['MAD5'])} (raw {len(df)})", flush=True)

    dfm = pd.DataFrame(rows); dfj = pd.DataFrame(jac_rows); dfc = pd.DataFrame(conc_rows)
    dfm.to_csv(BASE / "per_dataset_metrics.tsv", sep="\t", index=False)
    dfj.to_csv(BASE / "per_dataset_jaccards.tsv", sep="\t", index=False)
    dfc.to_csv(BASE / "disagreement_concentration.tsv", sep="\t", index=False)

    # bootstrap dataset-level (all 8) per pair + high/low-mito split
    rng = np.random.default_rng(42)
    hi_mito = {"liver","small_intestine","pancreas","large_intestine","lung","blood"}  # group; see findings note
    boot = []
    for pair in dfj["pair"].unique():
        sub = dfj[dfj["pair"] == pair]
        for stratum, mask in [("all", sub.index),
                              ("high_mito", sub[sub.tissue.isin(["small intestine","pancreas","large intestine"])].index),
                              ("low_mito", sub[sub.tissue.isin(["heart","bone marrow","liver","blood","lung"])].index)]:
            vals = sub.loc[mask, "jaccard"].dropna().to_numpy()
            if len(vals) == 0: continue
            draws = rng.choice(vals, size=(1000, len(vals)), replace=True).mean(axis=1)
            boot.append(dict(stratum=stratum, pair=pair, bootstrap_b=1000, n_datasets=len(vals),
                             jaccard_mean=round(float(vals.mean()),4),
                             jaccard_ci_lower=round(float(np.percentile(draws,2.5)),4),
                             jaccard_ci_upper=round(float(np.percentile(draws,97.5)),4)))
    pd.DataFrame(boot).to_csv(BASE / "per_stratum_bootstrap.tsv", sep="\t", index=False)

    fm = BASE / "tool_failure_modes.md"
    with fm.open("w") as fh:
        fh.write("# Audit QC-MAD CP1 — tool failure modes\n\n")
        if failures:
            for x in failures: fh.write(f"- {x}\n")
        else:
            fh.write("No method failed on any dataset. All 8 datasets produced 4-method outputs.\n")
        fh.write("\n## MAD log-transform note\n\nscuttle::isOutlier(log=TRUE) on n_genes_by_counts + "
                 "total_counts (lower tail); pct_counts_mt raw (upper tail). nmads 3 and 5.\n")

    for p, k in [("per_dataset_metrics.tsv","audit_qcmad_per_dataset_metrics"),
                 ("per_dataset_jaccards.tsv","audit_qcmad_jaccards"),
                 ("per_stratum_bootstrap.tsv","audit_qcmad_bootstrap"),
                 ("disagreement_concentration.tsv","audit_qcmad_concentration"),
                 ("tool_failure_modes.md","audit_qcmad_failure_modes")]:
        register_output(BASE / p, kind=k)
    print("compare done; outputs + registered", flush=True)

if __name__ == "__main__":
    main()
