#!/usr/bin/env python3
"""Audit 1 main CP4 — E2 metrics + bootstrap CIs.

E2 is structurally different from E1: pathway sets differ between
databases (Hallmark has 50 pathways, C5 GO:BP has 7000+), so direct
name-based Jaccard between databases is not meaningful. The comparison
question is: "do different databases produce the same biological
story?"

Per-database metrics (per input):
  - n_pathways_returned
  - n_sig_padj_05, n_sig_padj_01
  - frac_sig_padj_05 = n_sig / n_pathways
  - median_abs_NES (effect size of significant pathways)

Cross-database metrics (per (database_pair × input)) — story alignment
not pathway-name alignment:
  - leading_edge_gene_jaccard — Jaccard of the union of leading-edge
    genes from significant (padj<0.05) pathways. If two databases
    capture the same biology, the genes driving their significant
    pathways should overlap heavily.
  - direction_concordance — fraction of pathways significant in BOTH
    databases (by name match where possible) that have the same NES
    sign. When direct name matching is sparse, fall back to "fraction
    of significant pathways with NES > 0" in each, and report
    consistency of that sign pattern.

Within-database metrics (one per database):
  - mean_pathway_size
  - mean_pairwise_jaccard_5pct — pathway redundancy: mean Jaccard
    of gene-set memberships between a 5% random sample of pathway
    pairs (full O(N²) is too expensive for GO:BP)

Database coverage gap (per (database_pair × input)):
  - genes_in_sig_pathways_uniqueA = leading edge genes from sig
    pathways in DB-A that don't appear in any sig pathway of DB-B
  - genes_in_sig_pathways_uniqueB (symmetric)

Bootstrap CIs at the stratum (input_category × database_pair) level:
  - dataset-level resampling, B=1000

Outputs:
  e2/per_input_per_db_metrics.tsv
  e2/per_input_db_pair_metrics.tsv
  e2/per_database_within_metrics.tsv
  e2/per_stratum_bootstrap.tsv
  e2/E2_findings.md
"""
from __future__ import annotations
from itertools import combinations
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path("/mnt/nvme1/omics-audit/phase2/audit1_main")
RUNS = ROOT / "e2/runs"
INVENTORY = ROOT / "datasets/audit1_main_inputs.tsv"
DB_DIR = ROOT / "databases"
PER_INPUT_PER_DB = ROOT / "e2/per_input_per_db_metrics.tsv"
PER_INPUT_PAIR = ROOT / "e2/per_input_db_pair_metrics.tsv"
WITHIN_DB = ROOT / "e2/per_database_within_metrics.tsv"
PER_STRATUM = ROOT / "e2/per_stratum_bootstrap.tsv"
FINDINGS = ROOT / "e2/E2_findings.md"

DATABASES = ["hallmark", "c2_kegg", "c2_reactome", "c2_wikipathways", "c5_go_bp"]
B = 1000
RNG = np.random.default_rng(42)


def load(input_id: str, db: str) -> pd.DataFrame | None:
    p = RUNS / f"{input_id}__{db}.tsv"
    if not p.exists() or p.stat().st_size < 100:
        return None
    df = pd.read_csv(p, sep="\t")
    return df if len(df) > 0 else None


def sig_leading_edge_genes(df: pd.DataFrame, thr: float = 0.05) -> set[str]:
    sig = df[df["padj"].fillna(1) < thr]
    genes = set()
    for cell in sig["hit_genes"].fillna(""):
        if cell:
            genes.update(g.strip() for g in cell.split(";") if g.strip())
    return genes


def per_input_per_db(inv: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, ir in inv.iterrows():
        for db in DATABASES:
            d = load(ir["input_id"], db)
            if d is None:
                continue
            n_sig05 = int((d["padj"].fillna(1) < 0.05).sum())
            n_sig01 = int((d["padj"].fillna(1) < 0.01).sum())
            nes_abs = d.loc[d["padj"].fillna(1) < 0.05, "NES"].abs()
            rows.append({
                "input_id": ir["input_id"],
                "input_category": ir["input_category"],
                "database": db,
                "n_pathways_returned": len(d),
                "n_sig_padj_05": n_sig05,
                "n_sig_padj_01": n_sig01,
                "frac_sig_padj_05": round(n_sig05 / max(len(d), 1), 4),
                "median_abs_NES_sig": float(nes_abs.median()) if len(nes_abs) else float("nan"),
            })
    df = pd.DataFrame(rows)
    df.to_csv(PER_INPUT_PER_DB, sep="\t", index=False)
    print(f"per_input_per_db: {len(df)} rows")
    return df


def per_input_db_pair(inv: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for _, ir in inv.iterrows():
        loaded = {db: load(ir["input_id"], db) for db in DATABASES}
        sig_genes = {db: sig_leading_edge_genes(d) if d is not None else set()
                     for db, d in loaded.items()}
        for da, db in combinations(DATABASES, 2):
            if loaded[da] is None or loaded[db] is None:
                continue
            ga, gb = sig_genes[da], sig_genes[db]
            jac = (len(ga & gb) / len(ga | gb)) if (ga | gb) else float("nan")
            # Pathway-name overlap of sig pathways (rare across DBs)
            sig_a_names = set(loaded[da][loaded[da]["padj"].fillna(1) < 0.05]["pathway_name"])
            sig_b_names = set(loaded[db][loaded[db]["padj"].fillna(1) < 0.05]["pathway_name"])
            common_sig_names = sig_a_names & sig_b_names
            # Direction concordance on shared sig pathway names
            if common_sig_names:
                merged = (loaded[da].set_index("pathway_name").loc[list(common_sig_names), "NES"]
                          .to_frame("nes_a")
                          .join(loaded[db].set_index("pathway_name").loc[list(common_sig_names), "NES"]
                                .to_frame("nes_b")))
                dir_concord = float((np.sign(merged["nes_a"]) == np.sign(merged["nes_b"])).mean())
            else:
                dir_concord = float("nan")
            # Coverage gap
            uniq_a = ga - gb
            uniq_b = gb - ga
            rows.append({
                "input_id": ir["input_id"],
                "input_category": ir["input_category"],
                "db_a": da,
                "db_b": db,
                "leading_edge_gene_jaccard_sig05": round(jac, 4),
                "n_sig_pathway_names_shared": len(common_sig_names),
                "direction_concordance_on_shared": dir_concord,
                "n_genes_unique_to_a": len(uniq_a),
                "n_genes_unique_to_b": len(uniq_b),
                "n_sig_genes_a": len(ga),
                "n_sig_genes_b": len(gb),
            })
    df = pd.DataFrame(rows)
    df.to_csv(PER_INPUT_PAIR, sep="\t", index=False)
    print(f"per_input_db_pair: {len(df)} rows")
    return df


def within_database_redundancy() -> pd.DataFrame:
    """Mean pairwise gene-set Jaccard within each database (sampled to 5% of pairs)."""
    rows = []
    for db in DATABASES:
        # Use human DB for redundancy (mouse is structurally identical, slightly fewer pathways)
        p = DB_DIR / f"msigdb_{db}_human.tsv"
        if not p.exists():
            rows.append({"database": db, "n_pathways": 0,
                         "mean_pairwise_jaccard": float("nan"),
                         "mean_pathway_size": float("nan")})
            continue
        df = pd.read_csv(p, sep="\t")
        gs = {}
        for pwy, sub in df.groupby("pathway_name"):
            genes = set(sub["gene_symbol"].dropna().astype(str))
            if len(genes) >= 5:
                gs[pwy] = genes
        names = list(gs.keys())
        n = len(names)
        sizes = [len(gs[n_]) for n_ in names]
        # Sample pairs (cap at 1000 pairs for speed)
        all_pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        n_sample = min(1000, len(all_pairs))
        idx = RNG.choice(len(all_pairs), size=n_sample, replace=False) if len(all_pairs) > n_sample else range(len(all_pairs))
        jaccards = []
        for k in idx:
            i, j = all_pairs[k]
            a, b = gs[names[i]], gs[names[j]]
            u = a | b
            if u:
                jaccards.append(len(a & b) / len(u))
        rows.append({
            "database": db,
            "n_pathways": n,
            "mean_pathway_size": float(np.mean(sizes)),
            "median_pathway_size": float(np.median(sizes)),
            "mean_pairwise_jaccard_sample": float(np.mean(jaccards)) if jaccards else float("nan"),
            "n_pathway_pairs_sampled": len(jaccards),
        })
    df = pd.DataFrame(rows)
    df.to_csv(WITHIN_DB, sep="\t", index=False)
    print(f"within_database: {len(df)} rows")
    return df


def bootstrap_ci(values: np.ndarray, b: int = B) -> tuple[float, float, float]:
    vals = values[~np.isnan(values)]
    if len(vals) < 3:
        med = float(np.nanmedian(values)) if len(vals) > 0 else float("nan")
        return med, float("nan"), float("nan")
    n = len(vals)
    medians = np.empty(b)
    for i in range(b):
        medians[i] = np.nanmedian(vals[RNG.integers(0, n, n)])
    return float(np.median(vals)), float(np.percentile(medians, 2.5)), float(np.percentile(medians, 97.5))


def per_stratum_bootstrap(per_input_pair: pd.DataFrame) -> pd.DataFrame:
    rows = []
    metrics = ["leading_edge_gene_jaccard_sig05", "direction_concordance_on_shared"]
    for category, gcat in per_input_pair.groupby("input_category"):
        for (da, db), gpair in gcat.groupby(["db_a", "db_b"]):
            row = {"input_category": category, "db_a": da, "db_b": db,
                   "n_inputs": len(gpair)}
            for m in metrics:
                med, lo, hi = bootstrap_ci(gpair[m].values)
                row[f"{m}_median"] = med
                row[f"{m}_ci_lo"]  = lo
                row[f"{m}_ci_hi"]  = hi
            rows.append(row)
    df = pd.DataFrame(rows)
    df.to_csv(PER_STRATUM, sep="\t", index=False)
    print(f"per_stratum: {len(df)} rows (3 categories × 10 db pairs)")
    return df


def write_findings(per_input_per_db: pd.DataFrame,
                   per_input_pair: pd.DataFrame,
                   within_db: pd.DataFrame,
                   per_stratum: pd.DataFrame):
    lines = ["# E2 — Database choice effect (fgsea held constant)\n",
             "**Tool:** fgsea (GSEA-family)",
             "**Databases (5):** MSigDB Hallmark, C2 KEGG_LEGACY, C2 Reactome, C2 WikiPathways, C5 GO:BP",
             "**Bootstrap:** dataset-level, B=1000",
             "",
             "**Caveat:** E2 is GSEA-paradigm-specific. ORA tools may show",
             "different database sensitivity. If findings here are strong, an",
             "ORA-companion run could be warranted (see DEFERRED.md, decision",
             "deferred to CP7).",
             "",
             "## Within-database structure (redundancy + size)\n",
             "| Database | n pathways | mean pathway size | mean pairwise Jaccard (sampled) |",
             "|---|---:|---:|---:|"]
    for _, r in within_db.iterrows():
        lines.append(f"| {r['database']} | {int(r['n_pathways'])} | {r.get('mean_pathway_size', float('nan')):.1f} | {r.get('mean_pairwise_jaccard_sample', float('nan')):.4f} |")

    lines += ["",
              "## Per-database significant-pathway burden (medians across all 125 inputs)\n",
              "| Database | median n_sig (padj<0.05) | median frac_sig | median |NES| of sig |",
              "|---|---:|---:|---:|"]
    for db, g in per_input_per_db.groupby("database"):
        lines.append(f"| {db} | {g['n_sig_padj_05'].median():.0f} | {g['frac_sig_padj_05'].median():.3f} | {g['median_abs_NES_sig'].median():.3f} |")

    lines += ["",
              "## Cross-database leading-edge gene Jaccard (per stratum, medians + 95% bootstrap CI)\n",
              "| Category | DB-A | DB-B | n | sig-LE gene Jaccard | direction concordance (shared names) |",
              "|---|---|---|---:|---:|---:|"]

    def fmt(med, lo, hi):
        if np.isnan(med): return "—"
        if np.isnan(lo) or np.isnan(hi): return f"{med:.3f}"
        return f"{med:.3f} [{lo:.3f}, {hi:.3f}]"

    for _, r in per_stratum.iterrows():
        lines.append(
            f"| {r['input_category']} | {r['db_a']} | {r['db_b']} | {r['n_inputs']} | "
            f"{fmt(r['leading_edge_gene_jaccard_sig05_median'], r['leading_edge_gene_jaccard_sig05_ci_lo'], r['leading_edge_gene_jaccard_sig05_ci_hi'])} | "
            f"{fmt(r['direction_concordance_on_shared_median'], r['direction_concordance_on_shared_ci_lo'], r['direction_concordance_on_shared_ci_hi'])} |"
        )

    lines += ["",
              "## Interpretation",
              "",
              "**Leading-edge gene Jaccard** is the most meaningful cross-database",
              "metric here. Different databases test different pathway sets, so",
              "pathway-name overlap is rare. But the *genes that drive*",
              "significant pathways in each database can be compared — high",
              "Jaccard means the databases capture the same biology through",
              "different pathway lenses; low Jaccard means the databases",
              "fundamentally see different signals.",
              "",
              "**Direction concordance on shared pathway names** is included for",
              "completeness but is data-thin: only Hallmark and select C2 pathways",
              "have name overlap with other databases, so this metric has",
              "small N at the stratum level.",
              "",
              "**Within-database redundancy** (gene-set Jaccard between pathway",
              "pairs in the same database) shows how much each database",
              "double-counts the same biology under different pathway names.",
              "High redundancy means a 'significant pathway count' from that",
              "database is inflated by hierarchical / overlapping pathway",
              "definitions (notably C5 GO:BP).",
              "",
              "## Caveats",
              "",
              "1. **GSEA paradigm only.** ORA tools may show different database",
              "   sensitivity. ORA-companion deferred to CP7 decision per spec.",
              "2. **Reactome.db full extract not included here.** This audit uses",
              "   the MSigDB C2:CP:Reactome subset for cross-DB consistency",
              "   (msigdbr-curated, gene-mapped uniformly). reactome.db full",
              "   extract is available in CP2 cache for orthogonal validation.",
              "3. **WikiPathways MSigDB subset vs fresh GMT not contrasted here.**",
              "   The fresh May-2026 GMT is cached but not exercised by E2 since",
              "   comparing MSigDB-snapshot vs fresh-GMT is database freshness,",
              "   not database choice. Captured in DEFERRED.md if it becomes a",
              "   question for E1b or a separate sub-audit.",
              ""]

    FINDINGS.write_text("\n".join(lines) + "\n")
    print(f"Wrote {FINDINGS}")


def main():
    inv = pd.read_csv(INVENTORY, sep="\t")
    pid = per_input_per_db(inv)
    pip = per_input_db_pair(inv)
    wd  = within_database_redundancy()
    pst = per_stratum_bootstrap(pip)
    write_findings(pid, pip, wd, pst)
    sys.path.insert(0, "/mnt/nvme1/omics-audit/phase2/scripts")
    from dge_native import register_output
    for f in [PER_INPUT_PER_DB, PER_INPUT_PAIR, WITHIN_DB, PER_STRATUM, FINDINGS]:
        register_output(f, kind="audit1_main_e2_metric")
    print("Registered E2 metrics + findings.")


if __name__ == "__main__":
    main()
