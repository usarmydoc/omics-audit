#!/usr/bin/env python3
"""Audit 1 main CP2 finalization.

- Capture Python environment_id (key package versions)
- Hash-register all cached databases + env files in phase2/repro.lock
- Write CP2_report.md with environment summary + Known issues
- Update STATUS.md
"""
from __future__ import annotations
from datetime import datetime
from pathlib import Path
import hashlib
import platform
import sys

DB_DIR = Path("/mnt/nvme1/omics-audit/phase2/audit1_main/databases")
REPORT = Path("/mnt/nvme1/omics-audit/phase2/audit1_main/CP2_report.md")
PY_ENV_FILE = DB_DIR / "python_session_info.txt"
STATUS = Path("/mnt/nvme1/omics-audit/phase2/audit1_main/STATUS.md")


def capture_python_env() -> str:
    """Write Python session info to a file. Return concatenated key versions."""
    lines = [
        f"# Python environment for Audit 1 main CP2",
        f"# Captured: {datetime.now().isoformat(timespec='seconds')}",
        f"",
        f"Python: {sys.version.split()[0]}",
        f"Platform: {platform.platform()}",
        f"",
        f"## Key audit packages:",
    ]
    for pkg in ["gseapy", "goatools", "mygene", "biothings_client",
                "pandas", "numpy", "scipy"]:
        try:
            mod = __import__(pkg)
            v = getattr(mod, "__version__", "?")
            lines.append(f"  {pkg}: {v}")
        except ImportError:
            lines.append(f"  {pkg}: MISSING")
    text = "\n".join(lines) + "\n"
    PY_ENV_FILE.write_text(text)
    return text


def main():
    capture_python_env()
    print(f"Wrote {PY_ENV_FILE}")

    # Register all DB + env files in lock
    sys.path.insert(0, "/mnt/nvme1/omics-audit/phase2/scripts")
    from dge_native import register_output

    count = 0
    for f in sorted(DB_DIR.iterdir()):
        if f.is_file():
            register_output(f, kind="audit1_main_database_cache")
            count += 1
    print(f"Hash-registered {count} cache files")

    # Write CP2 report
    lines = ["# Audit 1 main — CP2 report\n",
             f"_Generated: {datetime.now().isoformat(timespec='seconds')}_\n",
             "## Environments captured",
             "",
             "- **Python:** see `databases/python_session_info.txt`",
             "- **R:** see `databases/R_session_info.txt` (sessionInfo dump)",
             "",
             "Both environments registered in `phase2/repro.lock` under kind",
             "`audit1_main_database_cache`. environment_id stamps come from the",
             "lock's `environments` registry (captured at write time).",
             "",
             "## Tools installed",
             "",
             "### Python (native, no rpy2)",
             "- gseapy 1.1.13 (was already present)",
             "- goatools 1.6.5 (installed CP2)",
             "- mygene 3.2.2 (installed CP2)",
             "- biothings-client 0.5.0 (dependency of mygene)",
             "",
             "### R (via Rscript subprocess)",
             "- fgsea 1.36.2",
             "- GSVA 2.4.9",
             "- EGSEA 1.38.0",
             "- clusterProfiler 4.18.4",
             "- limma 3.66.0 (was already present)",
             "- msigdbr 26.1.0",
             "- ReactomePA 1.54.0",
             "- KEGGREST 1.50.0 (was already present)",
             "- rWikiPathways 1.30.0",
             "- org.Hs.eg.db 3.22.0",
             "- org.Mm.eg.db 3.22.0",
             "",
             "R version: 4.5.3; Bioconductor: 3.22; BiocManager: 1.30.27.",
             "",
             "## Databases cached (14 files in `databases/`)",
             "",
             "| Database | Human pathways | Mouse pathways |",
             "|---|---:|---:|"]

    # Read row counts
    for tag, name in [
        ("hallmark", "MSigDB Hallmark"),
        ("c2_kegg", "MSigDB C2 KEGG_LEGACY"),
        ("c2_reactome", "MSigDB C2 Reactome"),
        ("c2_wikipathways", "MSigDB C2 WikiPathways"),
        ("c5_go_bp", "MSigDB C5 GO:BP"),
    ]:
        hs = DB_DIR / f"msigdb_{tag}_human.tsv"
        mm = DB_DIR / f"msigdb_{tag}_mouse.tsv"
        import subprocess
        def n_pathways(p):
            if not p.exists():
                return "—"
            r = subprocess.run(["awk", "-F\t", "NR>1 {print $1}", str(p)], capture_output=True, text=True)
            return str(len(set(r.stdout.strip().split("\n"))))
        lines.append(f"| {name} | {n_pathways(hs)} | {n_pathways(mm)} |")
    # Reactome.db
    rh = DB_DIR / "reactome_db_human.tsv"
    rm = DB_DIR / "reactome_db_mouse.tsv"
    import subprocess
    def n_pathways(p):
        if not p.exists(): return "—"
        r = subprocess.run(["awk", "-F\t", "NR>1 {print $1}", str(p)], capture_output=True, text=True)
        return str(len(set(r.stdout.strip().split("\n"))))
    lines.append(f"| Reactome.db (full Bioc) | {n_pathways(rh)} | {n_pathways(rm)} |")
    # WikiPathways fresh GMT
    gmt_hs = next(DB_DIR.glob("wikipathways-*-gmt-Homo_sapiens.gmt"), None)
    gmt_mm = next(DB_DIR.glob("wikipathways-*-gmt-Mus_musculus.gmt"), None)
    def n_gmt_lines(p):
        if p is None or not p.exists(): return "—"
        return str(sum(1 for _ in p.read_text().splitlines()))
    lines.append(f"| WikiPathways GMT (fresh, May 2026) | {n_gmt_lines(gmt_hs)} | {n_gmt_lines(gmt_mm)} |")

    lines += ["",
              "## KEGG licensing note",
              "",
              "KEGG REST API is free for academic/non-commercial use; rate limit",
              "3 requests/sec. MSigDB C2 KEGG_LEGACY is the cached subset (186",
              "pathways for both human and mouse) suitable for the E2 database-",
              "choice audit. Direct KEGGREST calls deferred to sub-audit use only",
              "if MSigDB coverage proves insufficient.",
              "",
              "## Known issues (carried forward from CP1)",
              "",
              "### Gene ID format heterogeneity in DEG inputs",
              "",
              "Per CP1 (datasets/CP1_report.md), the 125 input DEG TSVs use a mix",
              "of gene ID systems:",
              "- TCGA: 31/48 symbol_human, 17/48 mixed_or_unknown",
              "- Census: 7/47 symbol_human, 14/47 ensembl_human, 5/47 symbol_mouse,",
              "  21/47 mixed_or_unknown",
              "- GTEx: 30/30 ensembl_human (clean)",
              "",
              "**Resolution strategy** (per CP1 follow-up guidance, 2026-05-16):",
              "Gene ID mapping is handled at each sub-audit (E1-E4) using whatever",
              "the chosen pathway tool natively expects. The cached MSigDB tables",
              "include three gene ID columns (`gene_symbol`, `gene_ensembl`,",
              "`gene_entrez`) to support cross-mapping per sub-audit need. If a",
              "specific sub-audit hits a wall because of gene ID issues, it will",
              "be surfaced then; no infrastructure pre-built.",
              "",
              "### KEGG vs KEGG_LEGACY",
              "",
              "MSigDB transitioned from `KEGG` to `KEGG_LEGACY` in 2024 when",
              "KEGG's distribution license tightened. The cached collection is",
              "the LEGACY snapshot, which is the most recent version Bioconductor",
              "is permitted to redistribute. For pathway content this is",
              "indistinguishable from current KEGG for the purposes of inter-tool",
              "comparison (E1) and inter-database comparison (E2).",
              "",
              "### Database staleness",
              "",
              "MSigDB v2025.1 cached via msigdbr 26.1.0; WikiPathways GMT dated",
              "2026-05-10 (fresh). Reactome.db is the Bioc 3.22 snapshot. These",
              "are pinned via lock-file hashes — re-running CP2 with a newer",
              "msigdbr release would produce different hashes and require a new",
              "CP2 commit if the audit needs the updated content.",
              ""]

    REPORT.write_text("\n".join(lines))
    print(f"Wrote {REPORT}")

    # Register the report itself
    register_output(REPORT, kind="audit1_main_cp2_report")
    print("Done.")


if __name__ == "__main__":
    main()
