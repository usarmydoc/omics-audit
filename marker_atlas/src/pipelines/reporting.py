"""Report generation — TSV atlas and markdown summary."""

from __future__ import annotations

import hashlib
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from .scoring import AggregatedScore, MarkerScore

logger = logging.getLogger(__name__)


def _sha256(path: Path) -> str:
    """Compute SHA-256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(
    aggregated: list[AggregatedScore],
    output_dir: Path,
) -> Path:
    """Write marker_reliability_atlas.tsv.

    Columns: gene | cell_type | organism | tissue | n_datasets | mean_AUC | sd_AUC |
             ci_lower | ci_upper | mw_pval | reliability_score
    """
    records = []
    for s in aggregated:
        records.append({
            "gene": s.gene,
            "cell_type": s.cell_type,
            "organism": s.organism,
            "tissue": s.tissue,
            "n_datasets": s.n_datasets,
            "mean_AUC": round(s.mean_auc, 4),
            "sd_AUC": round(s.sd_auc, 4),
            "ci_lower": round(s.ci_lower, 4),
            "ci_upper": round(s.ci_upper, 4),
            "mw_pval": f"{s.mw_pval:.2e}",
            "reliability_score": round(s.reliability_score, 4),
        })

    df = pd.DataFrame(records)
    path = output_dir / "marker_reliability_atlas.tsv"
    df.to_csv(path, sep="\t", index=False)
    logger.info("Wrote %s (%d rows)", path, len(df))
    return path


def write_summary_markdown(
    aggregated: list[AggregatedScore],
    per_dataset_scores: list[MarkerScore],
    run_metadata: dict,
    output_dir: Path,
) -> Path:
    """Write summary_report.md with tables and per-cell-type breakdowns."""
    path = output_dir / "summary_report.md"

    high = [s for s in aggregated if s.mean_auc >= 0.80]
    moderate = [s for s in aggregated if 0.65 <= s.mean_auc < 0.80]
    low = [s for s in aggregated if s.mean_auc < 0.65]

    lines = [
        "# Cell Type Marker Reliability Atlas — Summary Report",
        "",
        "## Run Metadata",
        "",
        f"- **Date:** {run_metadata.get('date', datetime.now(timezone.utc).isoformat())}",
        f"- **Census version:** {run_metadata.get('census_version', 'unknown')}",
        f"- **Organisms:** {', '.join(run_metadata.get('organisms', ['Homo sapiens']))}",
        f"- **Datasets queried:** {run_metadata.get('n_datasets', '?')}",
        f"- **Datasets with data:** {run_metadata.get('n_datasets_with_data', '?')}",
        f"- **Total cells scored:** {run_metadata.get('n_cells_total', '?'):,}",
        f"- **Marker genes:** {run_metadata.get('n_genes', len(set(s.gene for s in aggregated)))}",
        f"- **Cell types:** {run_metadata.get('n_cell_types', len(set(s.cell_type for s in aggregated)))}",
        "",
        "## Reliability Tiers",
        "",
        f"- **High (AUC >= 0.80):** {len(high)} gene-cell_type pairs",
        f"- **Moderate (0.65 <= AUC < 0.80):** {len(moderate)} gene-cell_type pairs",
        f"- **Low (AUC < 0.65):** {len(low)} gene-cell_type pairs",
        "",
    ]

    # High reliability table
    if high:
        lines.extend([
            "## High-Reliability Markers (AUC >= 0.80)",
            "",
            "| Gene | Cell Type | Organism | mean AUC | sd AUC | 95% CI | N datasets | Reliability |",
            "|------|-----------|----------|----------|--------|--------|------------|-------------|",
        ])
        for s in sorted(high, key=lambda x: x.mean_auc, reverse=True):
            lines.append(
                f"| {s.gene} | {s.cell_type} | {s.organism} | {s.mean_auc:.3f} | {s.sd_auc:.3f} | "
                f"[{s.ci_lower:.3f}, {s.ci_upper:.3f}] | {s.n_datasets} | {s.reliability_score:.3f} |"
            )
        lines.append("")

    # Low reliability table
    if low:
        lines.extend([
            "## Low-Reliability Markers (AUC < 0.65)",
            "",
            "| Gene | Cell Type | Organism | mean AUC | sd AUC | 95% CI | N datasets | Reliability |",
            "|------|-----------|----------|----------|--------|--------|------------|-------------|",
        ])
        for s in sorted(low, key=lambda x: x.mean_auc):
            lines.append(
                f"| {s.gene} | {s.cell_type} | {s.organism} | {s.mean_auc:.3f} | {s.sd_auc:.3f} | "
                f"[{s.ci_lower:.3f}, {s.ci_upper:.3f}] | {s.n_datasets} | {s.reliability_score:.3f} |"
            )
        lines.append("")

    # Per-cell-type breakdown by organism
    lines.extend(["## Per-Cell-Type Summary", ""])
    from .census_access import MARKER_GENES
    organisms = sorted(set(s.organism for s in aggregated))
    for org in organisms:
        lines.append(f"### {org}")
        lines.append("")
        for ct, genes in MARKER_GENES.items():
            ct_scores = [s for s in aggregated if s.cell_type == ct and s.organism == org]
            if not ct_scores:
                continue
            lines.append(f"#### {ct}")
            lines.append("")
            lines.append(f"Markers tested: {', '.join(genes)}")
            lines.append("")
            best = max(ct_scores, key=lambda x: x.mean_auc)
            worst = min(ct_scores, key=lambda x: x.mean_auc)
            lines.append(f"- **Best:** {best.gene} (AUC={best.mean_auc:.3f})")
            lines.append(f"- **Worst:** {worst.gene} (AUC={worst.mean_auc:.3f})")
            lines.append("")

    # Tissue-specific divergence
    lines.extend(["## Tissue-Specific Divergence", ""])
    lines.append("Markers with sd_AUC > 0.10 across datasets (inconsistent performance):")
    lines.append("")
    divergent = [s for s in aggregated if s.sd_auc > 0.10]
    if divergent:
        for s in sorted(divergent, key=lambda x: x.sd_auc, reverse=True):
            lines.append(f"- **{s.gene}** for {s.cell_type} ({s.organism}): sd_AUC={s.sd_auc:.3f} across {s.n_datasets} datasets")
    else:
        lines.append("No markers showed significant tissue-specific divergence.")
    lines.append("")

    # Output file hashes
    lines.extend(["## Output File Checksums", ""])
    tsv_path = output_dir / "marker_reliability_atlas.tsv"
    if tsv_path.exists():
        lines.append(f"- `marker_reliability_atlas.tsv`: `{_sha256(tsv_path)}`")
    lock_path = output_dir / "repro.lock"
    if lock_path.exists():
        lines.append(f"- `repro.lock`: `{_sha256(lock_path)}`")
    lines.append("")

    content = "\n".join(lines)
    path.write_text(content)

    # Add summary_report.md hash
    lines.append(f"- `summary_report.md`: `{_sha256(path)}`")

    logger.info("Wrote %s", path)
    return path
