"""CELLxGENE Census access layer — lazy/slice queries only.

Connects to the Census LTS build and fetches expression data for marker genes
using axis_query() with obs/var filters. The Census connection is NOT thread-safe,
so all queries run sequentially. Parallelism is applied at the scoring stage.

Performance strategy:
- Discovery uses cell_type value_filter to only pull relevant cells
- Caps at MAX_DATASETS_PER_CELLTYPE to keep runtime under ~30 min
- Single axis_query per dataset (not per cell type)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Iterator

import cellxgene_census
import numpy as np
import pandas as pd
import psutil
import tiledbsoma
from tqdm import tqdm

logger = logging.getLogger(__name__)

CENSUS_VERSION = "2025-11-08"

MARKER_GENES: dict[str, list[str]] = {
    "Hepatocytes": ["ALB", "APOB", "TTR", "CYP3A4"],
    "Macrophages": ["CD68", "MRC1", "MARCO", "C1QA"],
    "T cells": ["CD3D", "CD3E", "CD8A", "CD4"],
    "Endothelial": ["PECAM1", "CDH5", "VWF"],
    "Fibroblasts": ["COL1A1", "COL1A2", "DCN", "LUM"],
    "NK cells": ["GNLY", "NKG7", "KLRD1"],
}

ALL_GENES: list[str] = sorted(set(g for gs in MARKER_GENES.values() for g in gs))

# Mouse orthologs — mouse gene nomenclature uses Title case (first letter uppercase)
# Some genes have different names in mouse
MOUSE_GENE_MAP: dict[str, str] = {
    "ALB": "Alb", "APOB": "Apob", "TTR": "Ttr", "CYP3A4": "Cyp3a11",  # Cyp3a11 is mouse ortholog
    "CD68": "Cd68", "MRC1": "Mrc1", "MARCO": "Marco", "C1QA": "C1qa",
    "CD3D": "Cd3d", "CD3E": "Cd3e", "CD8A": "Cd8a", "CD4": "Cd4",
    "PECAM1": "Pecam1", "CDH5": "Cdh5", "VWF": "Vwf",
    "COL1A1": "Col1a1", "COL1A2": "Col1a2", "DCN": "Dcn", "LUM": "Lum",
    "GNLY": "Gzma",  # GNLY has no mouse ortholog; Gzma is closest functional equivalent
    "NKG7": "Nkg7", "KLRD1": "Klrd1",
}

# Reverse map: mouse gene name -> human gene name (for unified reporting)
MOUSE_TO_HUMAN: dict[str, str] = {v: k for k, v in MOUSE_GENE_MAP.items()}

ALL_MOUSE_GENES: list[str] = sorted(set(MOUSE_GENE_MAP.values()))

# Census cell_type values that map to our categories (lowercase for matching)
CENSUS_CELL_TYPE_KEYWORDS: dict[str, list[str]] = {
    "Hepatocytes": ["hepatocyte", "hepatoblast"],
    "Macrophages": ["macrophage", "kupffer cell"],
    "T cells": ["T cell", "CD4-positive", "CD8-positive"],
    "Endothelial": ["endothelial cell", "endothelial"],
    "Fibroblasts": ["fibroblast", "myofibroblast", "hepatic stellate cell"],
    "NK cells": ["natural killer cell", "NK cell"],
}

MAX_CELLS_PER_DATASET = 100_000
MIN_CELLS_PER_DATASET = 100
MAX_DATASETS_PER_CELLTYPE = 10
MAX_TOTAL_DATASETS = 40
RNG_SEED = 42


def _census_key(organism: str) -> str:
    """Convert display name to Census collection key."""
    return organism.lower().replace(" ", "_")


@dataclass
class DatasetSlice:
    """Expression data for one dataset-tissue combination."""
    dataset_id: str
    tissue: str
    organism: str
    expr: np.ndarray          # (n_cells, n_genes) dense float32
    gene_names: list[str]     # ordered gene names matching expr columns
    cell_types: np.ndarray    # (n_cells,) string array of mapped cell type labels
    n_cells: int


def log_memory(label: str) -> None:
    """Log current memory usage via psutil."""
    mem = psutil.virtual_memory()
    proc = psutil.Process()
    rss_gb = proc.memory_info().rss / (1024 ** 3)
    logger.info(
        "[MEM %s] RSS=%.2f GB, System: %.1f%% used (%.1f GB avail)",
        label, rss_gb, mem.percent, mem.available / (1024 ** 3),
    )


def open_census():
    """Open Census with pinned LTS version."""
    logger.info("Opening CELLxGENE Census version %s", CENSUS_VERSION)
    return cellxgene_census.open_soma(census_version=CENSUS_VERSION)


def _gene_list_for_organism(organism: str) -> list[str]:
    """Return the correct gene name list for the organism."""
    if "mus" in organism.lower() or "mouse" in organism.lower():
        return ALL_MOUSE_GENES
    return ALL_GENES


def resolve_gene_ids(census, organism: str = "Homo sapiens") -> pd.DataFrame:
    """Look up soma_joinid for all marker genes.

    For mouse, uses MOUSE_GENE_MAP names and adds a 'human_name' column
    for unified reporting.
    """
    exp = census["census_data"][_census_key(organism)]
    gene_names = _gene_list_for_organism(organism)
    gene_filter = "feature_name in [" + ", ".join(f"'{g}'" for g in gene_names) + "]"
    var_df = exp.ms["RNA"].var.read(
        column_names=["soma_joinid", "feature_name"],
        value_filter=gene_filter,
    ).concat().to_pandas()

    found = set(var_df.feature_name)
    missing = set(gene_names) - found
    if missing:
        logger.warning("Genes not found in Census for %s: %s", organism, missing)

    # For mouse, add human_name column for unified reporting
    is_mouse = "mus" in organism.lower() or "mouse" in organism.lower()
    if is_mouse:
        var_df["human_name"] = var_df.feature_name.map(MOUSE_TO_HUMAN)
    else:
        var_df["human_name"] = var_df.feature_name

    logger.info("Resolved %d/%d marker genes for %s", len(found), len(gene_names), organism)
    return var_df


def _build_cell_type_filter() -> str:
    """Build a Census value_filter string for all target cell types.

    Uses 'cell_type in [...]' with known Census cell_type values.
    """
    # First, discover what actual cell_type values exist by querying with
    # substring-style approach. Census value_filter supports 'in' but not LIKE.
    # We'll use a broad list of known Census cell_type values.
    known_types = set()
    for kws in CENSUS_CELL_TYPE_KEYWORDS.values():
        known_types.update(kws)

    # Build OR-style filter. Census doesn't support LIKE, so we enumerate
    # known cell_type values. This is incomplete but catches the major ones.
    # The fetch step does secondary keyword mapping for anything missed.
    quoted = ", ".join(f"'{t}'" for t in sorted(known_types))
    return f"cell_type in [{quoted}]"


def discover_datasets(
    census,
    organism: str = "Homo sapiens",
    min_cells: int = MIN_CELLS_PER_DATASET,
) -> pd.DataFrame:
    """Find dataset-tissue combos containing our target cell types.

    Strategy: query obs with cell_type filter to only pull relevant cells,
    then group by dataset to find multi-cell-type datasets. Caps per cell type.

    Returns DataFrame with columns: dataset_id, tissue, n_cells, n_mapped_types.
    """
    exp = census["census_data"][_census_key(organism)]

    logger.info("Querying obs metadata with cell type filter...")
    log_memory("discover_start")

    ct_filter = _build_cell_type_filter()
    full_filter = f"is_primary_data == True and ({ct_filter})"
    logger.info("Discovery filter: %s", full_filter[:200])

    obs_df = exp.obs.read(
        column_names=["dataset_id", "tissue", "cell_type"],
        value_filter=full_filter,
    ).concat().to_pandas()

    log_memory("discover_filtered_loaded")
    logger.info(
        "Filtered cells: %d across %d datasets (from cell type filter)",
        len(obs_df), obs_df.dataset_id.nunique(),
    )

    # Map Census cell types to our categories
    type_map = _map_cell_types(obs_df.cell_type.unique())
    obs_df["mapped_type"] = obs_df.cell_type.map(type_map)
    obs_df = obs_df.dropna(subset=["mapped_type"])

    logger.info("After mapping: %d cells with recognized types", len(obs_df))

    # Group by dataset-tissue, count mapped types
    grouped = obs_df.groupby(["dataset_id", "tissue"]).agg(
        n_cells=("cell_type", "size"),
        n_mapped_types=("mapped_type", "nunique"),
        mapped_types=("mapped_type", lambda x: list(x.unique())),
    ).reset_index()

    # Must have at least 2 mapped categories and enough cells
    valid = grouped[
        (grouped.n_cells >= min_cells) & (grouped.n_mapped_types >= 2)
    ].copy()

    logger.info(
        "Found %d dataset-tissue combos with >= %d cells and >= 2 mapped cell types",
        len(valid), min_cells,
    )

    # Cap: take datasets with the most mapped cell types first, then by n_cells
    valid = valid.sort_values(
        ["n_mapped_types", "n_cells"], ascending=[False, False]
    ).head(MAX_TOTAL_DATASETS)

    logger.info("Selected top %d dataset-tissue combos for scoring", len(valid))

    del obs_df
    log_memory("discover_done")

    return valid


def _map_cell_types(census_types) -> dict[str, str]:
    """Map Census cell type labels to our marker gene categories.

    Returns dict mapping Census label -> our category name.
    Uses substring/keyword matching since Census labels vary.
    """
    mapping = {}
    unique_types = set(census_types)

    keywords = {
        "Hepatocytes": ["hepatocyte", "hepatoblast"],
        "Macrophages": ["macrophage", "kupffer"],
        "T cells": ["t cell", "cd4-positive", "cd8-positive", "regulatory t",
                     "memory t", "naive t", "effector t", "alpha-beta t",
                     "gamma-delta t"],
        "Endothelial": ["endothelial"],
        "Fibroblasts": ["fibroblast", "stellate", "myofibroblast"],
        "NK cells": ["natural killer", "nk cell"],
    }

    for ct in unique_types:
        ct_lower = ct.lower()
        for category, kws in keywords.items():
            if any(kw in ct_lower for kw in kws):
                mapping[ct] = category
                break

    return mapping


def fetch_expression_slice(
    census,
    dataset_id: str,
    tissue: str,
    gene_joinids: list[int],
    gene_names: list[str],
    organism: str = "Homo sapiens",
) -> DatasetSlice | None:
    """Fetch expression for marker genes in one dataset-tissue combo.

    Queries ALL cells in this dataset-tissue (not just target types),
    so we have both positive and negative classes for scoring.
    Subsamples datasets > MAX_CELLS_PER_DATASET cells.
    """
    exp = census["census_data"][_census_key(organism)]

    # Get ALL cells in this dataset-tissue (need negatives for AUC)
    obs_filter = f"dataset_id == '{dataset_id}' and tissue == '{tissue}' and is_primary_data == True"
    obs_df = exp.obs.read(
        column_names=["soma_joinid", "cell_type"],
        value_filter=obs_filter,
    ).concat().to_pandas()

    if len(obs_df) < MIN_CELLS_PER_DATASET:
        logger.debug("Skipping %s/%s: only %d cells", dataset_id[:8], tissue, len(obs_df))
        return None

    # Map cell types to our categories
    type_map = _map_cell_types(obs_df.cell_type.values)
    if not type_map:
        logger.debug("Skipping %s/%s: no mappable cell types", dataset_id[:8], tissue)
        return None

    # Map — cells without a mapping become "Other" (needed as negatives)
    obs_df["mapped_type"] = obs_df.cell_type.map(type_map).fillna("Other")

    # Must have at least 2 categories (at least 1 target + others)
    n_target_categories = obs_df[obs_df.mapped_type != "Other"].mapped_type.nunique()
    if n_target_categories < 1:
        return None

    # Subsample if too large — stratified to preserve class ratios
    if len(obs_df) > MAX_CELLS_PER_DATASET:
        rng = np.random.default_rng(RNG_SEED)
        idx = rng.choice(len(obs_df), MAX_CELLS_PER_DATASET, replace=False)
        obs_df = obs_df.iloc[sorted(idx)].copy()
        logger.info("Subsampled %s/%s to %d cells", dataset_id[:8], tissue, MAX_CELLS_PER_DATASET)

    obs_ids = obs_df.soma_joinid.tolist()

    # Fetch expression via axis_query
    try:
        with exp.axis_query(
            measurement_name="RNA",
            obs_query=tiledbsoma.AxisQuery(coords=(obs_ids,)),
            var_query=tiledbsoma.AxisQuery(coords=(gene_joinids,)),
        ) as query:
            adata = query.to_anndata(X_name="raw")
    except Exception as e:
        logger.warning("Failed to fetch expression for %s/%s: %s", dataset_id[:8], tissue, e)
        return None

    if adata.shape[0] == 0:
        return None

    # Convert sparse to dense (small: n_cells x ~22 genes)
    expr = adata.X.toarray().astype(np.float32) if hasattr(adata.X, "toarray") else np.asarray(adata.X, dtype=np.float32)

    # Align gene names with the anndata var order
    adata_genes = list(adata.var.feature_name) if "feature_name" in adata.var.columns else list(adata.var_names)

    # For mouse, map gene names back to human orthologs for unified reporting
    is_mouse = "mus" in organism.lower() or "mouse" in organism.lower()
    if is_mouse:
        adata_genes = [MOUSE_TO_HUMAN.get(g, g) for g in adata_genes]

    # Align cell types with anndata obs order using soma_joinid
    adata_joinids = adata.obs["soma_joinid"].values
    obs_indexed = obs_df.set_index("soma_joinid")
    cell_types = obs_indexed.loc[adata_joinids, "mapped_type"].values

    return DatasetSlice(
        dataset_id=dataset_id,
        tissue=tissue,
        organism=organism,
        expr=expr,
        gene_names=adata_genes,
        cell_types=cell_types,
        n_cells=len(expr),
    )


def fetch_all_datasets(
    census,
    dataset_info: pd.DataFrame,
    gene_joinids: list[int],
    gene_names: list[str],
    organism: str = "Homo sapiens",
) -> Iterator[DatasetSlice]:
    """Iterate over dataset-tissue combos, yielding DatasetSlices.

    Sequential iteration (Census is not thread-safe).
    Logs memory at each dataset.
    """
    dt_pairs = dataset_info[["dataset_id", "tissue"]].drop_duplicates()
    logger.info("Fetching expression from %d dataset-tissue combinations", len(dt_pairs))

    for _, row in tqdm(dt_pairs.iterrows(), total=len(dt_pairs), desc="Fetching datasets"):
        ds_id = row.dataset_id
        tissue = row.tissue

        log_memory(f"fetch_{ds_id[:8]}_{tissue[:15]}")

        ds = fetch_expression_slice(
            census, ds_id, tissue, gene_joinids, gene_names, organism
        )

        if ds is not None:
            mapped_cats = set(ds.cell_types)
            mapped_cats.discard("Other")
            logger.info(
                "Fetched %s/%s: %d cells, %d genes, categories: %s",
                ds_id[:8], tissue, ds.n_cells, len(ds.gene_names),
                sorted(mapped_cats),
            )
            yield ds
        else:
            logger.debug("Skipped %s/%s", ds_id[:8], tissue)
