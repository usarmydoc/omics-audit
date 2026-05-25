# CP1 — Tool failure modes, anomalies, and methodological notes

_Deliverable A, 2026-05-24. 9 datasets × 3 tools. All runs exit 0; no aborts after the gene-ID-check fix._

## Gene-ID integrity (C1 USA-mode lesson, applied prospectively)
- **PASS.** All three tools emit the identical STARsolo Ensembl gene-ID universe (63,241 IDs, same version suffixes). CellBender barcode-match and gene-id-match to the raw matrix = 1.000 on every dataset. No convention mismatch (ENSG vs symbol vs USA-suffix). The comparison correlates like-for-like IDs.
- The compare script's first pass aborted on lung at "gene-id overlap 0.988 < 0.99" — this was a **false alarm in the check, not the data**: it compared *nonzero-expression* gene sets (which legitimately differ because CellBender uses its own called cells vs STARsolo filtered), not the ID universe. Fixed: the C1 check now compares the full gene-id universe (1.000 overlap everywhere); correlations run on genes nonzero in both tools.

## Contamination-magnitude sanity (CP4 §3.5 check)
Per-dataset global fraction removed, checked against published tissue expectations (PBMC 1–10%; higher for organoid/intestine):

| dataset | SoupX | CellBender | DecontX | within expectation? |
|---|---|---|---|---|
| PBMC (1k/5k/10k) | 0.010 | 0.011–0.016 | 0.026–0.078 | yes |
| neuron_1k | 0.010 | 0.012 | 0.040 | yes |
| t_3k_v2 / pbmc_4k_v2 | 0.085 / 0.056 | 0.018 | 0.029 / 0.102 | yes (v2 higher ambient) |
| lung organoid | 0.010 | 0.047 | 0.257 | plausible (organoid) |
| intestine (high-ambient) | **0.015** | **0.352** | 0.211 | CB/DecontX plausible for C3 high-ambient tissue; **SoupX = under-detection** |

- **No degenerate output** (no tool reports 0% or 100% across all genes). Degeneracy guard passed.
- The magnitudes are plausible per §3.5 **except SoupX's floor behavior**, which is a genuine tool limitation (below), not a config error.

## SoupX — floors at rho≈0.01, under-detects high ambient
- `autoEstCont` returned **global rho ≈ 0.01 (its soft floor) on 6/9 datasets**, including the high-ambient **intestine (1.5%)** where CellBender found 35% and DecontX 21%. Logs confirm genuine estimation (e.g. intestine: "Using 534 independent estimates of rho. Estimated global rho of 0.01"), not a crash or hard fallback.
- SoupX *did* estimate higher contamination on the two v2 datasets (8.5%, 5.6%) — so it detects *moderate* ambient but **fails to detect extreme ambient**. This is the single most consequential tool behavior in CP1 and drives the agreement collapse on intestine.
- Per-gene structure still varies even at the rho floor (via the soup profile), so Spearman remains meaningful; but SoupX's *magnitude* is not comparable to the other two on high-ambient data.

## CellBender — symbol column + low-expression zeroing
- CellBender's `_filtered.h5` var has **no gene-symbol column**; `anndata_from_h5` falls back to the Ensembl ID as "symbol." Mitigation: the compare script sources a **canonical gene_id→symbol map from SoupX/DecontX** (identical gene_ids) for all mito/symbol logic. Without this, CellBender-first pairs would have classified **zero** mito genes.
- CellBender **zeros ambient for many low-expression genes** (only ~34% of genes get >0 removed on clean data). Consequence: in the lowest expression deciles its contamination is constant → Spearman undefined (`spearman_lowexpr` = NaN for several CellBender pairs). This is real tool behavior, reported as NaN, not imputed.

## DecontX — most aggressive, weak mito agreement
- Consistently removes the most on clean/moderate data and is the relative outlier in rank agreement with both other tools (mean ρ ≈ 0.39–0.41).
- **Mito-gene agreement with the others is poor and sometimes negative** (kidney CellBender_vs_DecontX mito ρ = −0.35; lung −0.12) — DecontX ranks mitochondrial-gene contamination differently from SoupX/CellBender.

## Cell-set caveat (Deliverable A scope)
- SoupX/DecontX operate on STARsolo's filtered cells; CellBender on its own called cells (raw input). Per-gene contamination is a gene-level proportion (robust to moderate cell-set differences), but cross-tool comparisons involving CellBender mix slightly different cell populations. This is inherent to running each tool in its native mode and is documented rather than forced (cell-calling is Deliverable B's axis).
