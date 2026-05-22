# Audit 3 C1 findings — per-gene count agreement

_Generated: 2026-05-18T19:51:04_
_B = 1000 (dataset-level), stratified by chemistry × tool_pair_

## Headlines

_Synthesis written after inspecting per_dataset_metrics + per_stratum_bootstrap_

## Per-stratum bootstrap (rho_gene_total point + 95% CI)

| chemistry | tool_a | tool_b | n_datasets | point | 95% CI |
|---|---|---|---|---|---|
| 3p_v3 | star_default | star_cr_mimic | 1 | 0.987 | [0.987, 0.987] (n<3, wide CI) |
| 3p_v3 | star_default | alevin_fry | 1 | 0.945 | [0.945, 0.945] (n<3, wide CI) |
| 3p_v3 | star_default | kb_count | 1 | 0.949 | [0.949, 0.949] (n<3, wide CI) |
| 3p_v3 | star_cr_mimic | alevin_fry | 1 | 0.952 | [0.952, 0.952] (n<3, wide CI) |
| 3p_v3 | star_cr_mimic | kb_count | 1 | 0.950 | [0.950, 0.950] (n<3, wide CI) |
| 3p_v3 | alevin_fry | kb_count | 1 | 0.949 | [0.949, 0.949] (n<3, wide CI) |

## Per-stratum bootstrap (rho_cell_total_umi)

| chemistry | tool_a | tool_b | n_datasets | point | 95% CI |
|---|---|---|---|---|---|
| 3p_v3 | star_default | star_cr_mimic | 1 | 0.939 | [0.939, 0.939] (n<3, wide CI) |
| 3p_v3 | star_default | alevin_fry | 1 | 0.871 | [0.871, 0.871] (n<3, wide CI) |
| 3p_v3 | star_default | kb_count | 1 | 0.933 | [0.933, 0.933] (n<3, wide CI) |
| 3p_v3 | star_cr_mimic | alevin_fry | 1 | 0.862 | [0.862, 0.862] (n<3, wide CI) |
| 3p_v3 | star_cr_mimic | kb_count | 1 | 0.965 | [0.965, 0.965] (n<3, wide CI) |
| 3p_v3 | alevin_fry | kb_count | 1 | 0.895 | [0.895, 0.895] (n<3, wide CI) |

## By-gene-category log2 ratio summary (across all 9 datasets)

| tool_a | tool_b | category | n_datasets | med(log2_median) | med(p05) | med(p95) |
|---|---|---|---|---|---|---|
| alevin_fry | kb_count | mitochondrial | 1 | nan | nan | nan |
| alevin_fry | kb_count | other | 1 | 0.000 | -0.556 | 0.130 |
| alevin_fry | kb_count | overlapping | 1 | 0.000 | -1.000 | 0.023 |
| alevin_fry | kb_count | pseudogene | 1 | 0.000 | -1.524 | 0.000 |
| star_cr_mimic | alevin_fry | mitochondrial | 1 | nan | nan | nan |
| star_cr_mimic | alevin_fry | other | 1 | 0.000 | -0.485 | 0.170 |
| star_cr_mimic | alevin_fry | overlapping | 1 | 0.000 | -0.322 | 0.415 |
| star_cr_mimic | alevin_fry | pseudogene | 1 | 0.000 | 0.000 | 1.000 |
| star_cr_mimic | kb_count | mitochondrial | 1 | nan | nan | nan |
| star_cr_mimic | kb_count | other | 1 | 0.000 | -0.585 | 0.034 |
| star_cr_mimic | kb_count | overlapping | 1 | 0.000 | -0.848 | 0.071 |
| star_cr_mimic | kb_count | pseudogene | 1 | 0.000 | -1.000 | 0.000 |
| star_default | alevin_fry | mitochondrial | 1 | nan | nan | nan |
| star_default | alevin_fry | other | 1 | 0.000 | -0.585 | 0.126 |
| star_default | alevin_fry | overlapping | 1 | 0.000 | -0.348 | 0.415 |
| star_default | alevin_fry | pseudogene | 1 | 0.000 | 0.000 | 1.263 |
| star_default | kb_count | mitochondrial | 1 | nan | nan | nan |
| star_default | kb_count | other | 1 | 0.000 | -0.585 | 0.000 |
| star_default | kb_count | overlapping | 1 | 0.000 | -0.807 | 0.014 |
| star_default | kb_count | pseudogene | 1 | 0.000 | -1.000 | 1.000 |
| star_default | star_cr_mimic | mitochondrial | 1 | nan | nan | nan |
| star_default | star_cr_mimic | other | 1 | 0.000 | -0.077 | 0.000 |
| star_default | star_cr_mimic | overlapping | 1 | 0.000 | -0.082 | 0.000 |
| star_default | star_cr_mimic | pseudogene | 1 | 0.000 | 0.000 | 1.000 |

## Sample size limitations

- 3p_v2: 0 datasets (below n=3 → bootstrap CIs are wide; interpret per-chemistry comparisons with this caveat).
- 3p_v3: 1 datasets.
- Mouse n=3 / Human n=6 — species stratification underpowered for mouse.
- Tissue: PBMC n=5 (incl. T cells), non-PBMC n=4.
- Per AUDIT_STANDARDS §3.1, these are headline findings, not caveats.
