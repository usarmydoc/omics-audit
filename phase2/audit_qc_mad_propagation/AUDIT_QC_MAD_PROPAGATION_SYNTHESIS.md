# Audit QC-MAD Propagation — Downstream Biological Consequences

_Phase 2. CP0–CP2 complete 2026-05-25. Standards: AUDIT_STANDARDS.md v1.0.3 + §5.3.2._

## Question
Does QC filtering method choice (quantile-data-driven, quantile-fixed-floors,
MAD k=3, MAD k=5) cascade to biological conclusions in downstream scRNA-seq
analysis? Mirrors Audit 3 C3 and Ambient Correction CP3 methodology — the third
propagation test in the corpus.

## Working set
3 datasets from QC-MAD (soma_joinid-matched, identical cells; Census 2025-11-08):
- **blood** (clean control, Census tissue ≈ PBMC in analyst terms)
- **liver** (QC-MAD second-worst cell-set Jaccard 0.817, C1 vs MAD3)
- **small_intestine** (QC-MAD worst Jaccard 0.774, C2 vs MAD5)
50,000 cells per dataset. C2 (fixed-floors) is the comparison reference.

## Headline finding
QC filtering method choice produces a **modest, tissue-INDEPENDENT** effect on
downstream biology. ARI vs C2 ranges **0.80–0.91** across all 3 datasets and 3
method comparisons; annotation agreement **88–98%**. The edge-case tissue
(small_intestine, where QC-MAD found the largest cell-set disagreement) is **no
more affected** than the clean control (blood). Bootstrap CIs ±0.005.

## Specific findings
1. ARI on cells called by both methods is 0.80–0.91 across all comparisons, with
   no tissue dependence (small_intestine 0.84–0.87 ≈ blood 0.80–0.88; liver
   highest at 0.88–0.91).
2. Annotation agreement 88–98%; cell-type label counts nearly identical across
   methods (blood 19 all; liver 10–11; small_intestine 23–24).
3. Where methods disagree on which cells to retain (MAD5 keeps 11,300 low-gene
   cells C2 removes on small_intestine; MAD3 +10,534; C1 +6,099), the shared
   core co-clusters at ARI 0.84 — borderline cells add at the periphery without
   reorganizing core analysis.

## Mechanism
QC method choice filters cells without reshaping counts. The cells removed by
stricter methods are mostly low-information borderline cells (low n_genes, low
total_counts). When retained by permissive methods, they cluster at the edges of
existing clusters or form small standalone clusters without affecting the core
biology. This contrasts cell-calling and ambient correction choices, which change
*what counts exist* (admitting different barcodes, removing different counts) and
thus reshape downstream clustering structurally.

## Corpus-level observation (cross-audit pattern)
Three propagation tests now exist:
| audit | technical choice | edge-case tissue | clean tissue |
|-------|------------------|------------------|--------------|
| Audit 3 C3 | cell-calling | ARI 0.50–0.70 (high-ambient intestine) | ~0.86–0.89 |
| Ambient Correction CP3 | ambient correction | ARI 0.50–0.70 (intestine) | 0.85–0.90 |
| **QC-MAD Propagation** | **QC filtering method** | **ARI 0.84–0.87 (small_intestine)** | **0.80–0.88** |

**Pattern:** scRNA technical choices that **RESHAPE COUNTS** (cell-calling,
ambient correction) propagate to biology on edge-case tissue. Choices that
**FILTER CELLS without reshaping counts** (QC method) do not propagate in the
same way — the effect is modest and tissue-independent.

This is a corpus-level finding that emerges from three audits, not any single
one. Future technical-choice-variance audits can use the distinction as a
framing heuristic: ask whether the choice reshapes counts or filters cells, then
expect propagation behavior accordingly. **Deliberately NOT encoded as a rule
yet** (3 audits is suggestive, not load-bearing); revisit after a 4th
propagation test confirms or breaks it.

## Heumos 2023 positioning
Relationship: **extends**. Heumos recommends MAD-based filtering (Germain 2020)
as more robust without testing biological consequences. This audit shows QC
method choice is **methodological hygiene rather than a substantive analytical
decision**: the methods produce slightly different cell sets (the QC-MAD finding)
but switching between them is unlikely to change biological conclusions (ARI
0.80–0.91). Extends Heumos by quantifying that the MAD-vs-default choice is
biologically low-stakes — unlike cell-calling/ambient correction.

## §5.3.2 tiering
Modest-effect existence finding: 3 datasets, multiple tissues, all ARI ≥0.80,
tight CIs, tissue-independent uniform behavior → **hard_default** for the
modest/uniform claim. The notable result (QC method does NOT show edge-case
amplification) is the headline; broad tissue-independence generalization is
flag_and_warn (n=3).

## Limitations
- 3 datasets (subset of QC-MAD's 8); replication would strengthen tissue-independence.
- C2 reference; alternative references give the same pairwise comparisons.
- High-ambient tissues beyond small_intestine not tested; "no propagation" might
  not generalize to extreme ambient burden.
- Default downstream parameters; clustering-resolution / HVG sensitivity untested.

## Prior audit relationship
**extends_prior.** Extends QC-MAD (tests whether its modest cell-set disagreement
propagates to biology) and extends Audit 3 C3 + Ambient CP3 (adds the third
propagation test, surfaces the cross-audit pattern).

## Audit corpus state after this audit
- 4 scRNA audits closed (Audit 3 counting, QC-MAD, Ambient Correction, QC-MAD Propagation).
- 3 propagation tests; pattern: counts-reshaping choices propagate on edge-case
  tissue, cell-filtering choices don't.
- 11 scRNA rules as reference documentation.
- §5.3.2 validated across 4 audits. BioOrchestrator tabled throughout.

## Rule drafted (CP2)
`draft_rules/scrna_qc_method_choice_modest_propagation.yaml` (hard_default, info;
validates under `--strict-steps`).
