# Audit 3 — scRNA-seq counting tools: synthesis

_Closed: 2026-05-23. Canonical drive `/mnt/nvme1/omics-audit/`. Lock: `phase2/repro.lock`._
_prior_audit_relationship: novel — first systematic audit of scRNA-seq counting-tool agreement in the corpus._

## The question

Phase 1 audits treated the count matrix as fixed input. In practice the matrix
depends on which tool counted the reads. Audit 3 closed that loop: do
STARsolo, alevin-fry, and kb-python (kallisto|bustools) — at default settings,
on identical 10x FASTQs — produce the same counts, the same cells, and the
same biology? (CellRanger gated → Audit 3c; 5' v2 gated → Audit 3e; multiome
gated → Audit 3d.)

Working set: **9 datasets**, 2 chemistries (3' v2 / v3), human + mouse, PBMC
+ T-cells + brain + lung + kidney + intestine. 4 tool configurations
(STARsolo default, STARsolo CR-mimic, alevin-fry, kb-python).

## The arc: counts agree → cells nest → biology diverges (on hard tissue)

**C1 — counts converge.** Cross-tool per-gene Spearman ~0.96 (9 datasets,
overlapping bootstrap CIs, chemistry-independent). Per-cell UMI ranking
agrees 0.93–0.97. *Counting-tool choice is low-risk for count-level analysis.*
(One methodological scar: a first-pass bug — not collapsing alevin-fry's
splici S/U/A buckets — produced a spurious "alevin-fry is the outlier" at
ρ~0.58; fixed by summing S+A. The corrected finding is convergence. This is
why §3.5 below was proposed.)

**C2 — cells form a permissiveness chain.** Native cell-callers nest
**STARsolo ⊂ alevin-fry ⊂ kb-python**, universal across 9/9 datasets
(containment 0.99–1.00). Negligible on clean data (cell-set Jaccard ~0.95+),
but up to **3× cell-count spread on high-ambient intestine** (22K / 42K / 67K),
the extra cells being low-UMI ambient barcodes. Crucially: applying **one
caller (EmptyDrops_CR) uniformly to every tool's raw counts** collapses the
divergence — mean Jaccard 0.88 (native) → **0.99 (uniform)**, 94% of
cross-tool cell-calling divergence removed, even on intestine (0.33 → 0.99).
*The disagreement is in the cell-calling algorithm, not the counts.*

**C3 — biology diverges on high-ambient tissue, washes out on clean tissue.**
Running a held-constant scanpy pipeline (QC → scDblFinder → scry-deviance HVG
→ PCA → Leiden → markers → CellTypist) per tool: on **intestine**, the four
native callers yield different clustering (ARI 0.61), different markers
(Jaccard 0.48), and different annotations (label agreement 0.34 STAR-vs-kb) —
*for the cells all tools agree are real*. On the **clean PBMC control** the
same comparison gives ARI 0.88 / markers 0.89 / labels 0.93. The contested
low-UMI cells survive standard QC (median 2,086 genes), so QC does not rescue
it. *On high-ambient tissue, counting-tool choice (via its cell-caller
default) cascades into biological conclusions.*

## The four rules (CP7, all hard_default under §5.3.2, all `novel`)

| rule | severity | what it tells the user |
|---|---|---|
| `scrna_counting_tool_per_gene_count_convergence` | info | reassurance: pick a counter on operational grounds; counts agree (~0.96) |
| `scrna_cell_calling_permissiveness_chain` | warn | awareness: native callers nest STAR⊂alevin⊂kb; benign on clean data, 3× on high-ambient (magnitude flag_and_warn, n=1) |
| `scrna_uniform_cell_caller_eliminates_disagreement` | warn | **action**: pin the caller (uniform EmptyDrops_CR), not the counter → Jaccard 0.99 |
| `scrna_cell_calling_biological_propagation_high_ambient` | warn | stakes: on high-ambient tissue, caller choice changes the biology you publish (existence hard_default; generalization flag_and_warn, n=1) |

One-line takeaway the rules encode: **counting-tool choice matters less than
discourse assumes; cell-calling algorithm + parameters matter more; the
reproducibility lever is pinning the caller with documented parameters,
especially on high-ambient tissue.**

## Honest limitations (carried into the rules per §3.1)

- **n=9 datasets** — below the §5.3.1 selection-audit floors; tiered under the
  new §5.3.2 equivalence criteria (≥8). The amendment exists *because* this
  audit surfaced the gap.
- **High-ambient propagation rests on n=1** (intestine) + n=1 clean control.
  Existence is solid (non-overlapping contrast CIs); generalization is
  flag_and_warn pending more high-ambient tissues.
- 5' v2, CellRanger, multiome excluded (gated) → Audits 3e / 3c / 3d.
- Whether intestine's contested cells are real rare types vs ambient artifact
  is unresolved.

## Standards contributions

- **§5.3.2 equivalence-finding tier criteria** (adopted) — tier boundaries for
  agreement/convergence audits, which §5.3.1 (tool-selection) could not express.
- **Clean-control-for-stress-test** closeout-amendment candidate (queued) —
  from C3's PBMC-vs-intestine design.
- **§3.5 "bootstrap CIs reflect sampling variance, not correctness"** candidate
  (queued for batched pass) — from C1's tight-CIs-on-buggy-data episode.
- 2 pipeline_step registry names: `scrnaseq_counting`, `scrnaseq_cell_calling`.

## Deferred follow-ups (see DEFERRED.md)

- **3c** CellRanger (needs license) · **3d** multiome (gated whitelist) ·
  **3e** 5' v2 (gated FASTQ) · **3f** Tabula Muris re-pull (SRA submission gap).
- **3g** (alevin-fry knee follow-up) — RESOLVED in-audit by CP5 Deliverable C.

## Provenance

All Audit 3 outputs hash-registered in `phase2/repro.lock` (91 entries,
verify 91/91, 0 drift at close). Superseded buggy CP4 outputs retained under
`c1/superseded_2026-05-18_buggy_usa_strip/` per §1.5. Native runtimes
throughout (scry / scDblFinder / EmptyDrops_CR via Rscript subprocess; no
rpy2). BioOrchestrator integration deferred to the next batched update.
