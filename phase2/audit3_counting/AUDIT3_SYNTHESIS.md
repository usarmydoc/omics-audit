# Audit 3 — scRNA-seq Counting Tools

_Closed: 2026-05-23. Canonical drive `/mnt/nvme1/omics-audit/`. Lock: `phase2/repro.lock`._
_prior_audit_relationship: novel — first systematic audit of scRNA-seq counting-tool agreement in the corpus._

## Question

Do scRNA-seq counting tools (STARsolo default, STARsolo CR-mimic, alevin-fry,
kb_count) produce equivalent outputs at default settings, and does any
disagreement propagate to biological conclusions?

## Working set

9 datasets, 2 chemistries (3' v2, 3' v3), human + mouse, tissues PBMC /
T-cells / neuron / lung / kidney / intestine. 4 tool configurations.
CellRanger excluded (gated → Audit 3c). 5' v2 dropped (URLs gated → 3e).
Multiome dropped (whitelist gated → 3d). Tabula Muris dropped (SRA submission
gap → 3f).

## Finding 1: Counts converge

Per-gene Spearman ~0.96 across all tool pairs, 9 datasets, both chemistries;
dataset-level bootstrap CIs overlap; chemistry-independent. Per-cell UMI
ranking agrees 0.93–0.97. Counting-tool choice is low-risk for count-level
analysis. (A first-pass bug — not collapsing alevin-fry's splici S/U/A buckets
— produced a spurious ρ~0.58 "outlier"; fixed by summing S+A. Corrected
finding is convergence.) **Rule:**
`scrna_counting_tool_per_gene_count_convergence` (hard_default, info).
Detail: `c1/C1_findings.md`.

## Finding 2: Cell-callers form a permissiveness chain

Native cell-callers nest universally: **STARsolo ⊂ alevin-fry ⊂ kb_count
bustools-knee** (9/9 datasets, containment 0.99–1.00). Negligible on clean
data (~0.95+ cell-set Jaccard); on high-ambient data a 3× cell-count spread
(intestine: 22K / 42K / 67K). Sanity check: kb's extra barcodes are low-UMI
ambient (median 3,414 vs 23,582 for jointly-called; 100% below the 10th
percentile of jointly-called cells). Companion: per-cell UMI Spearman 0.999 on
jointly-called cells (divergence confined to the calling boundary). **Rule:**
`scrna_cell_calling_permissiveness_chain` (hard_default on direction,
flag_and_warn on magnitude). Detail: `c2/C2_findings.md`.

## Finding 3: Uniform caller eliminates divergence

Applying EmptyDrops_CR uniformly to all 4 tools' raw counts collapses
cross-tool cell-set Jaccard from ~0.88 mean (native) to ~0.99 mean (uniform)
— a 94% reduction in divergence — including intestine (0.326 native → 0.987
uniform). Because the caller runs on each tool's raw counts, this isolates the
divergence to the calling algorithm, not the counts (which already agree,
Finding 1). **Rule:** `scrna_uniform_cell_caller_eliminates_disagreement`
(hard_default, warn) — the actionable recommendation. Detail:
`c2/C2_findings.md`, `c2/b_common_caller_4tool/`.

## Finding 4: Biological propagation on high-ambient tissue

Native cell-calling differences cascade to clustering, marker genes, and
cell-type annotation on high-ambient tissue. Intestine: ARI 0.60–0.65, marker
Jaccard 0.47–0.49, annotation agreement 0.34 (STAR vs kb). Clean PBMC control:
ARI 0.86–0.89, marker Jaccard 0.87–0.91, annotation agreement 0.90–0.97.
21,476 contested cells survive QC on intestine (median 2,086 genes/cell) vs 79
on PBMC. The clean-control contrast (non-overlapping bootstrap ARI CIs)
establishes ambient burden as the driver. **Rule:**
`scrna_cell_calling_biological_propagation_high_ambient` (hard_default on
existence, flag_and_warn on generalization). Detail: `c3/C3_findings.md`.

## Synthesis

The arc: counting tools agree on counts (Finding 1) → cell-calling
implementations diverge in permissiveness (Finding 2) → that divergence
propagates to biological conclusions on high-ambient tissue (Finding 4) →
a uniform downstream caller eliminates the divergence (Finding 3, the fix).

**Counting-tool choice matters less than current bioinformatics discourse
suggests. Cell-calling algorithm choice matters more than discourse suggests,
and the magnitude is tissue-dependent. The path to reproducibility is pinning
the downstream cell-caller (EmptyDrops_CR with documented parameters), not
pinning the counting tool.**

## Methodological observations queued

- **§3.5 candidate** — "bootstrap CIs reflect sampling variance, not
  correctness" (from the CP4 USA-mode bug: tight CIs on misconfigured data).
  Disposition pending (CP8 Step 5; see PENDING_AMENDMENTS.md).
- **§5.3.2** — equivalence-finding tier criteria. ADOPTED (this audit forced it).
- **Clean-control-for-stress-test** — the CP6 PBMC-alongside-intestine pattern,
  as a closeout-amendment candidate. Queued.
- **Intermediate-file retention** — CP5 needed alevin-fry rad files that had
  been cleaned; future audits should retain (or cheaply regenerate) the
  intermediates a downstream checkpoint may depend on.
- **Normalized vs raw for expression tests** — CP6 surfaced that audit prompts
  should state normalized-matrix vs raw-count explicitly for Wilcoxon-type
  tests (we used log-normalized, documented).

## Limitations

- 3 chemistries collapsed to 2 (5' v2 URLs gated).
- 3 binaries / 4 configs (CellRanger gated → Audit 3c).
- Mouse n=3 (Tabula Muris failed: SRA submission gap).
- High-ambient propagation rests on **n=1 tissue** (intestine); the clean
  control establishes the effect exists, but generalization needs replication.
- 9 datasets total; PBMC-dominant (5/9).

## Audits queued in DEFERRED.md

- **Audit 3c** — CellRanger re-test on existing datasets (bounded). Trigger: 10x license.
- **Audit 3d** — multiome chemistry (bounded). Trigger: gated 737K-arc-v1 whitelist access.
- **Audit 3e** — 5' v2 chemistry (bounded). Trigger: gated FASTQ URL access.
- **Audit 3f** — Tabula Muris re-pull from alternative source (bounded). Trigger: alt source with CB+UMI.
- **Audit 3g** — RESOLVED in CP5 Deliverable C (alevin-fry native knee cell-calling); not a future audit.
