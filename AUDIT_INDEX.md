# Audit Index

Registry of bioinformatics audits in this workspace. Each audit tests whether
conventional defaults hold under tool/parameter variation, on real data, with
hash-registered provenance (see `AUDIT_STANDARDS.md`). Rules drafted by audits
feed BioOrchestrator via batched updates.

Canonical drive: `/mnt/nvme1/omics-audit/`. Lock: `phase2/repro.lock`.

| Audit | Modality / question | Status | Rules | Location |
|---|---|---|---|---|
| **Audit 1 — pathway enrichment** | bulk + scRNA DEG → pathway tools (ORA/GSEA), databases, MT-correction, backgrounds | COMPLETE (2026-05-16) | 5 (+A6 muscat from 2a) | `phase2/audit1_main/`, rules in `phase2/draft_rules/` |
| **Phase 2 / 2a — version & method sensitivity** | DESeq2/edgeR/limma version drift; clustering metric; mito threshold; muscat dream | COMPLETE | (folded into Audit 1 rule set) | `phase2/version_sensitivity/`, `phase2a/` |
| **Audit 3 — scRNA-seq counting tools** | STARsolo / alevin-fry / kb-python: counts, cell-calling, biological propagation | **COMPLETE (2026-05-23)** | **4** (`audit3_counting/draft_rules/`) | `phase2/audit3_counting/` |
| **Audit QC-MAD — low-quality cell filtering** | MAD vs quantile vs fixed-floor QC filtering on scRNA-seq | **COMPLETE (2026-05-23)** | **2** (`audit_qc_mad/draft_rules/`) | `phase2/audit_qc_mad/` |
| **Audit Ambient Correction — ambient RNA removal** | SoupX / CellBender / DecontX: per-gene contamination equivalence, correction-vs-QC ordering, biological propagation | **COMPLETE (2026-05-25)** | **4** (`audit_ambient_correction/draft_rules/`) | `phase2/audit_ambient_correction/` |
| **Audit QC-MAD Propagation — QC method → biology** | does QC filtering-method choice (quantile / fixed-floor / MAD) cascade to downstream clustering/annotation? | **COMPLETE (2026-05-25)** | **1** (`audit_qc_mad_propagation/draft_rules/`) | `phase2/audit_qc_mad_propagation/` |

## Audit QC-MAD Propagation — at a glance

**Audit ID:** Audit QC-MAD Propagation (does QC filtering-method choice propagate to biology?). **Status:** COMPLETE / CLOSED.
**Work window:** 2026-05-25 (extends Audit QC-MAD + Audit 3 C3 + Ambient CP3).
**Checkpoints:** CP0 (inventory + Census re-pull feasibility) · CP1 (3 datasets ×
4 QC methods × downstream pipeline, 12 runs) · CP2 (synthesis + 1 rule + closeout)
— **all complete.**
**Standards contributed:** none (used existing §5.3.2 + `scrnaseq_qc_filtering` step).

**Key finding:** QC filtering-method choice is a **modest, tissue-INDEPENDENT**
downstream effect — ARI vs C2 reference 0.80–0.91 across blood/liver/small_intestine
(annotation 88–98%); the edge-case tissue is no worse than the clean control. The
edge-case amplification seen for cell-calling (C3) and ambient correction (CP3)
does NOT replicate — QC method is operational, not analytical.

**Cross-audit pattern (surfaced, not yet a rule):** technical choices that
**reshape counts** (cell-calling, ambient correction) propagate on edge-case
tissue; choices that **filter cells** (QC method) don't. Emerges from 3 propagation
tests; revisit as a rule after a 4th.

**1 rule:** `scrna_qc_method_choice_modest_propagation` (hard_default / info).
**Synthesis:** `phase2/audit_qc_mad_propagation/AUDIT_QC_MAD_PROPAGATION_SYNTHESIS.md`. **Heumos:** extends.
**Deferred follow-ups:** broader replication, extreme high-ambient, downstream-param sensitivity — `phase2/DEFERRED.md`.

## Audit Ambient Correction — at a glance

**Audit ID:** Audit Ambient Correction (scRNA-seq ambient RNA correction). **Status:** COMPLETE / CLOSED.
**Work window:** 2026-05-24 → 2026-05-25 (extends Audit 3 C3).
**Checkpoints:** CP0 (inventory + tool install; CellBender checkpoint patch) ·
CP1 (per-gene contamination, 9 datasets × 3 tools) · CP2 (ordering analysis,
2 orderings) · CP3 (biological propagation, intestine + PBMC) · CP4 (synthesis
+ 4 rules + closeout) — **all complete.**
**Standards contributed:** none to AUDIT_STANDARDS.md (used existing §5.3.2);
one registry extension `scrnaseq_ambient_correction` (no schema change).

**Arc:** tools produce non-equivalent per-gene corrections (CP1, ρ 0.39–0.57;
SoupX under-detects high ambient, CellBender tracks burden, DecontX most
aggressive) → only DecontX is ordering-sensitive while correct→QC is
universally stricter (CP2) → these differences propagate to clustering/markers/
annotation on high-ambient intestine (ARI 0.50–0.70 vs baseline) and wash out
on clean PBMC (0.85–0.90) (CP3).

**4 rules** (`audit_ambient_correction/draft_rules/`, all `extends_prior`, all
validate under `--strict-steps`):
1. `scrna_ambient_correction_tool_non_equivalence` (hard_default / info)
2. `scrna_ambient_correction_decontx_ordering_sensitivity` (hard_default / warn)
3. `scrna_ambient_correction_correct_then_qc_stricter` (hard_default / info)
4. `scrna_ambient_correction_high_ambient_biology_propagation` (flag_and_warn / warn)

**Synthesis:** `phase2/audit_ambient_correction/AUDIT_AMBIENT_CORRECTION_SYNTHESIS.md`
**Status detail:** `phase2/audit_ambient_correction/STATUS.md`. **Heumos:** extends.
**Deferred follow-ups:** parameter sensitivity, tool–counting interaction,
high-ambient replication — `phase2/DEFERRED.md`.

## Audit QC-MAD — at a glance

**Audit ID:** Audit QC-MAD (scRNA-seq low-quality cell filtering). **Status:** COMPLETE / CLOSED.
**Work window:** 2026-05-23 (extends Phase 1 p1).
**Checkpoints:** CP0 (inventory + Census re-pull) · CP1 (4-method comparison, 8 datasets) ·
CP2 (synthesis + §5.3.2 tiers) · CP3 (2 rules) · CP4 (synthesis + closeout) — **all complete.**
**Standards contributed:** none (used existing §5.3.2 amendment).

**Key finding:** the 4 QC filtering methods (pure-quantile, fixed-floor, MAD k=3/5)
produce largely equivalent cell sets (pair Jaccard 0.90–0.97); typical fixed-floor
defaults ≈ Heumos-recommended MAD k=3 (0.969). Disagreement is driven by dataset
gene-count distribution, not tissue/mito.

**2 rules:** `scrna_qc_filtering_method_equivalence` (conditional/info),
`scrna_qc_low_gene_dataset_caution` (conditional/warn). Both §5.3.2-tiered.
**Synthesis:** `phase2/audit_qc_mad/AUDIT_QC_MAD_SYNTHESIS.md`. **Heumos:** extends.

## Audit 3 — at a glance

**Audit ID:** Audit 3 (scRNA-seq counting tools). **Status:** COMPLETE / CLOSED.
**Work window:** 2026-05-16 → 2026-05-23.
**Checkpoints:** CP0 (env) · CP1 (inventory) · CP2 (refs/installs) · CP3
(acquisition + 4-tool counting) · CP4/C1 (count agreement) · CP5/C2
(cell-calling agreement) · CP6/C3 (biological propagation) · CP7 (4 rules) ·
CP8 (synthesis + closeout) — **all complete.**
**Standards amendment contributed:** §5.3.2 equivalence-finding tier criteria.

**Arc:** counts converge (C1, ρ~0.96) → native cell-callers nest STAR ⊂
alevin-fry ⊂ kb (C2; uniform caller → Jaccard 0.99) → on high-ambient tissue
the caller choice changes downstream biology (C3; washes out on clean tissue).

**4 rules** (all `hard_default` under §5.3.2 equivalence criteria, all `novel`):
1. `scrna_counting_tool_per_gene_count_convergence` (info)
2. `scrna_cell_calling_permissiveness_chain` (warn)
3. `scrna_uniform_cell_caller_eliminates_disagreement` (warn — the action)
4. `scrna_cell_calling_biological_propagation_high_ambient` (warn)

**Synthesis:** `phase2/audit3_counting/AUDIT3_SYNTHESIS.md`
**Status detail:** `phase2/audit3_counting/STATUS.md`
**Deferred follow-ups:** Audit 3c (CellRanger), 3d (multiome), 3e (5' v2),
3f (Tabula Muris re-pull) — `phase2/audit3_counting/DEFERRED.md`.

**Standards contributions:** §5.3.2 equivalence-finding tier criteria (adopted);
§3.5 "tight CIs ≠ correctness" + clean-control-for-stress-test (queued in
`standards/PENDING_AMENDMENTS.md`).

## Notes

- Rules are staged in audit-local `draft_rules/` and promoted to BioOrchestrator
  (`bioorchestrator/src/bioorchestrator/knowledge/rules/`) only via batched
  updates with a single version bump (not per-rule).
- BioOrchestrator is personal/internal; this index covers the audit corpus,
  which is kept public-ready (provenance, honest reporting, prior_audit tags)
  for a possible future public corpus.
