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

## Audit 3 — at a glance

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
