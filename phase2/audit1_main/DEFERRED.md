# Audit 1 main — DEFERRED scope-creep candidates

Per the spec: scope-creep candidates land here without action during this
work. Reviewed only at user discretion between checkpoints.

---

## Queued future audits (real audits, scoped, sized — NOT acted on here)

### E1b — Expression-matrix-based pathway methods on a representative subset

**Queued during:** CP3 scope reconciliation, 2026-05-16

**Why this is its own audit, not part of E1:**
The original E1 spec listed 6 tools (fgsea, GSVA, EGSEA, enrichr,
clusterProfiler ORA, camera). Three of those — GSVA, EGSEA, camera —
fundamentally cannot consume DEG TSVs as input; they require the
underlying expression matrix + design. E1 as originally scoped tested
"tool agreement on identical DEG inputs," so the three matrix-input
tools never fit cleanly. Rather than re-fetching expression matrices
for all 125 input comparisons (a 1-2 week task that would push Audit 1
timeline by 50%), the matrix-input tools are split into E1b as a
separate audit with its own scope, inputs, and methodology.

**Pattern reference:** Same approach Phase 2a took with A1-A7
amendments — closeout work that didn't fit cleanly into the original
audit got captured as separate scoped items, executed with their own
methodology, and produced clean output. E1b is Audit 1's equivalent.

**Proposed E1b scope (locked at audit time, not now):**

- **Tools:** GSVA, EGSEA (ensemble: camera, roast, fry, safe, padog,
  gage, plage, zscore, ssgsea, globaltest), camera (limma)
- **Inputs:** representative subset, e.g., 5 TCGA cancers + 5 Census
  datasets + 5 GTEx pairs = 15 expression matrices. Pre-fetched
  count matrices + design metadata via cellxgene_census / recount3.
- **Database:** MSigDB Hallmark (50 pathways), held constant across
  tools — direct analog to E1's tool-agreement methodology.
- **Metrics:** same as E1 — pairwise Jaccard, FDR overlap, NES
  direction (where applicable), Spearman/Pearson on full rankings.
- **Bootstrap:** dataset-level, B≥1000, per AUDIT_STANDARDS.md §2.4.
- **Stratify by:** input category (TCGA / Census / GTEx) AND by
  whether the tool is GSEA-style (camera) vs single-sample (GSVA,
  ssgsea via EGSEA) vs ensemble (EGSEA).
- **Estimated effort:** 1.5-2.5 weeks (data fetch + 15 datasets × 3
  tools × pairwise analysis + bootstrap + findings).
- **Output:** findings.md (E1b section), draft rule YAMLs if
  warranted (rule corpus location TBD — likely a new
  `audit1_main/e1b/` directory or a separate `audit1_supplement/`).

**Sequencing:** E1b runs only after Audit 1 main is fully complete
(all 7 checkpoints, all rules drafted, lock state stable). Not a
Phase 2b schema candidate — it's a real audit deliverable with its
own scope.

---
