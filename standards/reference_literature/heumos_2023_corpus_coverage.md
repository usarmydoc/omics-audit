# Heumos 2023 Corpus Coverage Summary

_Working artifact for future audit selection. Resolved 2026-05-23 against
Phase 1 (scrnaseq_audit p1–p5), Phase 2a (A2/A5/A6/A7), Audit 1 (E2),
Audit 3 (C1–C3). Pull from the "candidates" section when scoping the next
audit — not from the full DEFERRED list._

## What's covered

| Heumos area | Audit(s) | Relationship |
|---|---|---|
| Doublet detection | Phase 1 p2 (scDblFinder vs Scrublet, 19 datasets) + Phase 2a A5 (Demuxlet ground truth) | confirms (+extends with ground truth) |
| Clustering resolution selection | Phase 1 p9 + Phase 2a A2 (ARI vs V-measure 5× divergence) | extends |
| DGE (pseudobulk vs cell-level) | Phase 1 p4 + Phase 2a A6 + pseudobulk_vs_wilcoxon (26 datasets) | extends |
| Pathway enrichment (database choice) | Audit 1 E2 (24–124× across MSigDB) | extends |
| Read counting agreement | Audit 3 C1 (ρ~0.96, 9 datasets) | fills_gap |
| Cell calling | Audit 3 C2/C3 (permissiveness chain + biological propagation) | extends |
| Batch integration | Phase 1 p3 (bbknn/harmony/scanorama/scvi via ARI/LISI/kBET) | confirms |
| QC filtering method (MAD vs quantile) | Audit QC-MAD (8 Census datasets; pair Jaccard 0.90–0.97; C2≈MAD3 0.969) | extends |

## What's partially covered

Corpus touches the area; Heumos's recommendation has un-addressed dimensions.
All captured in `phase2/DEFERRED.md` → "Partially-covered gaps."

| Area | Covered | Heumos dimension missed | DEFERRED? |
|---|---|---|---|
| QC — ambient RNA | mito threshold (p1) + filtering-method equivalence (Audit QC-MAD) | ambient RNA correction (SoupX/CellBender/DecontX) still uncovered | yes (ambient = high-value; MAD-vs-quantile now RESOLVED by Audit QC-MAD) |
| Doublet | scDblFinder>Scrublet (p2/A5) | multi-method ensemble | yes (small) |
| Clustering | resolution + metric (p9/A2) | Leiden-vs-Louvain algorithm itself | yes (small) |
| Annotation | MarkerScore vs SingleR accuracy (p5) | CellTypist; 3-step workflow | yes (medium) |
| Batch | scRNA integration methods (p3) | cell-cycle regression | implied |
| Pathway | database + paradigm (Audit 1 E2) | decoupleR multi-method ensemble | yes (small-medium) |

## What's not covered

No corpus coverage; captured in `phase2/DEFERRED.md`.

- Normalization-by-purpose (shifted-log / Pearson-residual / scran) — durable.
- Feature selection (deviance HVG benchmark) — durable; used as CP6 method only.
- Compositional / differential-abundance (scCODA/MILO/DA-seq) — distinct step.
- Trajectory inference; RNA velocity; cell–cell communication; perturbation
  modeling — fast-moving subdomains, partly churned since 2023.
- scATAC-seq; CITE-seq; AIRR; spatial — modalities outside current corpus.

## Audit selection criteria for future work

- **Prioritize** gaps where Heumos's recommendation is qualitative and the
  corpus can quantify it (the corpus's repeated winning move: B1/E2/A6/Audit3).
- **Prioritize** gaps that tie into an existing finding (compounding value).
- **Prioritize** durable/architectural gaps over fast-moving ones (the
  finding survives longer).
- **Deprioritize** fast-moving subdomains where 2023 tool lists have churned
  (velocity, CCC, perturbation, spatial) — they need a current-methods scan
  first and the finding ages fast.
- **Deprioritize** modalities without hands-on operator context (scATAC, AIRR).

## Specific candidates worth considering when bandwidth allows

1. **Ambient RNA correction** (SoupX vs CellBender vs DecontX). *Highest
   value:* directly extends Audit 3 C3's ambient-burden finding (does
   correcting ambient before cell-calling change the C3 propagation?),
   qualitative-in-Heumos → quantifiable, durable, small-medium. Compounds an
   existing headline.
2. **Normalization-by-purpose** (shifted-log vs Pearson-residual vs scran →
   downstream clustering/DE). Architectural, durable, the kind of "does the
   default choice matter" question the corpus answers well. Medium.
3. **Annotation tool + workflow** (CellTypist vs SingleR vs MarkerScore;
   does the 3-step workflow change calls). Audit 3 CP6 already leaned on
   CellTypist un-benchmarked — closes that loop. Medium.
4. **Compositional analysis** (scCODA vs MILO vs DA-seq). A whole analysis
   step the corpus skips; bounded; relevant wherever condition contrasts
   matter. Medium.
5. **decoupleR enrichment ensemble** (does consensus beat single method).
   Small-medium extension of Audit 1; tests Heumos's specific decoupleR
   recommendation the E2 audit didn't.

Lower priority: Leiden-vs-Louvain (small but low stakes — Leiden is already
near-universal); doublet ensemble (small, scDblFinder already wins);
trajectory/velocity/CCC/perturbation (fast-moving, need current-methods scan);
modality-specific (no operator context).
