# E1 — Tool agreement on identical DEG inputs

**Database held constant:** MSigDB Hallmark (50 pathways, Hs + Mm)
**Tools:** fgsea (GSEA-family), gseapy.enrichr (ORA), clusterProfiler ORA
**Bootstrap:** dataset-level resampling, B=1000

## Pairwise agreement summary (medians with 95% bootstrap CI)

| Category | Tool pair | n | top_10 Jaccard | FDR<0.05 Jaccard | Spearman full | Direction agreement |
|---|---|---:|---:|---:|---:|---:|
| census_scrna | fgsea × clusterProfiler_ORA | 40 | 0.176 [0.176, 0.292] | 0.062 [0.000, 0.158] | 0.295 [0.223, 0.380] | 0.540 [0.500, 0.596] |
| census_scrna | fgsea × gseapy_enrichr | 40 | 0.250 [0.176, 0.333] | 0.056 [0.000, 0.167] | 0.316 [0.267, 0.435] | 0.591 [0.520, 0.670] |
| census_scrna | gseapy_enrichr × clusterProfiler_ORA | 40 | 0.818 [0.667, 0.909] | 0.500 [0.364, 0.600] | 0.952 [0.930, 0.972] | 0.860 [0.800, 0.903] |
| gtex_tissue_pair | fgsea × clusterProfiler_ORA | 30 | 0.333 [0.291, 0.429] | 0.000 [0.000, 0.057] | 0.333 [0.293, 0.368] | 0.480 [0.460, 0.500] |
| gtex_tissue_pair | fgsea × gseapy_enrichr | 30 | 0.333 [0.250, 0.429] | 0.049 [0.000, 0.077] | 0.354 [0.320, 0.398] | 0.480 [0.440, 0.490] |
| gtex_tissue_pair | gseapy_enrichr × clusterProfiler_ORA | 30 | 0.818 [0.603, 0.818] | 0.500 [0.000, 1.000] | 0.940 [0.834, 0.971] | 0.920 [0.890, 0.930] |
| tcga_cancer | fgsea × clusterProfiler_ORA | 48 | 0.333 [0.333, 0.429] | 0.127 [0.095, 0.155] | 0.424 [0.388, 0.464] | 0.560 [0.520, 0.600] |
| tcga_cancer | fgsea × gseapy_enrichr | 48 | 0.333 [0.333, 0.429] | 0.243 [0.185, 0.306] | 0.439 [0.418, 0.492] | 0.560 [0.520, 0.600] |
| tcga_cancer | gseapy_enrichr × clusterProfiler_ORA | 48 | 0.818 [0.818, 1.000] | 0.571 [0.422, 0.794] | 0.925 [0.916, 0.953] | 0.840 [0.800, 0.900] |

## Interpretation (per AUDIT_STANDARDS.md §3.4 — companion metrics framing)

Jaccard answers 'do the tools agree on which pathways pass FDR'.
Spearman answers 'do the tools rank pathways the same way'.
Direction agreement answers 'when both tools agree, do they assign
the same biological direction (up/down)'.

If Jaccard is low but Spearman + direction are high: tools largely
agree on the underlying biology and rank order, but disagree on
WHICH pathways pass the FDR cutoff — a multiple-testing-correction
and threshold artifact, not a fundamental biological disagreement.

If Jaccard AND Spearman are both low: tools fundamentally disagree;
user must pick one and document the choice.

## Caveats and gaps

1. **Direction agreement is fgsea × clusterProfiler_ORA only.** 
   gseapy.enrichr's log2_odds_ratio is currently uniformly positive 
   (Enrichr-style OR > 1 for any enriched term), so its 'direction' 
   carries no discrimination. Direction-agreement metrics with 
   gseapy as one tool collapse to a constant.
2. **top_50_jaccard omitted.** Hallmark has 50 pathways total, so 
   top-50 == full pathway set; Jaccard would be 1.0 by construction. 
   Replaced with top-10 / top-25.
3. **GSVA, EGSEA, camera not included.** Those require expression 
   matrices, not DEG TSVs; queued as E1b in DEFERRED.md.

