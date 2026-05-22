# CP4 / C1 — Reading order for results

Drafted before CP4 finishes so the post-completion read is fast and structured.
Goal: extract the four decisions needed to scope CP5 (C2 cell-barcode calling).

---

## Read in this order

### 1. `tool_failure_modes.md` — first

**Why:** if any (dataset × tool) failed or produced degenerate output, the
per-dataset numbers downstream are biased. Confirm what's missing before
interpreting headlines.

**Look for:**
- Datasets where fewer than 4 tools produced matrices
- kb_count detection-rate red flags (already pre-scanned for in the script)
- Any per-tool weirdness (alevin-fry's ~190K-gene salmon index, etc.)

### 2. `per_dataset_metrics.tsv` — second

**Sort by `rho_gene_total` ascending.** The worst (lowest correlation) pairs
across all 54 rows tell you where disagreement concentrates.

**Decisions to make:**
- Is there one tool that consistently disagrees with the other three?
  (Look at row-wise tool_a and tool_b frequencies in the bottom decile)
- Or is disagreement diffused across pairs?

**Companion columns to scan alongside rho_gene_total:**
- `rho_cell_total_umi` — do tools agree on per-cell depth even when per-gene rho is low?
- `pct_A_higher` / `pct_B_higher` — directional bias check (close to 0.5 → balanced; skewed → one tool systematically higher)
- `log2_median` — if non-zero, that's a systematic offset (e.g., kb_count consistently 0.3 log2-units below STARsolo)
- `det_rate_A` vs `det_rate_B` — do tools detect the same fraction of genes?

### 3. `per_gene_category_summary.tsv` — third

**The audit-decision payload.** This tells us whether disagreement is
gene-category-dependent.

**Decisions to make:**
- **If pseudogenes/overlapping show much larger log2_iqr than 'other'**:
  CP5 should stratify cell-calling agreement by gene-category fingerprint
  (high-ambient-RNA cells will pick up different pseudogenes per tool).
- **If 'other' is just as bad**: tools disagree on regular genes too;
  CP5 doesn't need to stratify by category; instead investigate per-cell-rank
  agreement at low-count cells.
- **If mitochondrial differs sharply**: a flag for downstream QC (mito% varies
  by tool, affecting cell-filtering decisions). Surface in findings.

**Compare medians across pairs:** if log2_median for (star_default vs kb_count)
is consistently > 0 across categories and datasets, STARsolo systematically
reports higher per-gene counts than kb_count.

### 4. `per_stratum_bootstrap.tsv` — fourth

**Whether stratified findings are robust to dataset resampling.**

**Decisions to make:**
- Is the chemistry split (v2 vs v3) meaningful, or are CIs overlapping?
- For v2 with n=2, the CI will be very wide (per AUDIT_STANDARDS §3.1 this is
  a headline finding, not a caveat). Surface it.
- Do `rho_cell_total_umi` CIs straddle 1.0 or are they tight?

### 5. `C1_findings.md` — fifth (last, then expand)

The script writes a stub findings.md with tables but no synthesis. After
reading 1-4, write the synthesis section in `C1_findings.md` covering:

- **Headlines (3-5 bullet points):** the punchiest summaries that survive
  bootstrap CIs. Don't soft-pedal weak findings — surface them clearly.
- **Where tools agree:** companion-metric framing per §3.4 — direction
  agreement, magnitude correlation, full-rank Spearman, gene-detection-rate
  overlap. Make this section as substantial as the disagreement section.
- **Where tools disagree:** specific gene categories, specific tools, specific
  chemistries.
- **Sample size limitations:** mouse n=3, v2 n=2 — these are the constraints,
  not afterthoughts.
- **Tool failure modes:** synthesize from §1's narrative.
- **Implications for CP5 (C2):** one short paragraph. What does the C1
  finding mean for how to scope cell-barcode-calling agreement?

---

## After writing findings: the four decisions that gate CP5 scope

These determine the CP5 prompt:

1. **Is one tool the outlier, or is disagreement diffuse?**
   - One outlier → CP5 can stratify "agreement among 3 trusted tools" vs
     "agreement that includes the outlier"
   - Diffuse → CP5 reports all 6 pair-wise Jaccards equally

2. **Is disagreement concentrated in specific gene categories?**
   - Yes → CP5's barcode-agreement analysis can use category-stratified
     "ambient signature" as a stratification axis
   - No → CP5 stays as scoped: pure Jaccard on called barcodes per pair

3. **Are tools "essentially equivalent" or "substantively different"?**
   - Essentially equivalent (rho > 0.95 typical) → CP6 (Phase 1 robustness)
     becomes a confirmation pass; CP5 may also be lighter
   - Substantively different → CP6 becomes a deeper investigation; CP5 needs
     more attention to ambient/low-UMI barcode regimes

4. **Is chemistry the dominant axis of variation, or species/tissue?**
   - Chemistry → CP5 stratifies by chemistry
   - Species/tissue → CP5 stratifies by species/tissue (and acknowledges n)
   - Neither — disagreement is dataset-specific → CP5 reports per-dataset
     without grouping

---

## Anti-checklist (things NOT to do during CP4 readout)

- Don't infer mechanism from per-gene correlation alone. Mechanism candidates
  go in C1_findings.md hedged per §3.5; the rule wording (CP8) stays on
  observations.
- Don't draft CP5+ prompts until findings.md is written. The CP5 scope flows
  from the four decisions above, not from CP4 numbers in isolation.
- Don't draft BioOrchestrator rules yet. Wait for CP8 with all of C1+C2+C3
  evidence assembled, per [[feedback-bioorchestrator-batch-updates]].
- Don't fix any CP3 output even if it looks weird. Anomalies go in
  tool_failure_modes.md as documentation, not corrections.

---

## Quick-reference: where each output lives

- `c1/per_dataset_metrics.tsv` — 9 × 6 = 54 rows
- `c1/per_gene_category_summary.tsv` — 9 × 6 × 4 categories = up to 216 rows
- `c1/per_stratum_bootstrap.tsv` — chemistry × tool_pair × 4 metrics
- `c1/tool_failure_modes.md` — narrative
- `c1/C1_findings.md` — synthesis (stub → expanded)
- `c1/gene_categories_cache.pkl` — GTF parse cache (don't read; regenerated on demand)

All TSV/MD outputs are hash-registered in `/mnt/nvme1/omics-audit/phase2/repro.lock`.
