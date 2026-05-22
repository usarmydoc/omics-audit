# Audit 3 CP4 (C1) — Verification report

_Verification run: 2026-05-25 (read-only; no re-runs, no modifications)_
_Verifier scope: confirm CP4 outputs clean/complete/lock-registered; cold-read findings to inform CP5._

## Verification status: PASS with discrepancies surfaced

All 5 required outputs exist, are internally consistent, and are
hash-registered with zero drift. Five discrepancies vs the verification
prompt's expectations are surfaced below; **none invalidate the CP4
findings**, but three are genuine gaps worth a decision before CP5.

---

## Step-by-step results

### Step 1 — File existence: PASS
All present in `c1/`: `per_dataset_metrics.tsv`, `per_stratum_bootstrap.tsv`,
`per_gene_category_summary.tsv`, `tool_failure_modes.md`, `C1_findings.md`.

### Step 2 — Per-dataset metrics: PASS (with row-count discrepancy)
- **54 data rows, not 48.** ⚠️ The prompt expected 8 datasets × 6 pairs = 48.
  Actual is **9 datasets × 6 pairs = 54**. The 9th dataset is
  `gse288156_mouse_intestine_scrna`, added in CP3 as the Tabula Muris
  replacement under the user's "do not subset, pull all" directive.
  **This is correct, not an error** — the prompt's count was stale (it
  predates the intestine addition). The working set is:
  pbmc_1k_v3, pbmc_5k_v3.1, pbmc_10k_v3.1, pbmc_4k_v2, t_3k_v2,
  neuron_1k_v3, gse287209_lung_organoid, gse325955_mouse_kidney_E18_5,
  **gse288156_mouse_intestine_scrna**.
- All 6 tool pairs present, 9 datasets each. ✓
- All required metric columns present (per-gene Spearman, per-cell UMI
  Spearman, log2 median/IQR/p05/p95, det_rate A/B, direction
  agreement pct_A/B/tie). ✓
- Zero NaN in any metric column. ✓
- **Companion-metric note (§3.4):** direction-agreement (pct_A_higher /
  pct_B_higher / pct_tie) and magnitude correlation (rho_gene_total) are
  present. The literal "full-rank Spearman on the flattened gene×cell
  matrix" named in the CP4 scope was implemented as **Spearman on
  per-gene summed counts** (rho_gene_total), not on the flattened
  cell-level matrix. This is a reasonable and standard substitution
  but is a deviation from the literal scope wording — flagged for the
  record.

### Step 3 — Bootstrap CIs: PASS (with B-documentation gap)
- 48 rows = **12 strata × 4 metrics**. Strata = 2 chemistries (3p_v2,
  3p_v3) × 6 tool pairs = 12. ✓ Matches prompt's expectation.
- Metrics bootstrapped: rho_gene_total, rho_cell_total_umi, log2_median,
  pct_A_higher.
- CI columns `boot_ci_low` / `boot_ci_high` present; **0 violations** of
  low ≤ point ≤ high. ✓
- `insufficient_n` flag present; 24 rows flagged (all 3p_v2 strata, n=2).
- ⚠️ **B value is not a column in the TSV.** B=1000 is documented only
  in the `C1_findings.md` header ("B = 1000 (dataset-level)"). The
  bootstrap TSV itself does not carry B per-row. Minor provenance gap;
  the value is recoverable from findings but not self-contained in the
  data file.

### Step 4 — Gene category stratification: PARTIAL (one category missing)
- ⚠️ **4 categories present, not 5.** Present: `mitochondrial`,
  `overlapping`, `pseudogene`, `other`. **`multi-mappers` is absent.**
  The CP4 scope listed multi-mappers as a category, but it was never
  implemented: the GENCODE GTF does not directly annotate which genes
  receive multi-mapping reads, so there was no annotation source to
  build the category from. This was a known limitation noted in the
  analysis script but not carried forward into a documented caveat in
  `per_gene_category_summary.tsv` or `C1_findings.md`. **Genuine gap.**
- Fold rule (<100 genes → fold into other) applied correctly: **0
  violations**. `mitochondrial` (~37–111 genes) is flagged
  `fold_into_other=True` in 54 rows as expected.
- Per-category log2 distributions present for each tool pair. ✓

### Step 5 — Tool failure modes: THIN (gaps surfaced + one resolved here)
`tool_failure_modes.md` is sparse (11 lines). It documents only:
"all 9 datasets produced 4-of-4 tool outputs" and the barcode-universe
difference. It does **not** address:
- ⚠️ **kb_count smaller-output question (asked by the prompt).**
  Resolved during this verification: kb_count's smaller on-disk size
  is a **barcode-universe difference, not gene loss.** kb's
  `counts_unfiltered` holds only barcodes that received reads
  (~240K–1.4M) vs STARsolo's full whitelist raw output (6.8M barcodes).
  kb_count actually detects **more** genes than its partners
  (det_rate 0.42–0.69 vs partner 0.39–0.68 across the 9 datasets), so
  the smaller file is fewer barcodes + a different output format, not
  fewer genes. This resolution is documented here but **not yet in
  `tool_failure_modes.md`.**
- ⚠️ **The USA-mode bug** (the single biggest "tool produced
  systematically different output" event) is documented prominently in
  `C1_findings.md` and `superseded_.../SUPERSEDED.md`, but
  `tool_failure_modes.md` does not cross-reference it. A reader landing
  on failure_modes.md alone would miss the most important caveat.

### Step 6 — C1_findings.md cold read: STRONG (one §4.1 gap)
- Headlines clear and specific (corrected: "all four tools converge,
  rho ~0.96, no outlier"). ✓
- Stratified by chemistry (v2 vs v3) throughout. ✓
- Sample-size limitations surfaced per §3.1: v2 n=2 (CIs artificially
  tight), mouse n=3, tissue PBMC n=5 / non-PBMC n=4. ✓
  - Note: findings report **mouse n=3** (intestine restored it); the
    prompt anticipated mouse n=2. Positive deviation — more data than
    expected.
- Companion-metric framing per §3.4 present ("Where the tools agree"
  section). ✓
- ⚠️ **§4.1 prior-audit-relationship tag is MISSING** from the findings
  header. No prior counting-tool audit exists, so the correct tag is
  `novel`, but the explicit tag required by §4.1 (in findings.md section
  header + rule YAML metadata + phase summary) is not present.
- Scope honesty: the USA bug + supersession are documented prominently.
  However, the broader scope caveats (3 chemistries collapsed to 2 —
  5' v2 dropped; 3 binaries / 4 configs — CellRanger excluded) live in
  `CP4_HANDOFF.md` and `DEFERRED.md`, **not** restated in
  `C1_findings.md`. A reader of findings.md alone would not see the
  chemistry/tool-coverage caveats.

### Step 7 — Lock file verification: PASS (CLI-incompatibility noted)
- All 6 c1/ main outputs registered in `phase2/repro.lock` (5 required
  files + the `verification_pbmc_1k_v3_buggy_vs_fixed.tsv` table).
- **Zero hash drift** — manual SHA256 re-hash of every registered c1/
  file matches the lock. ✓
- ⚠️ **`repro verify` CLI fails** with "Missing required field:
  repro_schema_version / repro_version". The lock was written by the
  custom `dge_native.register_output()` helper, whose JSON format
  predates the `repro-lock` CLI's schema-version requirements. This is
  a **pre-existing, Phase-2-wide** condition (not CP4-specific) — the
  whole `phase2/repro.lock` uses this format. The substantive check
  (file hashes match) passes via manual re-hash; only the CLI wrapper
  is incompatible.

---

## Headline findings from C1 (cold read)

1. **All four counting tools converge on per-gene counts.** Cross-tool
   per-gene Spearman ≈ 0.96 (v3, n=7), bootstrap CIs [0.95, 0.97], all
   overlapping. The same-tool baseline (STAR default vs CR-mimic) is
   0.991. No tool is an outlier.
2. **Tools also agree on per-cell UMI ranking** (rho 0.93–0.97 across
   all 6 pairs). Cell depth is preserved regardless of per-gene
   attribution.
3. **No gene-category-specific disagreement** after the USA-mode fix —
   all categories show median log2 = 0.0, IQR ≤ 0.07.
4. **The original "alevin-fry is the outlier" headline was a CP4
   analysis bug** (USA-suffix collapse restricting alevin-fry to its
   ~10% Ambiguous bucket), found and fixed; buggy outputs preserved
   under `superseded_2026-05-18_buggy_usa_strip/`.
5. **Chemistry effect is small** (v3 mean 0.957, v2 mean 0.92, within
   the v2 n=2 noise floor).

## Tool disagreement magnitude: NEAR-EQUIVALENCE

Cross-tool rho ≈ 0.96 with overlapping CIs, near-zero log2 offsets, and
category-level agreement at noise floor. By the prompt's own decision
tree, this is the "near-equivalence" branch.

## Surprises that change CP5/CP6 expectations

- **The biggest surprise was procedural, not biological:** the headline
  flipped from "alevin-fry outlier" to "all tools equivalent" only after
  catching a preprocessing bug in the audit's *own* analysis code. CIs
  were tight on the buggy data — bootstrap machinery does not protect
  against systematic preprocessing error. This is a manuscript-grade
  methodological lesson (already captured in `CP4_HANDOFF.md`).
- **Per-cell UMI agreement is already high (0.93–0.97).** This pre-empts
  part of C2: tools already agree on cell *depth*. C2's real question
  narrows to the *barcode-calling* step at the low-UMI / ambient
  boundary.

## Recommendation for CP5 scope

Per the prompt's decision tree, near-equivalence on counts points to:
**C2 (CP5) is not a from-scratch disagreement hunt; it is a focused test
of whether cell-*calling* defaults diverge where counts do not.**

Concrete starting points (NOT commitments — CP5 prompt drafted next,
per standing directive):

1. **Drop the "3-trusted-vs-1-outlier" framing entirely** — it came from
   the buggy pass. Report all 6 pairwise barcode agreements symmetrically.
2. **Focus the high-information stratum on the low-UMI / ambient regime.**
   Tools agree on deep cells; divergence (if any) is at the knee.
3. **Lead with companion metrics** (per-tool cell count, knee-point,
   ambient estimate, barcode rank correlation), with Jaccard stratified
   by UMI regime rather than as a single headline number.
4. **CP6 (Phase 1 robustness) gains relative importance.** With counts
   near-equivalent, the higher-value question becomes "do any
   cell-calling differences CP5 finds propagate into clustering / marker
   genes / annotation, or wash out downstream?" — which may make CP6 the
   sharper audit of the two. Decide CP6 framing at CP5 closeout.

## Recommended fixes before/at CP5 (for user decision — not done here)

These are gaps surfaced above, left unmodified per the read-only rule:

1. Add the `multi-mappers` category caveat to `per_gene_category_summary`
   docs / findings (or document why it was dropped). — REAL GAP
2. Add `prior_audit_relationship: novel` tag to `C1_findings.md` header
   per §4.1. — REAL GAP
3. Expand `tool_failure_modes.md` to (a) record the kb_count
   barcode-universe resolution above and (b) cross-reference the USA-mode
   bug. — REAL GAP
4. Optionally add a `B` column to `per_stratum_bootstrap.tsv` for
   self-contained provenance. — MINOR
5. Optionally restate the chemistry/tool-coverage scope caveats (5' v2
   dropped, CellRanger excluded) in `C1_findings.md`. — MINOR
6. The `repro verify` CLI incompatibility is Phase-2-wide; if CLI
   verification is desired, the `dge_native` lock writer needs
   `repro_schema_version` fields. Out of Audit 3 scope. — INFRA

---

## Bottom line

CP4 outputs are **clean, complete, and lock-registered with zero drift.**
The corrected finding is **near-equivalence across all four tools** on
per-gene counting. Three documentation gaps (multi-mappers caveat, §4.1
tag, thin failure-modes file) are real but do not affect the findings.
CP5 should be scoped as a focused barcode-calling test at the low-UMI
boundary, not a broad disagreement hunt. **Stopping here for CP5 scope
decision.**
