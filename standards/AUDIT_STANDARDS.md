# AUDIT_STANDARDS.md

Standards for all bioinformatics audit work in this workspace. Every audit — current and future — meets these standards. Reference this document from each audit's local CLAUDE.md so it loads at session start.

These standards exist because Phase 1 produced unhashed outputs, undocumented sample-size gaps, and over-framed findings. Phase 2 and 2a produced provenance-clean, sample-size-honest, and refinement-disciplined work. The difference is the standards below. Maintain them.

Depth over length. Four sections (rule YAML schema, bootstrap unit specification, confidence tier boundaries, closeout amendment triggers) are fully specified rather than gestured at. The rest stays terse on purpose.

---

## 1. Provenance and reproducibility

### 1.1 Hash registration at write time

Every output file (TSV, parquet, HDF5, JSON, or other artifact) is hash-registered in the audit's lock file the moment it is written, not in a separate snapshot step.

- Use `repro.register_output()` or equivalent inline in the write function
- Never rely on a separate `repro snapshot` step that runs hours later
- The lock file's `verified_outputs` section must be populated when the audit completes; if it is empty, the audit is not complete

### 1.2 Environment stamping

Every output carries an environment ID stamp that resolves to a full session info dump.

- R + Bioconductor + tool versions
- Python + key package versions
- OS, architecture, BLAS library
- One environment ID per consistent toolchain; audits using multiple toolchains carry multiple IDs

The lock file contains an `environments` registry mapping env_id to full session info, including R `sessionInfo()` and `pip freeze` (or `conda list`) dumps.

### 1.3 Native runtimes only

No cross-language bridges (rpy2, reticulate, in-process language interop).

- Python tools run in Python
- R tools run via `Rscript` subprocess
- Data exchange via files (TSV, parquet, HDF5)
- Rationale: cross-language bridges have produced documented bugs (Phase 1 rpy2 column-major bug, harmonypy 0.2.0 tensor handling) that contaminate findings

### 1.4 Drive and path discipline

- All audit outputs under a single canonical path (e.g., `~/omics-audit/`)
- Symlinks preserve old paths when consolidating
- No hardcoded absolute paths that break when drives are reorganized
- Lock file paths are relative to the audit root, not absolute

### 1.5 Lock file completeness checklist

Before declaring an audit complete:

- [ ] Every output file hash-registered
- [ ] `verified_outputs` populated
- [ ] `environments` registry populated
- [ ] `repro verify` exits zero
- [ ] Session info dumps archived
- [ ] If any output is regenerated, old version retained as `superseded_<date>` with retention rationale

---

## 2. Audit design

### 2.1 Pre-specified sample size

Every audit declares its sample size target before execution.

- "≥10 datasets across ≥3 tissues" or similar
- Not "as many as we can grab"
- Document why the target is sufficient for the claims being made
- If target is not met, surface the limitation in findings.md and do not encode rules that require the unmet evidence

### 2.2 Per-dataset transparency

Aggregate metrics are reported alongside per-dataset metrics, never instead of them.

- Readers should be able to see whether findings hold across all datasets or are driven by outliers
- Per-dataset TSVs are first-class outputs, not supplementary
- Aggregate statistics report bootstrap CIs, not just means

### 2.3 Multiple complementary metrics

Single-metric claims are insufficient. Every audit reports multiple complementary metrics where the question admits them.

For tool comparison questions specifically:
- Top-K Jaccard at K=50, K=100, K=200 (or appropriate range)
- FDR-thresholded overlap (padj < 0.05, padj < 0.01)
- Direction of effect agreement
- Spearman or Pearson correlation on effect size estimates (log2FC, etc.)
- Full rank correlation, not just top-K

For clustering quality questions:
- ARI alongside V-measure, homogeneity, completeness, silhouette
- Document which metric is monotonic in the parameter being swept (and therefore useless as a standalone selector)

For classification questions:
- AUROC, F1, precision, recall — not just one
- Report in-sample vs cross-validated separately
- Compare against most-frequent-class baseline

### 2.4 Bootstrap on the unit of analysis — full specification

Confidence intervals are computed by bootstrapping the unit of analysis, defined precisely below. Resampling at the wrong unit produces artificially tight CIs.

**Unit selection rules:**

| Audit type | Unit of analysis | Resample procedure | Minimum B |
|---|---|---|---|
| Tool-vs-tool across N datasets (e.g., B1, A1) | Dataset | Resample N datasets with replacement; compute per-resample metric; aggregate | 1000 |
| Pseudobulk DE within a single dataset | Donor | Resample donors with replacement; recompute pseudobulk DE per resample | 1000 |
| Cell-level DE with mixed model (e.g., NEBULA) | Donor (random effect group) | Resample donors; cells within sampled donors are included as a block | 500 |
| Cell type marker scoring across datasets (e.g., P5) | Dataset | Resample datasets; do not resample within-dataset cells | 1000 |
| Threshold sensitivity across datasets (e.g., P1, A3) | Dataset | Resample datasets; do not resample cells within | 1000 |
| Per-tissue audit (n datasets per tissue) | Dataset within tissue OR tissue, depending on claim scope | If claim is per-tissue: dataset within tissue (with adequate n per tissue). If claim is across tissues: tissue. | 1000 |
| Doublet ground truth benchmark (e.g., A5) | Cell (only when comparing classifications against known labels, never for tool-disagreement CIs) | Bootstrap cells for AUROC CI; bootstrap datasets for cross-dataset claims | 1000 |
| Clustering quality across datasets | Dataset | Resample datasets; do not resample cells within | 1000 |
| Integration method comparison across datasets | Dataset | Resample datasets; do not resample within | 1000 |

**Hierarchical structure rules:**

When the audit has nested structure (cells in donors in datasets in tissues), bootstrap at the level that matches the claim's scope:

- Claim about "this method on this dataset": bootstrap donors or cells (whichever is the unit the model treats as independent)
- Claim about "this method across datasets": bootstrap datasets (not cells within)
- Claim about "this method across tissues": bootstrap tissues (when n per tissue ≥ 3) OR datasets across tissues (when n per tissue < 3, with the caveat that tissue-level inference is limited)

Never bootstrap at multiple levels simultaneously unless the bootstrap design is explicitly hierarchical (e.g., two-stage bootstrap for variance components). Most audits should use single-level bootstrap.

**Forbidden bootstrap procedures:**

- Resampling cells without preserving donor structure for any DE-related claim (introduces pseudoreplication into the CI itself)
- Resampling genes for tool-agreement claims (genes are not the unit; tools are compared on the gene-rank ordering)
- Subsampling instead of bootstrapping (subsampling produces variance estimates, not confidence intervals)
- B < 500 for any published CI; B < 100 for any preliminary CI

**Required reporting:**

Every CI in findings.md includes:
- Bootstrap unit explicitly named
- B value
- CI method (percentile, BCa, or basic)
- Whether B was sufficient for the precision claimed (rule of thumb: B ≥ 1000 for 95% CIs reported to two decimal places)

### 2.5 Diagnostic runs in parallel with main audit

The diagnostic that produced "tools agree on biology, disagree on padj-ranking" in Phase 2 emerged from running tool-agreement diagnostics alongside the main audit. Pattern this:

- Compute direction agreement, log2FC correlation, and full rank correlation alongside any Jaccard-style disagreement metric
- Diagnostics often produce sharper framings than the headline metric
- Diagnostic outputs are first-class, hash-registered, environment-stamped

### 2.6 Failure mode documentation

Tools fail. Datasets fail. Document them.

- Which datasets failed each tool, and why
- Which tools produced degenerate output and at what point (e.g., Scrublet's auto-threshold maxing below the actual doublet score on PBMC data)
- Whether failures are tool bugs, version issues, or dataset properties
- Failure rates per tool per modality

---

## 3. Honest reporting

### 3.1 Sample size limitations surfaced loudly

If a per-tissue analysis has n=1 or n=2 per tissue, that is the headline finding, not a caveat. The Phase 2a A3 pattern is the model: "the rule we wanted to encode cannot be encoded from this data; here is what would be needed."

- Do not encode rules that require evidence the audit does not have
- Quantile-on-data approaches (compute the threshold from the working dataset) often beat per-category lookups when sample sizes are inadequate
- Insufficient-data findings are valid findings, not failures

### 3.2 In-sample vs cross-validated metrics

When a model or classifier is trained as part of the audit, report:

- In-sample accuracy / F1 / AUROC
- Leave-one-out or k-fold cross-validated accuracy / F1 / AUROC
- Most-frequent-class baseline accuracy
- Whether the model beats baseline under cross-validation

If the model does not beat baseline under cross-validation, do not deploy it. Report features as correlates only, not as predictors.

### 3.3 K-sensitivity, threshold-sensitivity, preprocessing-sensitivity

Every finding that depends on a parameter cut should be tested for sensitivity to that cut.

- Top-K Jaccard: test at multiple K values
- FDR-thresholded findings: test at multiple FDR levels
- Preprocessing-dependent findings: test under method-native vs default preprocessing
- Document where findings are stable across cuts and where they are not

### 3.4 Direction agreement reported alongside disagreement

When tools or methods disagree, report what they agree on too.

- Direction of effect agreement
- Magnitude correlation (log2FC Pearson, Spearman)
- Where the disagreement is concentrated (padj ranking at the margin vs effect direction vs effect magnitude)
- This produces sharper, more useful framings than single-axis disagreement claims

### 3.5 Mechanism candidates appropriately hedged

When a finding has a likely mechanism, document it with appropriate hedging.

- "Most likely candidate per documentation review, not isolated experimentally" — acceptable
- "Caused by X" — only with experimental isolation
- Mechanism speculation stays in findings.md, not in rule wording
- Rule wording stays on the observation

---

## 4. Follow-up discipline

### 4.1 Prior audit relationship tagging

Every finding that touches a prior audit claim carries a `prior_audit_relationship` tag with one of:

- **as_original**: rule matches prior audit finding exactly
- **refines_prior**: original claim is correct but needs nuance
- **contradicts_prior**: specific claim does not hold under broader testing
- **extends_prior**: new finding not in prior audit
- **novel**: rule does not touch any prior audit finding

The tag appears in:
- Rule YAML metadata
- findings.md section header
- The phase summary document

### 4.2 Contradictions surfaced before integration

If a follow-up audit contradicts a prior claim, surface it immediately. Do not proceed to rule encoding until:

- The contradiction is documented in findings.md
- The original claim, the contradicting evidence, and the strength of the contradiction are spelled out
- The supervised reviewer has decided how to handle it

### 4.3 Closeout amendment triggers — full specification

Before any audit integrates into BioOrchestrator, run a closeout amendments pass. The triggers below identify findings that require amendment. Auto mode can run this list mechanically; human review verifies the calls.

**Mandatory amendment triggers (any one is sufficient):**

| Trigger | Action required |
|---|---|
| Finding refines a prior audit claim | Promote to top-level section in findings.md; tag `refines_prior`; draft dedicated rule YAML; note in phase summary |
| Finding contradicts a prior audit claim | Promote to top-level section; tag `contradicts_prior`; surface immediately; do not encode any related rule until reviewer decision |
| Finding currently in a sub-bullet but tested on ≥10 datasets across ≥3 tissues | Promote to top-level section; draft dedicated rule YAML if not already drafted |
| Headline metric and a complementary metric tell substantially different stories | Top-level section explaining the divergence; rule reports both metrics, not just one |
| Sample size limitation makes rule encoding inappropriate | Tag `insufficient_data`; document what would be needed; recommend quantile-on-data or comparable fallback |
| Subset finding being extended to full set produced different results | Stop; document contradiction; do not encode the subset finding as a rule |
| Tool failure rate > 20% on the audit's dataset corpus | Document as a finding in its own right; rule notes failure modes per tool |
| Preprocessing-sensitivity test shows a published claim is preprocessing-dependent | Top-level section; tag `refines_prior` or `contradicts_prior` based on magnitude; rule encodes the dependence |
| Cross-validated metric falls below most-frequent-class baseline | Do not deploy as predictor; tag features as correlates only; rule reflects this |

**Discretionary amendment triggers (judgment-based):**

- Finding is interesting but doesn't meet the mandatory triggers above
- Methodological observation worth surfacing for future audits (e.g., the "Bioconductor archive policy" observation from Phase 2)
- Rule wording that came out of auto mode that needs human refinement

**Output of the closeout pass:**

- A dedicated `closeout_amendments.md` document listing every triggered amendment, the trigger that fired, and the action taken
- Updates to findings.md and STATUS.md integrating the amendments
- New or revised rule YAMLs as triggered
- Cross-reference back to the original finding location in the main summary

The closeout pass is a distinct, documented step — not folded silently into the main audit work.

### 4.4 Subset findings re-tested at full scale

When a finding is established on a subset (because the full set was computationally expensive at the time), re-test on the full set when the compute becomes available.

- If the subset finding holds: confirm and proceed
- If it does not hold: stop, document the contradiction, do not encode the subset finding as a rule
- The Phase 2a A6 NEBULA/muscat extension from 5 to 21 datasets is the pattern

---

## 5. BioOrchestrator integration

### 5.1 Supervised review required

Draft rule YAMLs are never auto-merged into BioOrchestrator. Every rule requires a supervised review pass before integration.

The review checks:

- Trigger conditions match the evidence
- Evidence pointers resolve correctly
- Recommendation is actionable (specifies how, not just what)
- Out-of-scope acknowledgment is present
- Mechanism speculation has not crept into rule wording
- `prior_audit_relationship` tag is correct

### 5.2 Evidence pointers

Every rule cites specific output paths and lock file entries as evidence.

- File paths relative to the audit root
- Lock file entry by hash
- Section reference in findings.md
- An auditor reading the rule should be able to trace it to the empirical basis in one click

### 5.3 Confidence tiers — full specification

Rules are tiered by evidence strength. Boundaries below are numerical, not soft, so that two audits use each tier to mean the same thing.

The tier values are unchanged (`hard_default` / `conditional` / `flag_and_warn` / `literature_based` / `insufficient_data`). What differs by audit type is *which numerical boundaries apply*:

- **§5.3.1 (tool-selection criteria)** applies to audits whose finding is "tool/method/parameter X wins" — the default, used by most audits.
- **§5.3.2 (equivalence-finding criteria)** applies to audits whose finding is equivalence, agreement, or convergence between tools/parameters — where there is no winning tool, so §5.3.1's dominance criteria are structurally inapplicable.

When an audit produces both selection and equivalence findings, each finding is tiered by the applicable subsection.

#### 5.3.1 Tool-selection tier criteria

**Tier: hard_default**

All of the following must hold:

- Tested on ≥15 datasets
- Datasets span ≥3 tissues
- Finding is stable across cuts: top-K Jaccard agreement remains within ±0.1 across K ∈ {50, 100, 200}; if FDR-thresholded, stable across {0.01, 0.05, 0.1}
- A single tool or method dominates: wins in >60% of datasets OR has mean metric advantage > 0.1 (Jaccard) or > 0.05 (ARI) over the next-best with non-overlapping bootstrap 95% CIs
- Preprocessing-sensitivity tested where applicable; finding holds under method-native and default preprocessing
- Bootstrap CIs reported with B ≥ 1000 on dataset-level resampling
- No `contradicts_prior` from a follow-up audit

Rule wording: declarative, recommends the specific tool/method/parameter.

**Tier: conditional**

All of the following must hold:

- Tested on ≥10 datasets
- Performance depends on an identifiable dataset feature (e.g., n_donors, n_cells_per_donor, tissue type)
- The dependence is statistically supported: significant correlation at p < 0.05 with effect size r > 0.4 (Pearson) or ρ > 0.5 (Spearman), OR significant difference between feature-defined subgroups (Kruskal-Wallis p < 0.05)
- The dependence has interpretable mechanism (not just empirical correlation)
- In-sample vs cross-validated performance reported separately if a model is used; CV must beat most-frequent-class baseline

Rule wording: conditional, with the feature-based decision logic embedded.

**Tier: flag_and_warn**

Triggered by any of:

- Substantial tool disagreement: top-K Jaccard < 0.5 at K=100 AND no single tool wins > 60% of datasets
- Direction agreement < 0.95 across tools (note: this is rare; most tools agree on direction even when disagreeing on ranking)
- Bootstrap CIs across tools overlap substantially
- Mechanism unclear and audit cannot disambiguate
- Findings refine or contradict prior audit claims with strength that warrants user awareness but no specific corrective action

Rule wording: surfaces the uncertainty, does not pick sides, recommends user awareness and reporting of choice.

**Tier: literature_based**

Used when:

- This audit did not test the convention directly
- A peer-reviewed source provides the basis
- Rule cites the source with full bibliographic information

Rule wording: cites the literature, makes clear this is not an audit-derived rule.

**Tier: insufficient_data**

Triggered when:

- Audit attempted to test the convention but sample size per relevant category is < required threshold (typically n=3 per category for any per-category claim; n=10 overall for any aggregate claim)
- Audit identified the convention as worth testing but data was not available
- A subset analysis raised the question but full analysis is not yet possible

Rule wording: explicitly acknowledges insufficient evidence; recommends quantile-on-data approach OR comparable fallback that does not require the absent evidence; does not encode a specific value that the data does not support.

**Tier escalation and demotion:**

- A `flag_and_warn` rule can be promoted to `conditional` when additional audit work establishes the conditional logic
- A `conditional` rule can be promoted to `hard_default` when additional audit work meets the dominance criteria
- Any tier can be demoted to `insufficient_data` if subsequent audit work undercuts its evidence base
- All promotions and demotions are documented in the rule's revision history within the YAML metadata

#### 5.3.2 Equivalence-finding tier criteria

§5.3.1 above assumes tool-selection audits where the finding is "tool X wins."
For audits where the finding is **equivalence, agreement, or convergence**
between tools or parameters, §5.3.1's dominance criteria ("a single tool
dominates: wins >60%") are structurally inapplicable — there is no winning
tool. The following tiers apply instead. (Adopted in response to Audit 3 CP7,
2026-05-23, which surfaced that a robust "counting tools converge" finding had
no tier under §5.3.1.)

**Tier: hard_default (equivalence)**

All of the following must hold:

- ≥8 datasets, spanning ≥2 tissues or chemistries
- Bootstrap CIs overlap across all tool pairs (B ≥ 1000, dataset-level)
- Cross-tool Spearman ρ ≥ 0.90 median, with all per-pair lower bounds ≥ 0.80
- The equivalence holds uniformly across stratifications — no
  chemistry-dependent or tissue-dependent exception to the equivalence

Rule wording: declarative — states that the tools/parameters are
interchangeable for the audited purpose; recommends choice on operational
grounds.

**Tier: conditional (equivalence)**

All of the following must hold:

- ≥6 datasets, spanning ≥2 conditions
- Bootstrap CIs overlap on most tool pairs
- Cross-tool ρ ≥ 0.80, with documented feature-dependence (the degree of
  equivalence varies by an identifiable axis: chemistry, ambient burden, etc.)

Rule wording: conditional — equivalence holds, with the feature-based
qualifier embedded.

**Tier: flag_and_warn (equivalence)**

Triggered by any of:

- Bootstrap CIs overlap substantially but cross-tool ρ < 0.80
- Equivalence holds only on a subset of stratifications
- Equivalence is directional only (e.g., a nesting/ordering holds but the
  magnitude varies materially across datasets)

Rule wording: surfaces the partial/conditional equivalence; recommends user
awareness and reporting of the tool/parameter choice.

**Tier: insufficient_data (equivalence)**

- <6 datasets for any aggregate equivalence claim.

**Nested-effect findings** (ordering universal, magnitude varies — e.g., a
permissiveness chain STAR ⊂ alevin-fry ⊂ kb): tier the **direction/ordering
claim** by the equivalence criteria above; tier any **magnitude sub-claim**
("the effect is 3× on tissue X") separately by §5.3.1, since a specific
magnitude claim about a between-tool difference resembles a selection finding.

**Contrast-based findings** (evidence from comparing two conditions — e.g.,
intestine vs PBMC control demonstrating ambient-burden dependence): the
**existence-of-effect** claim may be tiered by §5.3.2 if the contrast is
statistically robust (non-overlapping bootstrap CIs between conditions),
*even when n=1 per condition*. The **generalization** claim (does the effect
appear on other high-ambient tissues?) is tiered separately by the number of
conditions tested.

### 5.4 Rule YAML schema — full specification

The schema below is authoritative. BioOrchestrator's rule engine validates against it. Two people writing rules to this schema produce compatible YAMLs.

**Schema version:** 1.0.3

### Versioning note

Schema versions to date:

- **1.0.0** (initial) — required `last_reviewed`, `reviewer`, `revision_history` on every rule; required `out_of_scope` as string; no `status` or `severity` fields. Findings from Phase 2 rule contact: forced placeholder pollution on draft-state rules; `out_of_scope` content was naturally a list more often than a string; `level` (info/warn/error) in drafts had no schema home.
- **1.0.1** (friction amendments) — added `status` enum (gates review-related field requirements); split `created_date` from `revision_history`; `out_of_scope` accepts string OR list. Backward-compatible with 1.0.0.
- **1.0.2** (severity) — added `severity` enum distinct from `confidence_tier`. Backward-compatible with 1.0.0 and 1.0.1; default `severity: warn` when absent.
- **1.0.3** (this version, registry + reject) — adopted a canonical pipeline-step name registry (see "Pipeline-step name registry" below); added `reject_pipeline` and `block_step` to the `action_type` enum so rules whose recommendation IS to halt the pipeline (rather than change a parameter) can express that directly. Backward-compatible with 1.0.0-1.0.2: rules at older `schema_version` continue to validate; pipeline_step name compliance is a WARNING by default (errors only under `--strict-steps`) so existing rule corpora aren't broken by the registry's introduction. Resolves the v1.0.2 workaround in which rules with `severity: reject` had to carry `action_type: flag_only` cosmetically.

NOT adopted: `rule_type` enum, structured `recommendation.parameters`,
expanded `condition_type` enum, `literature_note` / `data_table` rule
types. The Phase 2 rule-categorization pass (RULE_CATEGORIZATION.md)
established that 4 of the 11 artifacts in `draft_rules/` were never
actually rules — they were findings, data tables, or coding standards
that should live elsewhere. Once those non-rules are moved out, the
schema does not need the breadth that v1.1.0 would have introduced.
The schema stays narrow; the rule corpus stays focused on actionable
pipeline rules.

**Required top-level fields:**

```yaml
rule_id: <string>                   # snake_case slug, unique across corpus
schema_version: "1.0.3"              # this schema's version
status: <enum>                      # draft | review_ready | reviewed | deployed | deprecated
title: <string>                     # human-readable, ≤80 chars
description: <string>               # 1-3 sentences, what the rule does
trigger_conditions:                 # see below
  - <condition object>
recommendation:                     # see below
  text: <string>
  action_type: <enum>
  code_example: <string, optional>
confidence_tier: <enum>             # one of: hard_default, conditional, flag_and_warn, literature_based, insufficient_data
severity: <enum>                    # one of: info, warn, error, reject. Default: warn if absent. Distinct from confidence_tier.
prior_audit_relationship: <enum>    # one of: as_original, refines_prior, contradicts_prior, extends_prior, novel
evidence:                           # see below
  - <evidence object>
out_of_scope: <string or list of strings>  # explicit statement of what was not tested
created_date: <ISO 8601 date>       # when the rule was first drafted
```

**Required at `status: review_ready` or later (not required at `draft`):**

```yaml
last_reviewed: <ISO 8601 date>      # required at review_ready+
reviewer: <string>                  # required at review_ready+
```

**Required at `status: reviewed` or later, when the rule has been revised after creation:**

```yaml
revision_history:                   # required only on second-and-later revisions
  - <revision object>
```

**Optional top-level fields:**

```yaml
mechanism_notes: <string>           # hedged speculation, never declarative
applicability:                      # narrows when the rule fires
  modalities: [<enum>...]           # e.g., ["scRNA-seq", "bulk_RNA-seq"]
  tissues: [<string>...]            # specific tissues if applicable
  organisms: [<enum>...]            # ["human", "mouse", "rat", ...]
related_rules: [<rule_id>...]       # cross-references
deprecation:                        # only present if rule is deprecated
  deprecated_date: <ISO 8601 date>
  superseded_by: <rule_id>
  reason: <string>
```

**Field constraints and enumerations:**

`rule_id`:
- Regex: `^[a-z][a-z0-9_]{2,63}$`
- Must be unique within the rule corpus
- Examples: `bioc_version_sensitivity`, `mito_threshold_quantile`, `sva_preprocessing_sensitivity`

`schema_version`:
- Semantic version of this schema
- Rules using an older schema version must be migrated before integration

`title`:
- Plain text, ≤ 80 characters
- No markdown, no quotes

`description`:
- 1-3 sentences, plain text or limited markdown (no headers, no code blocks)
- States what the rule does and when it fires

`trigger_conditions`:
- Array of condition objects, ≥ 1 required
- Each condition object:
  ```yaml
  - condition_type: <enum>          # see condition types below
    parameters: <object>            # type-specific parameters
    description: <string>           # human-readable explanation
  ```
- Condition types:
  - `dataset_feature_threshold` — fires when a dataset characteristic crosses a threshold (e.g., n_donors < 5)
  - `tool_invocation` — fires when a specific tool is being used
  - `tool_version_constraint` — fires when a tool version is outside a specified range
  - `pipeline_step` — fires at a specific pipeline step (e.g., normalization, clustering, DE)
  - `parameter_value_check` — fires when a parameter is set to a specific value or not pinned
  - `analysis_context` — fires when the analysis matches a context (e.g., longitudinal comparison)
  - `composite` — logical combination of other conditions (AND, OR)

`recommendation`:
- Required object with `text` and `action_type`
- `text`: actionable, specifies how not just what, plain text
- `action_type` enum: `replace_parameter`, `add_step`, `remove_step`, `flag_only`, `require_documentation`, `compute_from_data`, `pin_version`, `report_additional_metric`, `reject_pipeline` (v1.0.3+), `block_step` (v1.0.3+)
- `reject_pipeline` — the rule's recommendation IS to halt the pipeline. Used when no parameter change can salvage the configuration (e.g., a monotonic-by-construction clustering metric used standalone). Cross-field validator [xref-6] requires paired `severity: reject`. Resolves the v1.0.2 workaround in which `reject`-severity rules carried `action_type: flag_only` cosmetically.
- `block_step` — the rule blocks the current pipeline step only; upstream and downstream steps unaffected. Used when a single step's configuration is invalid but the rest of the pipeline is recoverable by changing only that step. Cross-field validator [xref-7] requires `severity ∈ {error, reject}`.
- `code_example`: optional, a snippet showing how to apply the recommendation

`confidence_tier`:
- Enum: `hard_default`, `conditional`, `flag_and_warn`, `literature_based`, `insufficient_data`
- Boundaries defined in Section 5.3
- Cross-field validation: see below

`severity`:
- Enum: `info`, `warn`, `error`, `reject`. Default: `warn` if absent.
- Distinct from `confidence_tier`. `confidence_tier` is evidence strength; `severity` is user-impact urgency. They are orthogonal axes.
- `info` — informational, no action expected
- `warn` — flag and recommend, do not block
- `error` — block pipeline unless explicitly overridden
- `reject` — hard reject, pipeline cannot proceed
- A `flag_and_warn` confidence_tier rule may carry any severity; a `hard_default` rule typically carries `warn` or `error`; an `insufficient_data` rule typically carries `info` or `warn`.
- **Severity is authoritative when severity and `recommendation.action_type` imply different engine behaviors.** `severity` drives whether the pipeline proceeds; `action_type` describes how the rule's recommendation should be presented to the user (e.g., `replace_parameter` describes the user-facing form of the suggested change). When a rule's intent is to hard-reject but `action_type` lacks a matching value (current example: schema v1.0.2 has no `reject_pipeline` in `action_type`), use `severity: reject` plus the closest applicable `action_type` (commonly `flag_only`) and document the intent in the rule's title/description. See `phase2/MIGRATION_NOTES.md` Phase 2b candidate B2b-2 for the schema amendment that would close this gap.

`status`:
- Enum: `draft`, `review_ready`, `reviewed`, `deployed`, `deprecated`
- Required on every rule.
- BioOrchestrator's rule engine only loads rules at `status: deployed`. `reviewed` rules are integration-ready but not yet live.
- Field requirements gate on status — see "Required at status:..." sections above.

`created_date`:
- Required on every rule.
- ISO 8601 date (YYYY-MM-DD)
- Set once at rule creation; never modified.

`prior_audit_relationship`:
- Enum: `as_original`, `refines_prior`, `contradicts_prior`, `extends_prior`, `novel`
- `novel` only for rules that don't touch any prior audit finding (e.g., entirely new modality)
- Cross-field validation: if `refines_prior` or `contradicts_prior`, then `evidence` must include both the prior audit's evidence and the new audit's evidence

`evidence`:
- Array of evidence objects, ≥ 1 required (≥ 2 required if `prior_audit_relationship` is `refines_prior` or `contradicts_prior`)
- Each evidence object:
  ```yaml
  - audit_id: <string>              # e.g., "P4", "A7", "B4_native"
    audit_phase: <string>           # e.g., "Phase 1", "Phase 2a"
    output_paths: [<relative path>...]  # paths relative to audit root
    lock_file_entries: [<sha256>...]    # hashes from the audit's lock file
    findings_md_section: <string>   # anchor in findings.md, e.g., "phase2a/findings.md#a7-sva-preprocessing-sensitivity"
    summary: <string>               # 1-2 sentence summary of what this evidence shows
  ```

`out_of_scope`:
- Required, plain text
- Explicitly states what was not tested
- Future-readers should understand the rule's boundaries

`last_reviewed`:
- ISO 8601 date (YYYY-MM-DD)
- Updated every time the rule is reviewed (not just when changed)

`reviewer`:
- Plain text, human reviewer name
- Auto mode is never the reviewer

`revision_history`:
- Array of revision objects, ≥ 1 required (first entry is the initial creation)
- Each revision object:
  ```yaml
  - date: <ISO 8601 date>
    reviewer: <string>
    change_type: <enum>             # one of: create, refine_evidence, change_tier, change_recommendation, deprecate, restore
    summary: <string>               # what changed and why
    previous_tier: <enum, optional> # if change_type is change_tier
    previous_recommendation: <string, optional>
  ```

`mechanism_notes`:
- Optional, plain text
- Always hedged: "most likely candidate per documentation review, not isolated experimentally"
- Never used as the basis for recommendation; recommendation stays on the observation

`applicability`:
- Optional object
- Narrows when the rule fires
- If absent, the rule applies whenever trigger_conditions match
- If present, the rule only fires when both trigger_conditions match AND the analysis falls within the applicability scope

`related_rules`:
- Optional array of rule_ids
- For cross-referencing rules that interact

`deprecation`:
- Only present if the rule is no longer current
- Required fields if present: `deprecated_date`, `reason`
- `superseded_by`: optional, references the replacement rule

**Cross-field validation rules:**

1. If `confidence_tier: hard_default`, then evidence must include ≥1 audit entry with ≥15 datasets covering ≥3 tissues
2. If `confidence_tier: conditional`, then evidence must include ≥1 audit entry with ≥10 datasets and the conditional logic must be expressible in `trigger_conditions`
3. If `confidence_tier: insufficient_data`, then `recommendation.action_type` must not be `replace_parameter` with a specific value; must be `compute_from_data`, `flag_only`, or `require_documentation`
4. If `prior_audit_relationship: refines_prior` or `contradicts_prior`, then evidence must contain ≥2 entries: the prior audit's and the refining/contradicting audit's
5. If `deprecation` is present, the rule is not loaded by the BioOrchestrator engine but is retained for provenance
6. (v1.0.3+) If `recommendation.action_type: reject_pipeline`, then `severity` must be `reject`. Enforced as a hard error by the validator.
7. (v1.0.3+) If `recommendation.action_type: block_step`, then `severity` must be `error` or `reject`. Enforced as a hard error by the validator.

**Pipeline-step name registry (v1.0.3+):**

`trigger_conditions[].condition_type: pipeline_step` references named pipeline steps via `parameters.step`. Step names are drawn from a registry to prevent naming drift across the rule corpus (two rules using slightly different spellings for the same step both fire on different inputs without the user noticing).

The registry's authoritative source is this section. The validator's machine-readable mirror lives at `standards/pipeline_step_registry.yaml` and is loaded at validate time.

Initial registry (2026-05-17):

| Step name | Description | Modalities | Parameters typically referenced |
|---|---|---|---|
| `pathway_enrichment` | Pathway enrichment / over-representation analysis. Covers both GSEA and ORA paradigms; operates on DEG lists (fgsea, gseapy.enrichr, clusterProfiler::enricher) or expression matrices (GSVA, EGSEA, camera). | bulk_rnaseq, scRNA-seq DEG lists, cross-modal | enrichment_paradigm, enrichment_tool, database, background, multiple_testing_correction |
| `scrnaseq_qc_filtering` | scRNA-seq quality control / cell filtering: mitochondrial fraction thresholding, library size bounds, gene count bounds, doublet removal. | scRNA-seq | mito_threshold, min_counts, max_counts, n_genes_threshold |
| `scrnaseq_clustering_resolution_selection` | scRNA-seq clustering resolution choice via selection metric (ARI, homogeneity, completeness, V-measure, silhouette). Decides which Leiden / Louvain resolution to keep across a sweep. | scRNA-seq | selection_metric, resolution, algorithm |
| `scrnaseq_de_test` | scRNA-seq differential expression test. Covers pseudobulk (muscat::pbDS), cell-level mixed models (muscat::mmDS dream, NEBULA, MAST), per-cell Wilcoxon / t-test. | scRNA-seq | method, random_effect, aggregation_function |

**Registry extension protocol:**

- Adding a new step name requires a standards-document commit; not done ad-hoc during rule drafting.
- An audit's checkpoint scope MAY propose registry additions if the audit needs to reference a step not yet registered. Additions are reviewed at audit closeout.
- The validator warns (does not error) when a rule's `pipeline_step` parameter uses an unregistered step name. Use `--strict-steps` to promote to error in CI/release-blocking contexts.
- Renaming an existing step name requires updating every rule that references it; the validator's warning surfaces the drift but does not auto-fix.

**Evidence pointer resolution:**

- `output_paths` are relative to `<audit_root>/` (e.g., `phase2a/p9_resolution_metrics_full.tsv`)
- `audit_root` is the directory containing the lock file
- BioOrchestrator's rule loader prepends the configured audit corpus root path
- If a path does not resolve, the rule is flagged as broken and not loaded
- `lock_file_entries` are SHA-256 hashes (64-char hex strings) that must exist in the audit's lock file

**Schema versioning:**

- Schema changes are versioned (1.0.0, 1.1.0, 2.0.0)
- Minor version bumps (1.0.0 → 1.1.0) are backward-compatible additions
- Major version bumps (1.0.0 → 2.0.0) require migration of all existing rules
- The rule loader supports the current major version and one prior major version during migration windows

**Example minimal rule:**

```yaml
rule_id: mito_threshold_quantile
schema_version: "1.0.0"
title: "Compute mitochondrial threshold from working dataset, not fixed default"
description: "The 20% mitochondrial threshold is tissue-dependent. Per-tissue thresholds cannot currently be encoded due to insufficient per-tissue sample size in the P1 audit. Compute the 95th percentile of mitochondrial fraction in the working dataset instead."
trigger_conditions:
  - condition_type: pipeline_step
    parameters:
      step: scrnaseq_qc_filtering
      parameter: mito_threshold
    description: "Fires when a scRNA-seq QC step applies a mitochondrial fraction threshold"
recommendation:
  text: "Compute the 95th percentile of mitochondrial fraction across the working dataset; use that as the threshold. Fall back to 20% only when working dataset has < 100 cells."
  action_type: compute_from_data
  code_example: |
    import numpy as np
    threshold = float(np.percentile(adata.obs['pct_counts_mt'], 95)) if adata.n_obs >= 100 else 0.20
confidence_tier: conditional
prior_audit_relationship: refines_prior
evidence:
  - audit_id: P1
    audit_phase: "Phase 1"
    output_paths: ["bulk_audit/output/p1_mito_threshold.tsv"]
    lock_file_entries: ["<sha256 placeholder>"]
    findings_md_section: "phase1_findings.md#p1-mito-threshold"
    summary: "20% threshold removes 0 cells in 33/60 dataset-tissue combinations; high in gut/liver, near-zero in brain/blood"
  - audit_id: A3
    audit_phase: "Phase 2a"
    output_paths: ["phase2a/p1_per_tissue_thresholds.tsv"]
    lock_file_entries: ["<sha256 placeholder>"]
    findings_md_section: "phase2a/findings.md#a3-per-tissue-mito-threshold"
    summary: "P1's n=1-2 per tissue insufficient for per-tissue table; quantile-on-data preferred"
out_of_scope: "Not tested for non-RNA modalities. Not tested for tissues outside CELLxGENE Census coverage. Not validated on synthetic data."
last_reviewed: "2026-05-16"
reviewer: "Ross Meade"
revision_history:
  - date: "2026-05-16"
    reviewer: "Ross Meade"
    change_type: create
    summary: "Initial draft from Phase 2a closeout amendments"
applicability:
  modalities: ["scRNA-seq"]
  organisms: ["human", "mouse"]
```

---

## 6. Documentation structure

Each audit produces, at minimum:

### 6.1 Per-audit documents

- `STATUS.md` at the audit root — high-level summary of what the audit did and produced
- `findings.md` per sub-audit or major finding — detailed results, sample sizes, metrics, limitations
- `draft_rules/<rule_name>.yaml` — one file per draft rule, awaiting supervised review
- `datasets/<audit_inputs>.tsv` — inventory of every input dataset with provenance
- `closeout_amendments.md` — produced by the closeout pass when amendments are triggered

### 6.2 Per-phase documents

- A comprehensive summary document at the phase level (the PHASE2_FULL_SUMMARY pattern)
- The summary lists every audit, every finding, every draft rule, with explicit relationship tags
- The summary is updated as the phase progresses, not written only at the end

### 6.3 Closeout amendment documents

When a phase undergoes a closeout amendments pass:

- A separate document captures the amendments per Section 4.3
- Each amendment is numbered and linked to the trigger that fired
- The amendments are integrated into findings.md and STATUS.md after they are documented separately

---

## 7. Session-level discipline

### 7.1 Standards loaded at session start

Every audit's local CLAUDE.md references this document. Auto mode drifts from standards over long sessions; loading them at session start mitigates drift.

### 7.2 Stop and report patterns

Auto mode is acceptable for mechanical phases (data generation, hash registration, metric computation). Supervised mode is required for interpretive phases:

- Mechanism identification
- Rule encoding
- Closeout amendments
- Integration into BioOrchestrator

Audits stop and report at the boundaries between these phases. Auto mode does not bridge a mechanical phase into an interpretive phase without human review.

### 7.3 No silent failures

If a tool install fails, a dataset cannot be loaded, or a metric cannot be computed: surface it loudly in the status report. Do not silently drop and proceed with reduced scope.

### 7.4 Compute discipline

- Local desktop preferred; AWS only when local cannot handle the workload
- AWS use governed by `~/.claude/aws-operating-rules.md`
- Budget alerts active before any AWS resource is launched
- Spot instances default for audit workloads

---

## 8. The standards are slower

Maintaining these standards is slower than not maintaining them. The alternative is Phase 1: outputs that survived because the analyst remembered them, claims that needed retroactive refinement, sample size gaps surfaced only when someone asked.

Phase 2 and 2a produced findings defensible enough to refine prior claims rather than just rehash them, sample-size limitations surfaced as findings rather than hidden as caveats, and provenance clean enough that the audit corpus can be cited internally without anyone successfully arguing it overstates its evidence.

Worth the time.
