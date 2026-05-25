# Audit Ambient Correction — scRNA-seq Ambient RNA Correction Methods

_Phase 2. CP0–CP4 complete 2026-05-25. Standards: AUDIT_STANDARDS.md v1.0.3 + §5.3.2._

## Question
Do ambient RNA correction tools (SoupX, CellBender, DecontX) produce
equivalent corrections at default settings, and does the choice (tool
or ordering) propagate to biological conclusions?

## Working set
9 datasets from Audit 3 (cross-audit consistency), 2 chemistries (3' v2,
3' v3), human + mouse, multiple tissues; plus intestine + PBMC for
biological propagation. STARsolo held constant as the upstream counting
tool (Audit 3 reference). All ambient tools run at default settings.

## Finding 1: Tools produce non-equivalent per-gene contamination estimates (CP1)
Cross-tool Spearman on per-gene contamination is moderate at best (means
0.57 SoupX↔CellBender, 0.41 CellBender↔DecontX, 0.39 SoupX↔DecontX; all 27
dataset×pair values 0.24–0.69). Fails §5.3.2 equivalence (needs ρ≥0.90).
Three distinct, architecture-driven behaviors:
- CellBender tracks ambient burden (~1% PBMC → 35% intestine)
- DecontX systematically most aggressive (up to ~26%)
- SoupX floors at ~1% and fails to detect high ambient (intestine 1.5%)

Tier: non-equivalence finding, hard_default by inverse §5.3.2 criteria
(9 datasets, 2 chemistries, robust cross-tool disagreement, ambient-burden-
dependent magnitude, tight pooled CIs). Rule: `scrna_ambient_correction_tool_non_equivalence`.

## Finding 2: Ordering effects are tool-specific (CP2)
Only DecontX is ordering-sensitive (per-gene contamination shifts when fit
on QC-passed cells vs all cells; pooled ρ 0.81 v2 / 0.86 v3, as low as 0.38
on kidney). SoupX is ordering-invariant (ρ ~0.99). CellBender's correction is
ordering-invariant by design (requires the raw droplet distribution; cannot
be refit on a QC-passed subset). Maps cleanly onto tool architecture: raw-
droplet methods (CellBender) and stable-soup methods (SoupX) don't care about
order; cell-matrix mixture methods (DecontX) do.

Tier: hard_default (universal pattern across 9 datasets, both chemistries,
mechanistically explained). Rule: `scrna_ambient_correction_decontx_ordering_sensitivity`.

## Finding 3: "Correct → QC" is stricter than "QC → correct" (CP2)
Universal effect across all tools and datasets. Correction lowers counts
below the QC floor, so correct→QC removes more cells than QC→correct. The
gap scales with ambient burden (CellBender drops 2,347 more cells on the lung
organoid, 823 more on kidney, under correct→QC).

Tier: hard_default (universal across 9 datasets, both chemistries,
mechanistically explained). Rule: `scrna_ambient_correction_correct_then_qc_stricter`.

## Finding 4: Tool choice and ordering propagate to biology on high-ambient tissue (CP3)
PBMC (clean control): ARI 0.85–0.90 vs no-correction baseline across all 5
corrected conditions; cross-tool ARI 0.78–0.89; ordering effects negligible.
Near-null behavior. Intestine (high-ambient): ARI 0.50–0.70 vs baseline;
cross-tool ARI 0.46–0.60; DecontX O1-vs-baseline 0.57 vs O2 0.66 (≈10-point
gap from the negligible PBMC case); CellBender most disruptive (ARI 0.50, and
its 35% removal + correct→QC drops 23% of cells). Contested cells reassign
within related subtypes (sub-cluster granularity), not across major lineages.

Tier: existence-of-effect hard_default (PBMC-vs-intestine contrast robust,
tight non-overlapping CIs, mirrors Audit 3 C3 methodology); generalization
flag_and_warn (n=1 high-ambient + n=1 control). Rule:
`scrna_ambient_correction_high_ambient_biology_propagation`.

## Synthesis
The audit's arc:
- Tools produce non-equivalent corrections (CP1, Finding 1)
- One tool (DecontX) is ordering-sensitive while others aren't (CP2, Finding 2)
- "Correct → QC" universally removes more cells (CP2, Finding 3)
- These technical differences propagate to biological conclusions on high-
  ambient tissue and wash out on clean tissue (CP3, Finding 4)

Practical implication: ambient correction tool choice is methodologically
substantive on high-ambient tissues. SoupX under-detects severely, CellBender
tracks ambient burden, DecontX is most aggressive and ordering-sensitive. For
high-ambient analyses, document the tool + ordering choice and acknowledge
that biological conclusions may shift with different choices; re-running a key
analysis under a second tool is advisable. For clean tissues, the choice is
effectively neutral.

## Heumos 2023 positioning
Relationship: **extends**. Heumos names SoupX and CellBender as ambient-
correction options to "consider" without empirical comparison; DecontX isn't
covered; ordering isn't addressed. This audit:
- Quantifies SoupX vs CellBender disagreement (Heumos doesn't measure)
- Adds DecontX to the comparison (Heumos doesn't cover)
- Identifies SoupX's high-ambient under-detection (Heumos doesn't anticipate)
- Documents that ordering matters specifically for DecontX (Heumos doesn't address)
- Tests biological propagation on high-ambient tissue (Heumos doesn't test)

The Heumos qualitative framing is correct directionally (ambient correction
matters) but underspecifies critical methodological choices (which tool, what
order) that this audit shows have substantive biological consequences where
ambient burden is high.

## Methodological observations
- §3.5 candidate (bootstrap CIs reflect sampling variance not correctness)
  applied as a sanity check throughout; no anomalies surfaced (magnitudes
  within published tissue expectations except SoupX's documented under-detection)
- §3.4 companion-metrics requirement satisfied in every rule (log2-ratio,
  directional agreement, stratification, retained-cell counts, ARI/NMI/marker Jaccard)
- §5.3.2 equivalence criteria applied; the CP1 non-equivalence (inverse-
  equivalence) finding tiered hard_default by symmetry
- Clean-control contrast pattern (PBMC vs intestine) reused from Audit 3 C3;
  reproducible
- Tool-architecture taxonomy (raw-droplet vs stable-soup vs cell-matrix-
  mixture) emerged as a useful conceptual frame; not formally in standards
- One registry extension required and adopted: `scrnaseq_ambient_correction`
  added to `pipeline_step_registry.yaml` (no schema change; see PENDING_AMENDMENTS.md)

## Limitations
- 9 datasets technical; 2 datasets biological propagation
- Biological propagation rests on n=1 high-ambient (intestine); needs replication
- Default parameters only; tool parameter sensitivity not tested
- STARsolo as fixed upstream; behavior with other counting tools untested
- QC-parameter choice for biological propagation (cp6 QC for C3 comparability,
  not CP2's C2 fixed-floor) discussed in cp3/tool_failure_modes.md
- CellBender HTML-report generation failures (cosmetic) noted, non-blocking;
  CellBender 0.3.2 also required a documented checkpoint patch (CP0)

## Prior audit relationship
**extends_prior.** Extends Audit 3 C3 (which established that high-ambient
tissues produce biological divergence from cell-calling choice) by:
- Testing whether ambient correction mitigates or compounds C3's finding
- Identifying ambient correction tool choice as an additional variance source
  on high-ambient tissue
- Quantifying tool-specific behaviors not previously characterized

## Audits queued in DEFERRED.md
- Ambient correction parameter sensitivity (default vs tuned)
- Tool-counting interaction (ambient correction × counting tool choice; does
  behavior differ on alevin-fry / kb_count vs STARsolo outputs?)
- Additional high-ambient tissue replication for generalization

## Rules drafted (CP4)
4 rules in `draft_rules/`, all passing the validator under `--strict-steps`:
- `scrna_ambient_correction_tool_non_equivalence` (hard_default, info)
- `scrna_ambient_correction_decontx_ordering_sensitivity` (hard_default, warn)
- `scrna_ambient_correction_correct_then_qc_stricter` (hard_default, info)
- `scrna_ambient_correction_high_ambient_biology_propagation` (flag_and_warn, warn)
