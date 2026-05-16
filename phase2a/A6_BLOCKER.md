# A6 BLOCKER — Phase 2a §1.9a verification cannot proceed

Status: **STOP-AND-REPORT** per Issue 6 of Phase 2a closeout.
Date discovered: 2026-05-16
Run: `a6_nebula_muscat.R` × 26 datasets, 151 min wall.

## Summary

The A6 overnight run completed all 26 datasets with `status=ok` in
`a6_summary.tsv`, but two of the three tools failed in ways that
prevent §1.9a verification:

| Tool | Datasets with output | Verdict |
|---|---|---|
| `muscat_pb_DESeq2` | 21 / 26 | OK — reasonable sig counts |
| `nebula` | 19 / 26 | OUTPUTS PATHOLOGICAL |
| `muscat_mm_dream` | **0 / 26** | **silent R failure across all** |

## Issue 1 — muscat_mm_dream produced no outputs

The `mmDS(sce, method = "dream", n_threads = 2, verbose = FALSE)`
block in `a6_nebula_muscat.R` lines 122-138 (approx) is wrapped in
`tryCatch`. The orchestrator (`a6_nebula_muscat_run.py`) only captures
R stdout/stderr when subprocess returncode != 0; for all 26 datasets
returncode was 0, so we have no record of the error.

Per-dataset runtimes in `a6_summary.tsv` range from 22–260 sec.
muscat::mmDS with method="dream" on real datasets typically takes
10+ minutes per dataset, so the speed alone suggests mmDS never
actually executed the full fit.

**Net:** the §1.9a bimodal claim — that muscat-dream was inflated on
lymph node and skin fibroblast but conservative elsewhere — cannot be
verified from this run. There is no muscat_mm_dream output to inspect.

## Issue 2 — NEBULA outputs are pathological

For the 19 datasets where NEBULA wrote output, the results suggest
degenerate model fits, not real biological signal:

| Example dataset | NEBULA sig (padj<0.05) | NEBULA \|log2FC\|>5 | NEBULA p=0 |
|---|---|---|---|
| 218acb0f | 9479 / 9482 (99.97%) | 8818 | 9394 |
| 4a5b00e0 | 9992 / 10000 (99.92%) | 7868 | 9751 |
| 67b6b9ac | 10000 / 10000 (100%) | 7327 | 9647 |

In all three: nearly every gene is "significant" with extreme effect
sizes and exact-zero p-values. This is not biology; it looks like
NEBULA's negative-binomial mixed model is failing to estimate the
variance component (likely because most P4 datasets have ≤6 donors),
collapsing to a fixed-effect fit on cell-level data and reporting
pseudoreplication-inflated test statistics.

**Cross-check with §1.9a 5-dataset claim:**
- §1.9a 5-dataset subset: NEBULA inflation 0.85× vs pseudobulk (slightly conservative)
- 15-dataset full check (this run): median **5.5×**, mean **190×**, max **2499×**
- top-100 Jaccard NEBULA-vs-pseudobulk: median **0.000**, max 0.015
- log2FC Spearman NEBULA-vs-pseudobulk: median **0.021** (no correlation)

These three metrics together say NEBULA in this run is computing
something that bears no relationship to the pseudobulk DESeq2 reference.

It is unclear whether this is:
(a) a real audit finding — NEBULA mis-applied is inflated, which is
    exactly the kind of thing the audit corpus should encode; or
(b) a configuration bug in our R script — `as.numeric(donor)` or
    the offset spec or method="LN" produces this only because of how
    we set it up.

## Issue 3 — §1.9a verification is blocked

Issue 6 of Phase 2a closeout said:

> If the full 21 results contradict the original §1.9a framing,
> surface this immediately. Don't proceed with rule encoding until
> you've reported the contradiction and we've decided how to handle it.

This is that moment. Neither tool gives us trustworthy 21-dataset
metrics for the bimodal claim. The pseudobulk reference is fine but
it's not what §1.9a is about — §1.9a was specifically about mixed
models on cell-level data.

## Options for user direction

1. **Debug NEBULA configuration** — re-examine `a6_nebula_muscat.R`
   lines 68-93. Try `method = "HL"` instead of `"LN"`, or add explicit
   gene-level filtering, or use `nbglmm` directly. Re-run NEBULA only
   on the 21 P4 datasets. Cost: ~1 hour scripting + 1-2 hour run.

2. **Debug muscat_mm_dream silent failure** — add explicit stderr
   pipe-through and error logging in the orchestrator; re-run mmDS
   only on 5 datasets first to surface the real error. Cost: ~30 min
   scripting + 30-60 min run (dream is slow).

3. **Skip mixed-model arm entirely** — record §1.9a verification as
   inconclusive in the Phase 2a closeout. Do not encode a NEBULA or
   muscat-dream rule. The pseudobulk arm (muscat_pb_DESeq2) is fine
   and aligns with Phase 1 P4. Encode that part if needed.

4. **Contradict prior** — record the 5.5× NEBULA inflation finding as
   `prior_audit_relationship: contradicts_prior` with the caveat that
   the contradiction may itself be artifactual. Risky — we're flagging
   a finding that we ourselves don't trust.

My recommendation: **option 1 + option 2 in sequence**. The NEBULA
result is suspicious enough to warrant a configuration debug before
calling it an audit finding, and muscat_mm_dream's silent failure
needs to be understood before we can claim §1.9a is or isn't
contradicted. Total cost ~3-4 hours.

## What is not blocked

The other items in Phase 2a closeout (A1-A5, A7) are already complete
and not affected by this. The BioOrchestrator v0.3 implementation
work is paused at Checkpoint 4/5 boundary and also not affected — it
does not depend on A6.

## Files

- Outputs: `/mnt/nvme1/omics-audit/phase2a/a6_mixed_model/<short>/`
- Analysis: `/mnt/nvme1/omics-audit/phase2a/a6_analysis.tsv`
- Script: `/mnt/nvme1/omics-audit/phase2/scripts/a6_nebula_muscat.R`
- Orchestrator: `/mnt/nvme1/omics-audit/phase2/scripts/a6_nebula_muscat_run.py`
- Summary: `/mnt/nvme1/omics-audit/phase2a/a6_mixed_model/a6_summary.tsv`
- Log: most recent A6 log (Python orchestrator only — R stdout/stderr lost on success)
