# Pending AUDIT_STANDARDS.md amendment candidates

Methodological observations surfaced during audit work that are candidates
for the next AUDIT_STANDARDS.md amendment pass. Queued here; **not**
implemented mid-audit (schema changes mid-audit are forbidden).

---

## §3.5 candidate — "Bootstrap CIs reflect sampling variance, not correctness"

**Surfaced from:** Audit 3 CP4 verification, 2026-05-22.

**Content:** Tight bootstrap CIs measure internal consistency of the sampling
distribution. They do not detect systematic errors in tool configuration,
gene ID handling, or other upstream issues that affect all bootstrap
replicates equally. Before treating tight CIs as strong evidence, verify
that the underlying data does not have systematic configuration errors.

**Audit 3 CP4 example:** bootstrap CIs were tight (B=1000, narrow intervals)
on data with an unhandled gene-ID suffix issue (alevin-fry USA-mode -U/-A
suffixes collapsed incorrectly, leaving only the ~10% ambiguous bucket as
alevin-fry's per-gene count). The CIs appeared to strongly support a
"tools disagree" headline (cross-tool rho ~0.58). After fixing the gene-ID
collapse (sum S+A per gene), the headline became "tools converge"
(rho ~0.96) — the CIs were tight in both cases, but only one was correct.

**Recommended amendment:** when an audit produces high-confidence findings
(tight CIs, low p-values), require a "sanity check" step that verifies
the underlying data is correctly configured before treating the findings
as load-bearing. This is a pre-publication check, not a methodological
change to bootstrap procedure.

**Status:** queue for next AUDIT_STANDARDS.md amendment pass; do not
implement schema changes mid-audit.
