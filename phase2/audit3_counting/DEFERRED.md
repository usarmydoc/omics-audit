# Audit 3 — DEFERRED scope-creep candidates

Per the audit spec: scope-creep candidates land here without action
during this work. Triaged at audit closeout or before the next audit
starts.

---

## Queued future audits (real audits, scoped, sized — NOT acted on here)

### Audit 3e — 5' v2 chemistry sub-audit if 10x access acquired

**Queued during:** CP3 user direction (2026-05-17) — both 5' v2 datasets
in the CP1 inventory (`sc5p_v2_hs_PBMC_5k`, `sc5p_v2_hs_PBMC_10k`) were
dropped because their FASTQ tarball URLs are not findable via public
URL probing. The 10x catalog pages are JavaScript-rendered SPAs and
the download links are gated behind an authenticated session.

**Trigger to run:** user obtains 10x credentials (same trigger as
Audit 3c for CellRanger). At that point all the FASTQ download URLs
for current 10x datasets become accessible.

**Proposed Audit 3e scope:**

- **Datasets:** `sc5p_v2_hs_PBMC_5k` + `sc5p_v2_hs_PBMC_10k` (the two
  Audit 3 5' v2 candidates), plus optionally a 5' v2 mouse demo if
  one exists at that point
- **Tools:** same 4 configurations as Audit 3 (STARsolo default,
  STARsolo CR-mimic, alevin-fry, kb count). 5' v2 chemistry needs
  different STARsolo parameters than 3' (R1 still has CB+UMI but
  R2 represents the 5' end of the mRNA — reverse orientation)
- **Reference:** GRCh38 + GENCODE v45 from Audit 3 CP2
- **Whitelist:** 737K-august-2016 (same as 3' v2; already downloaded)
- **Metrics:** same as C1 + C2; if 5' v2 results differ qualitatively
  from 3' results, that's the headline finding
- **Estimated effort:** 3-4 days (data acquisition + processing +
  metrics + findings; all tools already installed and validated)

**Sequencing:** Audit 3e runs only after Audit 3 is fully complete
and only when 10x credentials become available.

---

### Audit 3d — Multiome chemistry sub-audit if credentials acquired

**Queued during:** CP3 user direction (2026-05-17) — `pbmc_10k_multiome`
dataset (slot for `chemistry: multiome_rna`) dropped from the CP3
working set because the 737K-arc-v1 cell-barcode whitelist required by
all 4 counting tools is gated behind 10x's cellranger-arc tarball,
which is not publicly downloadable without a 10x account. Half-
processing the dataset (e.g., STARsolo only with `--soloCBwhitelist
None`) would produce asymmetric coverage and contaminate the C1/C2
cross-tool comparison.

**Trigger to run:** either —
- User obtains 10x credentials and pulls cellranger-arc (which bundles
  the 737K-arc-v1.txt whitelist), OR
- C1/C2 findings on 3' and 5' chemistries surface a multiome-specific
  question worth pursuing

**Why this is its own audit, not part of Audit 3:**
Per CP1 feasibility, multiome was already n=1 / flag_and_warn-only
territory — adding it after-the-fact via a 3d sub-audit doesn't lose
anything the main audit was relying on. Multiome is also somewhat
distinct from 3'/5' single-cell because the protocol couples GEX with
ATAC, which means even after the counting comparison there's a
multimodal-integration question separable from pure counting-tool
comparison.

**Proposed Audit 3d scope (locked at audit time, not now):**

- **Dataset(s):** `pbmc_10k_multiome` (the original Audit 3 dataset)
  + 1-2 more multiome datasets from 10x demo for n=2-3 chemistry
  replication
- **Tools:** same 4 configurations as Audit 3 (STARsolo default,
  STARsolo CR-mimic, alevin-fry, kb count)
- **Reference:** GRCh38 + GENCODE v45 from Audit 3 CP2; no new
  reference build needed
- **Whitelist:** 737K-arc-v1.txt from cellranger-arc bundle
- **Metrics:** same as C1 + C2 + C3; new C4 if multimodal-integration
  question arises
- **Estimated effort:** 1 week (data acquisition + processing +
  metrics + findings; tool installs already in place from CP2)

**Sequencing:** Audit 3d runs only after Audit 3 is fully complete
and only if a trigger fires.

---

### Audit 3c — CellRanger re-test if/when user obtains license

**Queued during:** CP2 user direction (2026-05-17) — user opted not to
register for CellRanger access. Audit 3 proceeds with 3 binaries / 4
tool configurations (STARsolo-default, STARsolo-CR-mimic, alevin-fry,
kallisto-bustools). STARsolo-CR-mimic partially fills the gap by
exercising the STAR-authors' documented Cell-Ranger-compatible flag set,
but is not a substitute for running CellRanger itself.

**Trigger to run:** user obtains a CellRanger license + completes the
gated download from the 10x portal. Optional secondary trigger:
Audit 3 findings show STARsolo-CR-mimic deviates substantially from
the other tools, in which case running the actual CellRanger to
adjudicate becomes more valuable.

**Why this is its own audit, not part of Audit 3:**
Adding CellRanger now would block CP2 indefinitely on the user
performing a manual registration step they have declined. Continuing
without it preserves audit momentum, and the 3-binary/4-config
comparison is itself publishable. A subsequent CellRanger-only sub-audit
can compare CellRanger's output to the existing 4 configs without
re-running the others (count matrices already produced).

**Proposed Audit 3c scope (locked at audit time, not now):**

- **Tool:** CellRanger (latest version available at audit time)
- **Datasets:** same 11 as Audit 3 (or representative subset if time-constrained)
- **Reference:** same GRCh38 + mm39 builds Audit 3 used (re-formatted
  to CellRanger's reference format via `cellranger mkref` if needed)
- **Metrics:** per-gene Spearman correlation between CellRanger counts and
  each of Audit 3's 4 configs; barcode-call Jaccard with each. Tests:
  - Does STARsolo-CR-mimic actually match CellRanger?
  - How far apart are CellRanger and the open-source alternatives?
- **Estimated effort:** 1 week (data already counted by 4 configs;
  CellRanger is the only new pipeline + the comparison metrics)

**Sequencing:** Audit 3c runs only after Audit 3 is fully complete and
only when CellRanger is available.

---

### Audit 3b — Mouse expansion to species-symmetric working set

**Queued during:** CP1 user direction (2026-05-17).

**Trigger to run:** Audit 3 C1 or C2 findings show species-dependent
counting tool behavior — i.e., at least one tool-pair metric
(Spearman ρ on per-gene counts, Jaccard on barcode calls, etc.)
differs between mouse-stratum aggregate and human-stratum aggregate
at the 95% bootstrap CI level.

**Why this is its own audit, not part of Audit 3:**
Audit 3's working set is 8 human + 3 mouse (PBMC-dominant). Mouse n=3
is sufficient to detect "species-dependent counting tool behavior" as
a signal but is far short of demonstrating generalizability across
mouse tissues / chemistries. If the C1/C2 finding for species
differences is null, Audit 3b doesn't need to run; if the finding is
real, 3b expands the mouse arm to match the human arm.

**Proposed Audit 3b scope (locked at audit time, not now):**

- **Working set:** 9 + 9 species-symmetric. Add 6 mouse datasets across
  matched tissues + chemistries to bring mouse count to 9.
- **Mouse sources:**
  - Tabula Muris GSE109774 — additional droplet tissues
    (kidney, lung, marrow, spleen, thymus, etc.)
  - 10x Genomics mouse demo datasets — additional cells/tissues
  - HCA mouse projects if needed for tissue coverage
- **Reference:** unified mouse reference build (GENCODE mouse latest
  vs mm10 vs mm39 decision lifted from Audit 3 CP2)
- **Metrics:** identical to Audit 3 (C1/C2/C3) so results are directly
  comparable
- **Estimated effort:** 2-3 weeks (data acquisition + processing +
  metrics + findings)
- **Output:** findings.md (Audit 3b section), draft rule YAMLs only if
  warranted

**Sequencing:** Audit 3b runs only after Audit 3 is fully complete
(all 8 checkpoints, all rules drafted, lock state stable) AND only if
the species-dependence trigger fires.

---

## Open questions captured but not acted on

### `chemistry_exact` v3 vs v3.1 finding disposition

If the audit's `chemistry_exact` stratification surfaces a
v3-vs-v3.1 effect:
- Magnitude small or null → finding is "sub-version-invariant," may
  fold into the broader `pathway_database_choice`-style rule for
  pathway audits (no separate rule needed)
- Magnitude substantial → emit a separate rule specifically about
  chemistry version pinning

Decision deferred until CP4 metrics land.

### mm10 vs mm39 reference for mouse datasets

Audit 3 mouse datasets were originally processed by submitters against
mm10; CP0 found a cached mm39 reference set at `/mnt/nvme2/refs/mm39/`.
CP2 will either build matching mm10 references or reuse cached mm39
+ document the reference-version difference as a known variable. The
faithful-recreation option (mm10) is the more rigorous choice for
"do counts agree across tools?"; the parsimonious option (mm39) treats
reference-version as one of the experimental knobs. Audit 3 will use
the parsimonious option per CP2 decision; if findings depend on
reference version, that's a Phase 2c candidate (separate sub-audit on
reference-version sensitivity, queued here for future triage).

---

## Hard rejections (out of scope, do not re-litigate)

### Mereu HCA benchmark (6e177195-...)

**Rejected during:** CP1 user direction (2026-05-17).

**Reason:** Adds ~150 GB download and a 4-species reference build
(GRCh38 + mm10 + canFam3 + cat) for a single dataset. Methodological
strength is real (within-sample cross-chemistry comparison) but
operational cost is disproportionate to the C1 anchor it provides.

**Re-surface trigger:** if Audit 3 C1 findings prove sensitive to
biological variance (cross-dataset variance dominates cross-tool
variance), Mereu becomes a Phase 2b candidate as the cleanest available
within-sample anchor. Otherwise stays rejected.

---

## How to use this file

- **Queued future audits** are real follow-ups with bounded scope and
  clear trigger conditions. Reviewed at Audit 3 closeout.
- **Open questions** are mid-audit decisions where the audit will pick
  the parsimonious option and re-evaluate at closeout.
- **Hard rejections** are explicitly out-of-scope decisions with
  re-surface triggers documented. Don't reopen without external trigger.
