# Audit 3 CP1 — Dataset feasibility analysis

**Date:** 2026-05-17
**Inventory file:** `audit3_inputs.tsv` (10 datasets)

## Working set summary

| Chemistry | n datasets | n cells | est. size |
|---|---:|---:|---:|
| 3' v3 / v3.1   | 4 | ~18,778 |  ~115 GB |
| 3' v2          | 2 |  ~7,040 |   ~43 GB |
| 5' v2          | 2 | ~15,000 |  ~105 GB |
| multiome RNA   | 1 | ~10,970 |   ~50 GB |
| mixed (Mereu)  | 1 (subset) | varies |  ~150 GB |
| **Total**      | **10** | **~52k** | **~463 GB** |

All datasets at `acquisition_status: needs_download` per CP0 (zero local
FASTQs after the deeper scan).

## Sample size adequacy per AUDIT_STANDARDS.md §5.3

| Confidence tier | Requirement | Status |
|---|---|---|
| `hard_default`     | ≥15 datasets, ≥3 tissues, +dominance | **NOT MET** (10 datasets, 2 tissues without Mereu's mixed) |
| `conditional`      | ≥10 datasets, dependency expressible in trigger_conditions | **MET** (10 datasets, dependencies likely chemistry-stratified) |
| `flag_and_warn`    | substantial disagreement; no size floor | **available as fallback** |
| `literature_based` | no size floor | n/a |
| `insufficient_data`| no size floor | only if findings collapse |

**Tissue diversity is the binding constraint.** PBMC dominates (7 of 10
datasets); brain (1) and Mereu's mixed reference (1) add some breadth
but the audit cannot claim tissue-level generalization. Rules drafted
from C1-C3 findings should be at most `conditional`, never
`hard_default`, until a follow-up audit covers more tissues. This
limitation will be surfaced loudly per §3.1 in CP7 synthesis.

The user-locked spec says "≥8 datasets across the four major 10x
chemistries" — that bar is met (10 ≥ 8; all four chemistries present
with ≥1 dataset each, three chemistries with ≥2 for within-chemistry
replication). The tissue under-representation is a separate concern
from the spec's literal sample-size requirement.

## Storage feasibility

| Resource | Est. size |
|---|---:|
| FASTQs (working set, all 10 datasets at once) | ~463 GB |
| Tool references (CellRanger + STAR + alevin-fry + kallisto, partly cached) | ~80-120 GB |
| STARsolo transient files (per-dataset, peak) | ~50-100 GB |
| Tool installations | ~10-20 GB |
| Per-tool outputs (4 × 10 datasets, count matrices + metadata) | ~5-10 GB |
| **Total max projected** | **~700 GB** |
| **Free across nvme1 + nvme2 + /** | **~7.1 TB** |
| **Working set vs free space** | **~10%** |

Working set is ~10% of free space — well below the spec's 50%
phased-download threshold. **No phased download-and-delete batching
required.** All 10 datasets can coexist on disk through CP3-CP6.

Per-drive recommendation:
- **`/mnt/nvme2`** (3.1 TB free) — primary FASTQ landing. Already hosts
  `refs/mm39/` so mouse reference reuse is local.
- **`/mnt/nvme1`** (3.2 TB free) — audit outputs, lock-tracked artifacts.
- **`/` (808 GB free)** — leave for OS + R/Python envs; do NOT use for
  large FASTQ landing.

## Disqualifying access barriers

None identified. All 10 candidates are open-access (10x public CC-BY-4.0
or HCA open).

## Notes on chemistry coverage

- **3' v2:** uses the classic `pbmc_3k` and `pbmc_4k` datasets. These
  are the de-facto benchmark set for v2 chemistry in published
  benchmarks; absent if not included. Within-chemistry replication (n=2)
  is the minimum for measuring tool-pair variance vs dataset-pair variance.
- **3' v3 / v3.1:** v3 (`pbmc_1k_v3`, `neurons_900_v3`) and v3.1
  (`5k_pbmc_v3_nextgem`, `10k_pbmc_v3_nextgem`) lumped under one
  chemistry category for stratification. Tools should not treat these
  as distinct, but the audit can flag if they show different
  patterns.
- **5' v2:** `sc5p_v2_hs_PBMC_5k` and `sc5p_v2_hs_PBMC_10k`. Both
  human PBMC — chemistry is what differs from 3' v2 / v3 PBMCs, useful
  for isolating chemistry effects from tissue effects.
- **Multiome RNA:** `pbmc_granulocyte_sorted_10k_arc` is the canonical
  10x multiome demo. Audit 3 uses ONLY the GEX library; the paired ATAC
  library is downloaded only if total download stays within budget,
  otherwise skipped (ATAC ≈ 50 GB additional).

## Mereu benchmark — why it's in the inventory

The Mereu et al. 2020 HCA project applied **13 scRNA-seq protocols** to
the same biological sample (5-species mixed reference: 60% human PBMC,
30% mouse colon, 6% HEK293T, 3% NIH3T3, 1% MDCK). Restricting to the 10x
subset (3' v2, 3' v3, 5' v2) gives a within-sample cross-chemistry
comparison that no other public dataset offers. This is the strongest
anchor for the C1 (per-gene count agreement) sub-audit because biological
variance is held constant; any tool disagreement is purely
methodological.

Audit 3 scope-limits the Mereu data to the 10x protocols only — the
non-10x submissions (CEL-seq2, Drop-seq, Quartz-seq2, Smart-seq2,
ddSEQ, SCRB-seq, etc.) are out of scope per the audit's
counting-tool-comparison framing (the counting tools tested don't apply
to those protocols).

## CellRanger access — flagged risk (CP2)

CellRanger requires a license-agreement download from
`https://www.10xgenomics.com/support/software/cell-ranger/downloads`
with registration. The license is free for academic use but the
download is gated.

Per the audit spec: "If CellRanger access is blocked, document and
surface before proceeding. The audit can proceed with three tools but
the scope of the C1/C2 findings changes." CP2 will exercise the
download path before assuming access is available.

## Reference reuse from CP0

| Reference | Path | Saves |
|---|---|---|
| mouse mm39 (GRCm39 + GENCODE vM34) | `/mnt/nvme2/refs/mm39/` | ~30 GB download + index build time for `neurons_900_v3` (mouse) |
| mouse mm39 salmon index           | `/mnt/nvme2/refs/mm39/salmon_index/` | alevin-fry index build |
| mouse mm39 STAR index (sjdb=150)  | `/mnt/nvme2/refs/mm39/star_index_sjdb150/` | STARsolo index build (verify sjdb matches dataset read length) |

Human reference (GRCh38) will need to be downloaded and indexed for each
of the four tools in CP2. ~30 GB FASTA + GTF + 4 × index = ~80-100 GB
total for human references.

## Phased download plan: NOT REQUIRED

Working set ≪ 50% of free space. CP3 can download all 10 datasets
upfront and retain them through CP4-CP6 analysis. Phased
download-then-delete only becomes necessary if:

- Working set grows beyond ~3.5 TB (currently ~0.5 TB → ~7× headroom)
- Per-dataset STAR temp files balloon unexpectedly
- A new tool with a larger reference footprint is added to the audit

If any of these conditions arise mid-audit, switch to phased mode and
process datasets in 2-3 batches.

## Disqualifications: None in this CP1 inventory

All 10 candidates are open-access, in-scope chemistry, and within
storage budget. No `access_status: blocked` entries.

## Deliverables

- `inventory/audit3_inputs.tsv` (10 datasets, 20 columns including
  CP0-derived `acquisition_status` and `chemistry_in_scope`)
- `inventory/feasibility.md` (this file)

## Open questions for user review before CP2

1. **Mereu inclusion:** the Mereu 10x subset is methodologically the
   strongest C1 anchor but adds ~150 GB download and one more genome
   build (canFam3 for dog MDCK). Worth the cost? Or skip Mereu and
   accept the weaker cross-chemistry comparison from disjoint 10x demo
   sets?
2. **Mouse coverage:** only 1 dedicated mouse dataset (`neurons_900_v3`).
   Add another mouse 10x demo (e.g., 1k mouse heart cells) to improve
   the species split, or accept the imbalance and stratify findings by
   species in the analysis?
3. **Reference for Mereu:** building a 4-species combined reference
   (human + mouse + dog + cat-from-MDCK = actually canine) requires
   GENCODE + Ensembl pulls for each. ~10 GB additional ref data + 4
   parallel index builds. Out of normal CP2 scope?
4. **5' v2 versus 5' v1:** the audit spec includes 5' v2 but not 5' v1.
   Both `5k_5p_v2` and `10k_5p_v2` are confirmed v2; OK.
5. **`pbmc_1k_v3` and `neurons_900_v3` use v3 (not v3.1):** stratify
   v3 from v3.1 in analysis, or treat as one chemistry? The spec lists
   "3' v3" as one chemistry — current inventory treats v3 and v3.1 as
   the same.
