# Audit 3 CP1 — Dataset feasibility analysis

**Date:** 2026-05-17 (revised same-day after user direction)
**Inventory file:** `audit3_inputs.tsv` (11 datasets)
**Revision history:**
- v1 (initial): 10 datasets including Mereu HCA benchmark
- v2 (this version, applied user direction 2026-05-17): Mereu removed,
  2 mouse datasets added (heart_1k_v2, tabula_muris_liver_droplet);
  added `chemistry_exact` column for v3 vs v3.1 stratification

## Working set summary

| Chemistry | n | n cells | est. size | mouse / human |
|---|---:|---:|---:|---|
| 3' v3 (v3 + v3.1)  | 4 | ~18,778 |  ~115 GB | 1 mouse / 3 human |
| 3' v2              | 4 |  ~9,051 |   ~74 GB | 2 mouse / 2 human |
| 5' v2              | 2 | ~15,000 |  ~105 GB | 0 mouse / 2 human |
| multiome RNA       | 1 | ~10,970 |   ~50 GB | 0 mouse / 1 human |
| **Total**          | **11** | **~53.8k** | **~344 GB** | **3 mouse / 8 human** |

(Note: user direction expected "12 datasets, 9 human + 3 mouse";
arithmetic gives 11 with my Mereu-removed +2 mouse additions. The
8-vs-9 human discrepancy was flagged at delivery — likely the user
double-counted neurons_900_v3 which is mouse, not human. Inventory is
8 human + 3 mouse = 11 as built; adjust if you want a 12th dataset
added.)

All 11 at `acquisition_status: needs_download` per CP0.

## v3 vs v3.1 stratification (per user direction)

New `chemistry_exact` column captures exact 10x sub-version while
`chemistry` retains the broad category for stratification:

| chemistry | chemistry_exact | n |
|---|---|---:|
| 3p_v3      | v3      | 2 (pbmc_1k_v3, neurons_900_v3) |
| 3p_v3      | v3.1    | 2 (pbmc_5k_v3.1, pbmc_10k_v3.1) |
| 3p_v2      | v2      | 4 |
| 5p_v2      | 5p_v2   | 2 |
| multiome_rna | multiome_v1.0 | 1 |

Analysis can stratify by `chemistry_exact` to test whether v3 → v3.1
chemistry changes affect counting tool behavior — finding either way
is informative.

## Mouse coverage (3 datasets, 3 distinct tissues, 2 chemistries)

| dataset | tissue | chemistry | n_cells | size |
|---|---|---|---:|---:|
| neurons_900_v3                 | brain | 3p_v3 (v3) | 931  |  ~5 GB |
| heart_1k_v2_mouse              | heart | 3p_v2 (v2) | 1011 |  ~6 GB |
| tabula_muris_liver_droplet     | liver | 3p_v2 (v2) | 1500 | ~25 GB |

n=3 spans 3 tissues and 2 chemistries. Meets user's "n=3 starts to show
whether patterns hold across mouse data" floor. Does NOT prove
generalizability — that's queued as Audit 3b in DEFERRED.md if C1/C2
findings show species-dependent counting tool behavior.

## Sample size adequacy per AUDIT_STANDARDS.md §5.3

| Confidence tier | Requirement | Status |
|---|---|---|
| `hard_default`     | ≥15 datasets, ≥3 tissues, +dominance | **NOT MET** (11 datasets; 4 tissues but PBMC dominates 8/11) |
| `conditional`      | ≥10 datasets, dependency expressible in trigger_conditions | **MET** (11 ≥ 10; chemistry-stratifiable) |
| `flag_and_warn`    | substantial disagreement; no size floor | fallback |
| `insufficient_data`| no size floor | last resort |

Tissue diversity increased from 2 → 4 distinct (PBMC, brain, heart,
liver) but PBMC dominance (8/11) still binds. Audit-3-derived rules cap
at `conditional` tier; will be surfaced loudly in CP7 per §3.1.

## Storage feasibility (revised totals)

| Resource | Est. size |
|---|---:|
| FASTQs (working set, all 11 datasets at once) | ~344 GB |
| Tool references (GRCh38 + mm39 + 4 tool indexes; mm39 partly cached) | ~80-100 GB |
| STARsolo transient files (per-dataset, peak) | ~50-100 GB |
| Tool installations | ~10-20 GB |
| Per-tool outputs (4 × 11 datasets, count matrices + metadata) | ~5-10 GB |
| **Total max projected** | **~550 GB** |
| **Free across nvme1 + nvme2 + /** | **~7.1 TB** |
| **Working set vs free space** | **~8%** |

Working set dropped from CP1-v1's ~700 GB to ~550 GB after Mereu
removal (-150 GB) + 2 small mouse adds (+30 GB). **Phased download
NOT required.**

## Disqualifying access barriers

None. All 11 candidates are open-access (10x public CC-BY-4.0 or
GEO/CZ Biohub open release for Tabula Muris).

## Chemistry coverage notes

- **3' v3 (v3 + v3.1) — 4 datasets:** v3 (pbmc_1k_v3,
  neurons_900_v3) and v3.1 (pbmc_5k_v3.1, pbmc_10k_v3.1).
  `chemistry_exact` column lets analysis stratify if v3 → v3.1
  matters; broad category enables pooled analysis if it doesn't.
- **3' v2 — 4 datasets (best-replicated chemistry):** pbmc_3k_v2,
  pbmc_4k_v2 (human), heart_1k_v2_mouse, tabula_muris_liver_droplet
  (mouse). Within-chemistry replication n=4 enables the cleanest
  tool-pair variance estimate.
- **5' v2 — 2 datasets:** sc5p_v2_hs_PBMC_5k, sc5p_v2_hs_PBMC_10k.
  Minimum n=2 for within-chemistry replication.
- **Multiome RNA — 1 dataset:** pbmc_10k_multiome (GEX library only).
  n=1 means findings here can only be flag_and_warn, not conditional,
  unless another multiome dataset is added. Flagged in CP7 limitations.

## Tabula Muris liver subset — implementation note

`tabula_muris_liver_droplet` references GSE109774 (the parent series for
all of Tabula Muris). The specific GSM accessions for the liver-droplet
samples must be confirmed at download time by parsing GEO for samples
where `source_name` includes 'Liver' AND `library_strategy_protocol`
indicates 10x Chromium droplet (vs FACS/Smart-seq2). The Tabula Muris
paper notes 3 mice per droplet-tissue combination; likely 3 GSMs total
for liver droplet. Cell count and size are estimated from the Tabula
Muris paper Table 1; precise figures captured at CP3.

## Reference reuse from CP0

| Reference | Path | Saves |
|---|---|---|
| mouse mm39 (GRCm39 + GENCODE vM34) | `/mnt/nvme2/refs/mm39/` | ~30 GB download + index build for all 3 mouse datasets |
| mouse mm39 salmon index           | `/mnt/nvme2/refs/mm39/salmon_index/` | alevin-fry index build |
| mouse mm39 STAR index (sjdb=150)  | `/mnt/nvme2/refs/mm39/star_index_sjdb150/` | STARsolo index build (verify sjdb matches read length per-dataset) |

Caveat: existing cached reference is `mm39` (GRCm39); the inventory
lists `reference_genome_expected: mm10` for all 3 mouse datasets (the
original counter version each submitter used). Decision needed at CP2:
either (a) build matching mm10 references from scratch to faithfully
re-create the submitters' configuration, or (b) reuse the cached mm39
and treat the genome-version difference as one of the experimental
variables. Recommend (b) for parsimony, with a note in CP4 findings
documenting the reference is mm39 not mm10.

Human reference (GRCh38) will need to be downloaded and indexed for
each of the four tools in CP2. ~30 GB FASTA + GTF + 4 × index = ~80-100
GB total for human references.

## CellRanger access — flagged risk (CP2)

Unchanged from v1: requires registration + license agreement at
`https://www.10xgenomics.com/support/software/cell-ranger/downloads`.
Free for academic use but gated. CP2 will exercise the download path
before assuming access.

## Limitations to surface in CP7 per §3.1

Per user direction, surface in CP7:

1. **Human tissue concentration in PBMC** — 7 of 8 human datasets are
   PBMC; only the 1 multiome PBMC sample varies any structural
   parameter from the others (nuclei vs whole cells). Rules derived
   from human data are PBMC-specific without out-of-domain validation.
2. **Mouse n=3 is a floor, not a generalizability proof** — Spans 3
   tissues but 2 of 3 are 3' v2 (one is 3' v3). Single-dataset chemistry
   slots for mouse mean species × chemistry interactions are
   under-powered.
3. **Audit 3 rules cap at `conditional` tier** — no `hard_default`
   tier possible without ≥15 datasets across ≥3 tissues (the 4 tissues
   are present but the count is below 15).
4. **Multiome n=1 is single-flag-only territory** — multiome findings
   can be flag_and_warn at best; not conditional.
5. **`chemistry_exact` v3 vs v3.1 is a measurement, not a hypothesis** —
   no prior expectation for direction or magnitude of difference; if
   the analysis finds a difference, that's a novel finding worth a
   rule; if not, that's also a publishable claim of
   sub-version-invariance.

## Deferred follow-up audit

If C1/C2 findings show species-dependent counting tool behavior,
**Audit 3b** is queued in `DEFERRED.md`:
- Scope: mouse expansion to species-symmetric 9 + 9 (9 mouse + 9 human,
  matched tissues where possible)
- Sources: Tabula Muris (mouse: kidney, lung, marrow, spleen, etc.) +
  HCA + 10x demo mouse extras
- Trigger to run: Audit 3 C1 or C2 finding that mouse vs human tool
  behavior differs at α=0.05 in any tool-pair metric

## Deliverables (revised)

- `inventory/audit3_inputs.tsv` (11 datasets × 21 columns; v2)
- `inventory/feasibility.md` (this file, v2)
- `inventory/_rebuild_audit3_inputs.py` (script that produced v2 TSV;
  re-runnable)
- `inventory/_build_inventory.py` (CP0 FASTQ scan builder, unchanged)
- `inventory/local_fastq_inventory.tsv` (CP0, unchanged)
- `inventory/local_fastq_summary.md` (CP0 + addendum, unchanged)

## Open questions remaining for CP2

(All from v1 except #1-3 + #5 resolved by user direction)

- ~~Mereu inclusion~~ → **resolved: skip**
- ~~Mouse coverage~~ → **resolved: +2 mouse, 3' v2 priority**
- ~~Mereu reference~~ → **moot**
- ~~v3 vs v3.1 stratification~~ → **resolved: stratify (chemistry_exact column)**
- ~~5' v2 only~~ → **resolved: accept**

**New questions surfaced by v2:**

1. **mm39 vs mm10 mouse reference** (see "Reference reuse from CP0").
   Recommend reuse cached mm39 + document; user override possible.
2. **Inventory size discrepancy** (8 human vs user's 9 human note).
   Flag for confirmation before CP2 starts downloads. Either accept
   11 datasets as-is, or add a 12th (probably another mouse-tissue
   3' v3 dataset to keep chemistry balance) before CP2.
