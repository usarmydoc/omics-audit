# Audit 3 CP0 — Local FASTQ discovery summary

**Date:** 2026-05-17
**Scope:** Read-only scan of all three storage drives for existing scRNA-seq
FASTQ files (and BAMs convertible via bamtofastq) that could be reused as
Audit 3 inputs.

## Drives scanned

| Mount | Size | Used | Free | Filesystem |
|---|---|---|---|---|
| `/` (vgmint-root) | 1.8 TB | 931 GB | 808 GB | ext4 |
| `/mnt/nvme1` | 3.7 TB | 498 GB | 3.2 TB | NTFS |
| `/mnt/nvme2` | 3.7 TB | 582 GB | 3.1 TB | NTFS |
| **Total** | **9.2 TB** | **2.0 TB** | **7.1 TB** | |

Skip filters applied: `__pycache__`, `node_modules`, `.git`, `.cache`,
`.conda`, `.npm`, `site-packages`, `.local`, `.venv`, plus the standard
system roots (`/proc`, `/sys`, `/dev`, `/snap`, `/var`, `/tmp`, `/run`).
`/` scan used `-xdev` to avoid double-counting nvme mounts.

Patterns matched: `*.fastq.gz`, `*.fq.gz`, `*.fastq`, `*.fq`, `*.bam`
(case-insensitive).

## Totals by category

| Category | n files | Total size | Notes |
|---|---:|---:|---|
| `test_fixture`   | 26 |   <1 MB | R/Bioconductor package test data + SPAdes assembler test FASTQs |
| `user_excluded`  | 42 | 318 GB | Folders on user's 2026-05-17 exclude list |
| **Total**        | **68** | **~318 GB** | |

## Usable for Audit 3: 0 files

No discovered files are usable for Audit 3. Net result: every dataset
in the CP1 working set will need to be downloaded fresh.

## User exclusion list (2026-05-17)

The user explicitly instructed not to use files from any of these folders:

| Folder | n files | Size | Modality (per [[MEMORY]] / file shape) |
|---|---:|---:|---|
| `/mnt/nvme1/12282025Jihoon`           | 0 (none found) | — | (no FASTQs in this folder) |
| `/mnt/nvme1/EGSEA_01122026Jiuen`      | 14 | 114 GB | Bulk RNA-seq mouse (`.fastqsanger.fastq`, Galaxy export) |
| `/mnt/nvme1/miRNA_01052026Ang`        | 18 |  6 GB | miRNA-seq (trimmed, paired duplicate files) |
| `/mnt/nvme2/Agnieszka 02-20-2026`     |  2 | 94 GB | scRNA — likely Parse chemistry, awaiting confirmation |
| `/mnt/nvme2/Katherine`                |  8 | 87 GB | Visium HD spatial transcriptomics (10x naming pattern) |

(The user's exclusion list also references the same paths via the
Desktop symlinks `~/Desktop/nvme1` and `~/Desktop/nvme2` which resolve
to `/mnt/nvme1` and `/mnt/nvme2` respectively per [[MEMORY]].)

These folders are off-limits for Audit 3 work regardless of modality
match. The categorization above is documentation of what's in them, not
a basis for using them.

## Test fixtures (categorized `test_fixture`, 26 files, <1 MB total)

Library-bundled test data, not biological inputs. Will not be touched
by Audit 3.

- **R/Bioconductor BAM fixtures** (18 files): `Rhtslib/testdata/`,
  `Rsamtools/extdata/`, `Rsamtools/unitTests/cases/`,
  `GenomicAlignments/extdata/`. All KB-range synthetic data shipped
  with the packages for unit testing.
- **SPAdes assembler test FASTQs** (8 files): `share/spades/test_dataset/`
  (E. coli 1K) and `share/spades/test_dataset_plasmid/` in both
  `~/miniforge3/pkgs/spades-4.2.0.../` and
  `~/miniforge3/envs/bio-linux/share/spades/`. All KB-range.

## Implications for CP1

- **All Audit 3 datasets need fresh download.** Plan ~50-200 GB per dataset
  × ≥8 datasets = ~400 GB – 1.6 TB working set on the FASTQ side, comfortably
  within the 7.1 TB of free space across the two nvme drives.
- **Storage strategy:** processing in batches with FASTQ deletion between
  batches is NOT yet required at the projected scale; defer batching plan
  to CP1 once exact dataset sizes are inventoried.
- **No `acquisition_status: local_complete` or `local_partial` entries**
  will appear in CP1's `audit3_inputs.tsv`. All candidates start at
  `needs_download` (or `unavailable` for blocked-access datasets).

## Surprises

- Smaller than expected. The two nvme drives are at 13-16% utilization;
  no abandoned scRNA-seq downloads from prior projects were uncovered.
  The drives' content is dominated by current-project subdirectories
  (the bulk and spatial work the user excluded) plus the existing
  omics-audit + bioorchestrator + marker_atlas tree.
- No `.partial`, `.tmp`, or other in-flight download markers anywhere.
- Zero BAM-only candidate datasets (the 18 BAMs found are all R-package
  test fixtures, not biological data).

## Deliverables

- `inventory/local_fastq_inventory.tsv` (68 rows × 12 cols, one row per
  discovered file)
- `inventory/local_fastq_summary.md` (this file)
- `inventory/_build_inventory.py` (the script that produced the TSV;
  re-runnable if the scan is repeated)
