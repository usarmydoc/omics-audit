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

## Addendum (2026-05-17) — deeper scan for accessions + references

User-requested follow-up scan to catch anything CP1 might re-download or
CP2 might rebuild unnecessarily. Patterns:

- Accession-prefixed directories: `GSE*`, `GSM*`, `SRR*`, `ERR*`, `DRR*`,
  `PRJNA*`, `PRJEB*`, `PRJDB*`, `E-MTAB-*`, `E-PROT-*`, `E-GEOD-*`
- Archive formats: `*.sra`, `*.tar`, `*.tar.gz`, `*.tgz`, multi-part
  FASTQs (`*.fastq.gz.[0-9]+`)
- Reference + index files: `*.fa.gz`, `*.fasta.gz`, `*.gtf.gz`, `*.gff.gz`,
  `gencode.*.annotation.gtf*`, `genome.fa`, `*.primary_assembly.fa*`,
  `refdata-gex-*` / `refdata-cellranger-*` dirs, STAR index files
  (`SAindex`, `Genome`), salmon/kallisto index files (`*.idx`,
  `t2g.txt`), `gencode*` / `GRCh3{7,8}*` / `GRCm3{8,9}*` /
  `hg{19,38}*` / `mm{10,39}*` directories.

### Reusable artifacts found

| Path | Size (approx) | Audit 3 relevance |
|---|---:|---|
| `/mnt/nvme2/refs/mm39/GRCm39.primary_assembly.genome.fa.gz`     | ~830 MB | **CP2 reuse** — saves Gencode mouse genome download |
| `/mnt/nvme2/refs/mm39/gencode.vM34.primary_assembly.annotation.gtf{,.gz}` | ~1.4 GB | **CP2 reuse** — paired GTF for the same |
| `/mnt/nvme2/refs/mm39/gencode.vM34.transcripts.fa.gz`           | varies | **CP2 reuse** — transcriptome FASTA (alevin-fry / kallisto input) |
| `/mnt/nvme2/refs/mm39/salmon_index/`                            | varies | **CP2 reuse** — pre-built salmon index (alevin-fry compatible) |
| `/mnt/nvme2/refs/mm39/star_index_sjdb150/`                      | ~25 GB | **CP2 reuse** — pre-built STAR index (sjdb=150; valid if any mouse dataset uses 150 nt reads) |

If any Audit 3 datasets are mouse, the mm39 reference set above is the
starting point for CP2's reference building. Verify `gencode.vM34` is
acceptable for the audit (current GENCODE Mouse release at time of writing
is vM34) and that `sjdb=150` matches the dataset's read length before
reusing the STAR index. Otherwise rebuild from the cached FASTA + GTF.

### Not reusable (verified contents)

GSE RAW archives surfaced but inspected with `tar tf` to confirm they
contain processed count matrices, NOT raw FASTQs. Cannot substitute for
Audit 3 FASTQ downloads:

| Archive | Contents | Notes |
|---|---|---|
| `/mnt/nvme1/omics-audit/phase2a/gse96583/GSE96583_RAW.tar`   | `*_barcodes.tsv.gz`, `*_*.mtx.gz`/`*.mat.gz` (10 entries) | Kang et al. PBMC IFN-β. Chemistry is 10x Chromium 3' v1 — **out of Audit 3 scope** anyway (v2/v3/5p_v2/multiome only). |
| `/home/ross/scrna_GSE145197/data/raw/GSE145197_RAW.tar`      | `*_UMI_tab_*.txt.gz` (10 entries) | Droin 2021 mouse liver circadian. Processed UMI tables, not FASTQ. |
| `/mnt/nvme1/omics-audit/scrnaseq/projects/p2_doublet_audit_followup/data/GSE108313_RAW.tar` | `*_Hashtag-*`/`*-RNA.umi.txt.gz` etc. (6 entries) | Stoeckius hashtag/multiplexed; processed UMI files, not FASTQ. |
| `/mnt/nvme2/liver_scvi/data_raw/GSE216584_RAW.tar`           | `*_*_RSEC_MolsPerCell.csv.gz` (22 entries) | Hautz 2023 NMP liver. BD Rhapsody chemistry — **out of Audit 3 scope** (10x only). |

### Not relevant

The deeper scan also surfaced ~80 false positives: Windows graphics
driver cache (`D3DSCache/*.idx`), Java/Adobe cache (`AppData/Local/...`),
EMBOSS test database (`entrynam.idx`), small R/Bioconductor `extdata`
FASTAs (`someORF.fa.gz`, `dm3_upstream2000.fa.gz`), and megahit/SPAdes
test data. Not enumerated in the inventory TSV; all confidently skip.

### Storage budget verification

| Resource (max projection) | Size |
|---|---:|
| FASTQ working set (≥8 datasets × 50-200 GB) | 0.4-1.6 TB |
| Tool references (4 × 10-30 GB; minus mm39 reuse) | 0.04-0.12 TB |
| STARsolo temp files (transient, per-dataset) | ~0.1 TB |
| Tool installations (CellRanger + STAR + salmon + kallisto) | ~0.02 TB |
| **Total max projected** | **~1.85 TB** |
| **Free across nvme drives + `/`** | **~7.1 TB** |
| **Headroom** | **~3.8×** |

Phased download-and-delete batching NOT required at this scale. Proceed
to CP1 with all-at-once download budget.

## Deliverables

- `inventory/local_fastq_inventory.tsv` (68 rows × 12 cols, one row per
  discovered file from the primary FASTQ/BAM scan)
- `inventory/local_fastq_summary.md` (this file, including deeper-scan
  addendum)
- `inventory/_build_inventory.py` (the script that produced the TSV;
  re-runnable if the scan is repeated)
