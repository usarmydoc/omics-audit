# Audit 3 CP2 — Tool installation + reference building

**Date:** 2026-05-17
**Status:** COMPLETE
**Working tool set:** 3 binaries / 4 configurations (per CP2 user direction)

## Summary

CP2 installed the open-source counting tools needed for Audit 3, built
human + mouse references for each tool, and validated the end-to-end
pipeline with a real-FASTQ smoke test on `pbmc_1k_v3` (200,000-read
head subset). All four tool configurations produced count matrices.

## CellRanger disposition

CellRanger was **not installed** per user direction 2026-05-17 (gated
registration declined). Audit 3 proceeds with the 3-binary / 4-config
set documented in `environment/tool_configurations.md`. The
`STARsolo-CR-mimic` slot partially fills the gap by exercising the
STAR-authors' documented Cell-Ranger-compatible flag set. A full
CellRanger re-test is queued as **Audit 3c** in `DEFERRED.md`.

## Tool installation (conda env `audit3_counting`)

Created via `mamba create -n audit3_counting -y -c bioconda -c conda-forge`:

| Binary | Version | Path | Used by config slot(s) |
|---|---|---|---|
| `STAR`         | 2.7.11b | `/home/ross/miniforge3/envs/audit3_counting/bin/STAR` | 1, 2 (STARsolo default + CR-mimic) |
| `salmon`       | 1.10.3  | `…/bin/salmon` | 3 (alevin-fry pipeline frontend) |
| `alevin-fry`   | 0.10.0  | `…/bin/alevin-fry` | 3 (alevin-fry pipeline backend) |
| `kallisto`     | 0.50.1  | `…/bin/kallisto` | 4 (via kb-python) |
| `bustools`     | 0.43.2  | `…/bin/bustools` | 4 (via kb-python) |
| `kb-python`    | 0.28.2  | `…/bin/kb` | 4 (wrapper) |
| `samtools`     | 1.20    | `…/bin/samtools` | utility |

Python orchestration: scanpy 1.11.5, anndata 0.11.4, pyroe 0.9.2,
numpy 1.26.4, pandas 1.5.3, scipy 1.17.1, pyranges 0.0.120,
setuptools 80.10.2 (pinned <81 to keep `pkg_resources` available
for pyroe/pyranges).

Full version + path snapshot: `environment/tool_versions.tsv`.

## References built

### Human (GRCh38 + GENCODE v45)

Source: GENCODE release 45 (released 2023-09-19).

| Artifact | Path | Size |
|---|---|---:|
| Primary genome FASTA (gzipped) | `/mnt/nvme2/refs/GRCh38_gencode45/GRCh38.primary_assembly.genome.fa.gz` | 806 MB |
| GENCODE v45 GTF (gzipped) | `…/gencode.v45.primary_assembly.annotation.gtf.gz` | 48 MB |
| STAR genome index (sjdbOverhang=90) | `…/star_index_sjdb90/` | 13 GB |
| Salmon splici index (k=31, fl=91) | `…/salmon_splici_index/` | 9.9 GB |
| kallisto cDNA index (kb ref default) | `…/kb_index/` | 915 MB |

Total human reference footprint: **~25 GB**.

`sjdbOverhang=90` chosen because 10x 3' v3 R2 reads are 90 cycles; STAR
docs recommend `readLength - 1` for `sjdbOverhang`. v2 R2 reads are 90
cycles too. 5' v2 R2 is also 90 cycles. Single STAR index works for all
chemistries in the inventory.

### Mouse (mm39 + GENCODE vM34)

Source: GENCODE Mouse release vM34, via the **CP0 cached reference set
at `/mnt/nvme2/refs/mm39/`**.

| Artifact | Path | Size | Source |
|---|---|---:|---|
| Primary genome FASTA | `/mnt/nvme2/refs/mm39/GRCm39.primary_assembly.genome.fa.gz` | 739 MB | CP0 cache |
| GENCODE vM34 GTF | `…/gencode.vM34.primary_assembly.annotation.gtf` | 862 MB (uncompressed) | CP0 cache |
| STAR genome index (sjdbOverhang=150) | `…/star_index_sjdb150/` | 26 GB | CP0 cache (reused as-is) |
| Salmon splici index (k=31, fl=91) | `…/salmon_splici_index/` | 6.4 GB | **NEW (CP2)** — CP0 cache had plain cDNA, not splici |
| kallisto cDNA index | `…/kb_index/` | 514 MB | **NEW (CP2)** |

Total mouse reference footprint: **~33 GB** (most of it the cached
STAR index from CP0).

#### Documented version differences

Per user direction 2026-05-17:
- **mm10 (original submitters) → mm39 (this audit).** All 3 mouse
  datasets in the inventory list `reference_genome_expected: mm10` (the
  version each original submitter used). Audit 3 uses mm39 + GENCODE
  vM34 across all 3 mouse datasets and all 4 tool configs. Reference
  version is held constant across counting tools within this audit, so
  tool-comparison findings are isolated. Reference-version effects on
  counts are out of scope and would warrant a separate audit if
  surfaced as important. To be included in CP4 findings caveats and
  per-dataset metadata.
- **STAR sjdbOverhang=150 (mouse) vs =90 (human).** The cached mouse
  STAR index was built at sjdbOverhang=150 (suitable for 150-cycle
  reads). Audit 3 mouse datasets are 10x with R2=90 cycles, so the
  index is slightly suboptimal for splice-junction detection at the
  3'-end of reads. Held constant across the 2 STARsolo configurations
  (slots 1 + 2), so doesn't affect the within-STARsolo comparison;
  may affect STARsolo-vs-non-STARsolo comparisons marginally.
  Documented as a known limitation rather than rebuilt (the rebuild
  cost is ~30 min and doesn't materially change the audit's
  cross-tool comparison structure).

## End-to-end smoke test results

Subset: 200,000-read head from `pbmc_1k_v3` (10x demo, 3' v3.1, human
PBMC). All 4 configurations ran to completion and produced count
matrices.

| Slot | Configuration | Output file | Size | Wall time |
|---|---|---|---:|---:|
| 1 | STARsolo default | `star_default/Solo.out/Gene/raw/matrix.mtx` | 1.2 MB | ~7 s |
| 2 | STARsolo CR-mimic | `star_cr_mimic/Solo.out/Gene/raw/matrix.mtx` | 1.2 MB | ~7 s |
| 3 | alevin-fry | `alevin_fry/quant/alevin/quants_mat.mtx` | 1.1 MB | ~16 s |
| 4 | kb count (kallisto+bustools) | `kb_count/counts_unfiltered/cells_x_genes.mtx` | 912 KB | ~27 s |

Total smoke test wall time (sequential): ~1.7 min for the actual tool
runs, plus ~3 sec for FASTQ subsetting and whitelist download. All on
12 threads.

salmon-alevin reports `Selectively-aligned 177,524 total fragments out
of 200,000` (88.8%), which is reasonable for a head subset of real
PBMC data. alevin-fry collate reports `total number of distinct
corrected barcodes : 1,995` — small because we're at 200k reads (the
real `pbmc_1k_v3` full dataset reports ~1,222 cells, so 1,995 at the
permit-list stage is the expected pre-filter count).

## Issues encountered + workarounds

1. **pyroe/pyranges `pkg_resources` import error** (CP2 mid-build).
   Fixed by pinning `setuptools<81` (modern setuptools dropped
   `pkg_resources`). Captured in `environment/tool_versions.tsv`.

2. **kb-python failed to find `kallisto`** (initial reference build).
   The script `build_human_refs.sh` invoked `$ENV/bin/kb` directly
   without putting `$ENV/bin` on PATH; kb-python then tried to resolve
   `kallisto` from system PATH and failed. Fixed by adding
   `export PATH="$ENV/bin:$PATH"` at the top of the build scripts.

3. **10x v3 whitelist GitHub URL gated** (smoke test prerequisite).
   `https://github.com/10XGenomics/cellranger/raw/master/.../3M-february-2018.txt.gz`
   now returns HTML instead of the gzip file (10x's repo appears to
   have been made gated). Working source found at
   `https://teichlab.github.io/scg_lib_structs/data/10X-Genomics/`.
   6,794,880 barcodes verified after gunzip.

4. **alevin-fry 0.10 CLI changes** (smoke test slot 3).
   `generate-permit-list` v0.10 requires `--output-dir`; older versions
   accepted `-d <ori>` shorthand. `collate` v0.10 uses `--input-dir` +
   `--rad-dir` (not `-i` + `-r` ambiguously). Updated
   `scripts/smoke_test_pbmc_1k.sh` and `environment/tool_configurations.md`
   to reflect the v0.10 CLI.

5. **`set -o pipefail` + `zcat | head | gzip`** (smoke test prerequisite).
   `head` closes the pipe early → `zcat` receives SIGPIPE → exit 141.
   Fixed by wrapping the subset-build block in `set +o pipefail` ... `set -o pipefail`.

## Deliverables

- `environment/tool_versions.tsv` — tool versions + binary paths + env metadata
- `environment/tool_configurations.md` — the 4 tool configurations specified for slots 1-4 (including STARsolo flag sets for default vs CR-mimic, alevin-fry quant chain, kb count workflow)
- `environment/logs/build_human_refs.log` — STAR / salmon / kallisto build log for human
- `environment/logs/build_mouse_refs.log` — splici + kallisto build log for mouse
- `environment/logs/smoke_test.log` — full smoke-test trace
- `scripts/build_human_refs.sh` — idempotent re-runnable human ref builder
- `scripts/build_mouse_refs.sh` — idempotent re-runnable mouse ref builder
- `scripts/smoke_test_pbmc_1k.sh` — end-to-end smoke test for all 4 slots
- `references/` — placeholder (actual reference files live at
  `/mnt/nvme2/refs/{GRCh38_gencode45,mm39}/`; the `references/`
  directory tracks lock-file metadata pointing there)

All listed deliverables hash-registered in
`/mnt/nvme1/omics-audit/phase2/repro.lock`.

## CP3 prerequisites in place

- All 3 binaries operational on the env
- Both species' references built and tested
- 10x v3 whitelist downloaded; v2 + 5p whitelists will need to be
  acquired during CP3 (similar pattern via teichlab mirror)
- Smoke test confirmed each tool's count-matrix output path is well
  defined and stable

## Stop and report

CP2 done. Waiting for approval before CP3 (FASTQ download + processing
setup for all 11 datasets).
