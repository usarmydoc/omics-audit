# Audit 3 — Tool configurations

**Date:** 2026-05-17
**Env:** `audit3_counting` (conda) at `/home/ross/miniforge3/envs/audit3_counting`

Audit 3 compares **4 distinct counting configurations** built from **3 binaries**:

| Slot | Binary | Configuration | Algorithm family |
|---|---|---|---|
| 1 | `STAR 2.7.11b` | STARsolo with vanilla defaults | Splice-aware full alignment + STARsolo default UMI dedup |
| 2 | `STAR 2.7.11b` | STARsolo with CellRanger-compatible flag set | Splice-aware full alignment + CR-style UMI dedup (replicates Cell Ranger 4-7 output) |
| 3 | `salmon 1.10.3 + alevin-fry 0.10.0` | Default alevin-fry pipeline (splici reference, knee-point cell filtering) | Selective alignment with k-mer hashing |
| 4 | `kallisto 0.50.1 + bustools 0.43.2 + kb-python 0.28.2` | kb count workflow with default options | Pseudoalignment with equivalence classes |

CellRanger is intentionally absent (user direction 2026-05-17 — registration
not pursued). The STARsolo-CR-mimic slot fills part of the gap by exercising
the official CellRanger-compatible flag set documented by the STAR authors.
A full CellRanger re-test is queued in `DEFERRED.md` (Audit 3c).

## Slot 1: STARsolo (default)

```
STAR \
  --runMode alignReads \
  --soloType CB_UMI_Simple \
  --soloCBwhitelist <10x_whitelist_per_chemistry> \
  --soloFeatures Gene \
  --genomeDir <STAR_index> \
  --readFilesIn <R2.fastq.gz> <R1.fastq.gz> \
  --readFilesCommand zcat \
  --soloCBstart 1 --soloCBlen <BC_LEN> \
  --soloUMIstart <UMI_START> --soloUMIlen <UMI_LEN> \
  --outSAMtype BAM SortedByCoordinate \
  --runThreadN <N>
```

Default UMI dedup and cell-calling: `--soloUMIdedup 1MM_All` (default),
`--soloCBmatchWLtype 1MM_All` (default). No CR-specific flags applied.

## Slot 2: STARsolo (CR-mimic)

Per STAR authors' recommendations for matching Cell Ranger output
(documented in STARsolo manual §3 "Matching CellRanger results"):

```
STAR \
  ... (same as Slot 1 base) ... \
  --soloUMIfiltering MultiGeneUMI_CR \
  --soloUMIdedup 1MM_CR \
  --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
  --clipAdapterType CellRanger4 \
  --outFilterScoreMin 30 \
  --soloFeatures Gene
```

The differences from Slot 1 collectively replicate Cell Ranger 4-7's
UMI deduplication, cell-barcode matching, and adapter clipping
behavior. Per benchmarks in Kaminow et al. 2021 (STARsolo paper),
this configuration matches Cell Ranger output to within ~1-2% on the
PBMC 5k v3 dataset.

## Slot 3: alevin-fry

Reference: splici (spliced + intronic) transcriptome built with pyroe.
Quant chain: `salmon alevin` → `alevin-fry generate-permit-list` →
`alevin-fry collate` → `alevin-fry quant`.

```
salmon alevin \
  -i <salmon_splici_index> \
  -l ISR \
  --chromiumV3 (or --chromium for v2) \
  -1 <R1.fastq.gz> -2 <R2.fastq.gz> \
  -p <N> --sketch \
  -o <output>

alevin-fry generate-permit-list \
  -i <output> -d fw \
  --expected-ori fw \
  --unfiltered-pl

alevin-fry collate -t <N> \
  -i <output> -r <output>

alevin-fry quant -t <N> \
  -i <output> \
  -m <t2g_3col.tsv> \
  -r cr-like \
  --use-mtx \
  -o <output>/quant
```

Resolution mode: `cr-like` mimics Cell Ranger's USA-mode quantification.
Permit list uses `--unfiltered-pl` (knee-point cell filtering at quant
stage; the alevin-fry equivalent of CR's cell calling).

## Slot 4: kallisto-bustools (kb count)

Reference: kallisto cDNA index + t2g mapping built with `kb ref`.

```
kb count \
  -i <kallisto_index> -g <t2g.txt> \
  -x 10xv3 (or 10xv2, 5p_v2 per chemistry) \
  -o <output> \
  --filter bustools \
  -t <N> \
  <R1.fastq.gz> <R2.fastq.gz>
```

`--filter bustools` invokes the knee-point filter from bustools.
Resolution: `--workflow standard` (default), the equivalent of
Cell Ranger's exon counting.

## Reference compatibility across slots

All 4 slots build from the **same** primary-genome FASTA + GTF pair
(GENCODE human v45 / GENCODE mouse vM34), so reference content is held
constant across tools. Within each tool, the per-tool index is the
natural artifact (STAR genome index, salmon splici index, kallisto
cDNA index) — different formats but encoding the same underlying
reference.

`reference_genome_expected: mm10` in the inventory refers to what the
**submitters** used. Audit 3 uses mm39 + GENCODE vM34 for mouse data
and documents the version difference per user direction (parsimonious
choice; ref-version effect on counts is out of scope and would warrant
a separate audit if surfaced).

## Multi-threading and resource budget

- STARsolo: ~30 GB RAM peak for GRCh38; ~16 GB for mm39. Use 12-16 threads.
- alevin-fry (`salmon alevin`): ~10 GB RAM. Use 12-16 threads.
- kallisto: ~5 GB RAM. Use 12-16 threads.
- bustools: ~2 GB RAM. Use 8-12 threads.

Concurrent dataset processing budget on 64 GB RAM machine: 1 STARsolo +
1 alevin-fry + 1 kallisto pipeline simultaneously stays under 50 GB.
