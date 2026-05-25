# Audit — Ambient RNA Correction — CP0 Inventory & Feasibility

_Generated: 2026-05-24_
_Standards: AUDIT_STANDARDS.md v1.0.3 + §5.3.2 equivalence amendment_
_Status: CP0 COMPLETE — all tools verified, CellBender GPU smoke test PASSED (end-to-end), scope decisions confirmed. Ready for CP1 on approval._

## Scope recap (locked)

Three-tool, three-deliverable audit of scRNA-seq ambient RNA correction.

| Tool | Lang | Runs on | Output |
|------|------|---------|--------|
| SoupX | R | filtered cells + clustering (soup from empty droplets) | corrected counts + per-cell ρ (contamination fraction) |
| CellBender | Python (GPU) | raw barcode×gene matrix | corrected matrix + per-cell/per-gene posteriors; does its own cell-calling |
| DecontX (celda) | R | filtered cells (+ optional raw background) | corrected counts + per-cell contamination estimate |

Deliverables: **A** per-gene/per-cell contamination estimate comparison across 3 tools × 9 datasets (mirror Audit 3 C1); **B** ordering analysis (correction vs QC vs cell-calling); **C** biological propagation on intestine + PBMC control (mirror Audit 3 C3/CP6).

---

## Step 1 — Working set status

**Cross-audit working set = the same 9 Audit 3 datasets. ALL PRESENT on disk. No re-counting needed.**

The recall that `processed/` was cleaned during Audit 3 wind-down is **incorrect** — `phase2/audit3_counting/processed/` is fully intact, with per-tool count matrices for all 9 datasets. Every dataset has a complete **raw (unfiltered) + filtered** matrix pair from `kb_count`, and STARsolo `Gene/{raw,filtered}` in CellRanger format.

| Dataset | Tissue/type | Raw barcodes | Filtered cells | Genes | Raw nnz |
|---------|-------------|-------------:|---------------:|------:|--------:|
| 10x_pbmc_1k_v3 | PBMC (v3) | 239,943 | 1,194 | 63,241 | 3.5M |
| 10x_neuron_1k_v3 | neuron (v3) | 381,166 | 1,321 | 57,180 | 6.4M |
| 10x_t_3k_v2 | T cell (v2) | 357,347 | 3,879 | 63,241 | 7.6M |
| 10x_pbmc_4k_v2 | PBMC (v2) | 285,756 | 4,680 | 63,241 | 8.8M |
| 10x_pbmc_5k_v3.1 | PBMC (v3.1) **[C-control]** | 533,414 | 4,911 | 63,241 | 11.1M |
| 10x_pbmc_10k_v3.1 | PBMC (v3.1) | 1,439,561 | 11,398 | 63,241 | 27.5M |
| gse325955_mouse_kidney_E18_5 | mouse kidney | 2,736,613 | 25,171 | 57,180 | 101M |
| gse287209_human_lung_organoid | lung organoid | 3,159,866 | 33,340 | 63,241 | 128M |
| gse288156_mouse_intestine_scrna | mouse intestine **[C-stress, high ambient]** | 4,073,122 | 67,132 | 57,180 | 237M |

- **CellBender** (needs raw): ✓ all 9 — `kb_count/counts_unfiltered/cells_x_genes.mtx` and `star_default/Solo.out/Gene/raw/{matrix.mtx,barcodes.tsv,features.tsv}`.
- **SoupX / DecontX** (need filtered cells): ✓ all 9 — `kb_count/counts_filtered/` and `star_default/Solo.out/Gene/filtered/`.
- FASTQs: not required (no re-generation). Per Audit 3 inventory they remain available if ever needed, but this audit operates entirely on existing matrices.

**No scope decision needed on re-counting** (unlike QC-MAD CP0) — matrices are complete.

### Scope decision A (surface for confirmation): fix ONE upstream counting tool

To avoid confounding the ambient-correction comparison with Audit 3's counting-tool variation, the audit should hold the **upstream counting tool constant** and apply all three correction methods on top of it. Recommendation: **STARsolo (`star_default`)** — it is (a) the Audit 3 C3 reference/most-conservative caller, (b) native CellRanger-format (`matrix.mtx`+`barcodes.tsv`+`features.tsv`, genes×cells) which is the **direct native input for all three tools**, removing a conversion step. `kb_count` matrices remain available as a sensitivity check if desired. _Awaiting confirmation._

---

## Step 2 — Tool installation status

| Tool | Status | Version | Location |
|------|--------|---------|----------|
| SoupX (R) | ✓ installed | 1.6.2 | user lib `~/R/x86_64-pc-linux-gnu-library/4.5` |
| DecontX via celda (R) | ✓ installed | celda 1.26.0 (`decontX()` exported & loads) | user lib |
| DropletUtils (R, support) | ✓ present | 1.30.0 | for `read10xCounts`/empty-droplet infra |
| Seurat / SingleCellExperiment | ✓ present | 5.5.0 / 1.32.0 | |
| CellBender (Python) | ✓ installed + **patched** | 0.3.2 (torch 2.3.1+cu121, pyro 1.8.6, numpy 1.26.4) | dedicated conda env `ambient_cb` |

Notes:
- **R packages installed by the user** (per workflow rule — Claude does not run R installs). SoupX came from CRAN; celda from Bioconductor 3.22 with `update=FALSE` (the `spatial`-in-`/usr/lib/R/library` unwritable warning is benign and was skipped).
- **DecontX is NOT a separate package** — `decontX()` is exported by `celda` 1.26.0. (A standalone `decontX` Bioc package also exists but is unnecessary.)
- **CellBender env is isolated** as `ambient_cb` (python 3.10) — deliberately NOT installed into `liver_scvi` to avoid disturbing that env's scVI/torch stack.
- **CellBender required a source patch** to produce any output (see Step 4 + `environment/cellbender_env.md`). The version pins are a conservative known-good stack left in place after debugging; the patch — not the versions — is the fix.

---

## Step 3 — Input-format compatibility per tool

Using STARsolo `Gene/{raw,filtered}` (CellRanger format) as the fixed upstream (Scope decision A):

| Tool | Native input | From STARsolo output | Conversion needed |
|------|--------------|----------------------|-------------------|
| SoupX | `load10X()` on CellRanger dir (raw + filtered), needs clustering | direct (raw + filtered dirs) | none; clustering computed in-pipeline |
| CellBender | CellRanger h5 / mtx dir / h5ad; raw matrix | direct (`Gene/raw` mtx dir) | none |
| DecontX | SingleCellExperiment of filtered counts (+ optional raw background) | `read10xCounts(Gene/filtered)` → SCE | trivial (DropletUtils) |

If `kb_count` is used instead, its `cells_x_genes.mtx` is **cells×genes (transposed)** vs CellRanger genes×cells — would require a transpose + h5ad/SCE build step. STARsolo avoids this. No information loss either way.

---

## Step 4 — GPU verification for CellBender

- `nvidia-smi`: RTX 4070 Ti, 12,282 MiB total, **idle** (214 MiB used, 3% util) — full VRAM available.
- PyTorch (verified in `liver_scvi`, torch 2.10+cu128): `cuda.is_available() == True`, device = "NVIDIA GeForce RTX 4070 Ti", CUDA 12.8. Compute capability 8.9 (Ada) supported.
- **CellBender GPU smoke test on `10x_pbmc_1k_v3` (STARsolo raw): PASSED end-to-end.** Trained 20 epochs on the 4070 Ti (~45 s, CUDA confirmed), checkpoint saved, posterior computed (11 chunks), corrected outputs written: `_filtered.h5` (1288 cells × 63241 genes), `.h5`, `_posterior.h5`, `_metrics.csv`, `.pdf`. Metrics sane (9.8% ambient removed; 1288 cells vs 1200 expected). Output loads via `cellbender.remove_background.downstream.anndata_from_h5`. Evidence in `cp0/smoke_test/` (`SMOKE_PASS.log`, `ROOT_CAUSE_probe.log`, metrics, PDF).

**Resolved blocker — CellBender 0.3.2 needed a source patch to produce any output.** Out of the box it trains on GPU but crashes at checkpoint save: `TypeError: cannot pickle 'weakref.ReferenceType'`. Diagnosed (recursive pickle probe) to pyro caching `param.unconstrained = weakref.ref(...)` on every constrained `nn.Parameter`; `torch.save(model_obj)` can't pickle it. Reproduced across torch 2.11/2.3.1, pyro 1.9.1/1.8.6, numpy 2.2.6/1.26.4, cellbender 0.3.0/0.3.2 — **not a version issue**. Fix: strip those cached weakrefs before save (pyro recreates lazily; posterior reloads fresh). Patch + idempotent reapply script + full diagnosis in `environment/cellbender_env.md` and `environment/apply_checkpoint_patch.py`.

Two **non-fatal** CellBender quirks (do not block — the audit consumes corrected matrices + metrics):
- HTML report step fails (`nbconvert`) *after* completion; the PDF report and all matrices are already written.
- Generic `anndata.read_h5ad` rejects CellBender output (`droplet_latents` kwarg); use CellBender's `anndata_from_h5`.

Two input-staging findings:
- STARsolo writes plain v3-named files → CellBender misreads as CellRanger v2 and wants `genes.tsv`. Stage a gzipped v3 layout (`matrix.mtx.gz`/`features.tsv.gz`/`barcodes.tsv.gz`) first. Verified.

---

## Step 5 — Ordering analysis (sensible-comparison matrix)

The task's 4 orderings assume cell-calling, correction, and QC are freely reorderable. They are **not**, because two of three tools structurally require cells to be called before they can run:

| Tool | Cell-calling | Reason |
|------|-------------|--------|
| SoupX | **must precede** correction | estimates the soup from empty droplets but needs filtered cells + clustering to assign per-cell ρ |
| DecontX | **must precede** correction | operates on an already-filtered SCE of called cells (+ optional raw background) |
| CellBender | **performed by the tool itself** | jointly infers ambient + which barcodes are cells from the raw matrix |

Consequence — the "correction before vs after cell-calling" axis is **only meaningful for CellBender** (the sole tool that runs pre-cell-calling, which is its native mode). For SoupX/DecontX, cell-calling is fixed upstream. So orderings 1 and 3 (correction *before* cell-calling) are **impossible for SoupX and DecontX** and excluded.

**The axis that IS testable across all three tools is correction-vs-QC:**

| Condition | SoupX | CellBender | DecontX |
|-----------|:-----:|:----------:|:-------:|
| **O1: correct → QC** (correct called cells, then QC-filter) | ✓ | ✓ (raw→CB(corr+CC)→QC) | ✓ |
| **O2: QC → correct** (QC-filter called cells, then correct) | ✓ | ✓ (QC raw barcodes→CB) | ✓ |

→ Minimum sensible set = **2 orderings × 3 tools** (CellBender's native raw mode is folded into O1). The "before/after cell-calling" framing is documented as collapsing to this. _Recommend Deliverable B run O1+O2 for all 3 tools; impossible orderings explicitly excluded with rationale above._

---

## Step 6 — Deliverable C reusability

Audit 3 C3/CP6 pipeline is **fully scripted and reusable**:
- `audit3_counting/scripts/cp6_pipeline.py` — scanpy: shifted-log norm → scry deviance HVG (top-2000) → PCA(50) → neighbors(15) → Leiden(res 1.0, igraph) → Wilcoxon markers → CellTypist majority-vote.
- `cp6_deviance_hvg.R` (scry), `cp6_scdblfinder.R` (doublets), `cp6_compare.py` (ARI/NMI/marker overlap), `cp6_run_all.sh`.
- Parameters held constant in `c3/pipeline_parameters.yaml` (documented, no silent defaults).
- Datasets: **control = 10x_pbmc_5k_v3.1**, **stress = gse288156 mouse intestine** (high ambient) — exactly the task's intestine + PBMC control.
- CellTypist models on hand: `Adult_Mouse_Gut.pkl` (intestine), `Immune_All_Low.pkl` (PBMC).
- **Existing uncorrected C3 baseline** (`c3/intestine/`, `per_tool_pipeline_outputs/`) serves as the no-correction reference — the audit adds correction arms on top.

Plan: for the 2 datasets, run {3 tools × O1/O2} corrected inputs through the cp6 pipeline, compare cluster/annotation/marker structure vs the uncorrected baseline. Tests whether ambient correction reduces Audit 3 C3's high-ambient permissiveness-chain divergence.

---

## Step 7 — Heumos positioning

- **Deliverable A** — *extends* Heumos: Heumos recommends ambient correction but does not quantify how much the three tools' contamination estimates differ on identical input.
- **Deliverable B** — *fills gap*: Heumos does not address ordering (correct vs QC; pre/post cell-calling).
- **Deliverable C** — *extends Audit 3 C3 + Heumos*: tests the practical/biological consequence of correction (does it close C3's high-ambient gap?).

Finding type is open (equivalence "tools agree" vs difference "they diverge") → **§5.3.2 equivalence-tier criteria apply** for tier assignment when results land, alongside §5.3.1 if a method dominates.

---

## Run-count & feasibility estimate (hardware: 9900X 12c/24t, 60 GB, RTX 4070 Ti 12 GB)

| Deliverable | Runs | Notes |
|-------------|------|-------|
| A (contamination estimates) | 3 tools × 9 datasets = **27** | SoupX/DecontX minutes each (R); CellBender GPU 5–20 min/dataset |
| B (ordering) | 2 orderings × 3 tools × 9 = **54** | shares correction runs with A where overlapping |
| C (biological propagation) | 3 tools × 2 orderings × 2 datasets = **12** downstream scanpy runs + existing uncorrected baseline | intestine 67k cells is the heaviest |
| **Total** | **~90–95 correction/pipeline runs** | |

CellBender on 4070 Ti will beat published V100/A100-era benchmarks; intestine (4M droplets, `--total-droplets-included` bounded) is the long pole at perhaps 20–40 min. SoupX/DecontX single-threaded R at high Zen5 clocks = minutes. **Total compute ≈ 2–3 days of background runs** — consistent with the stated envelope. No step is predicted to materially exceed this.

---

## Blocking issues

- **None.** The one risk (CellBender GPU output) is resolved: smoke test passes end-to-end with the documented patch. SoupX/DecontX are standard R and not yet smoke-tested, but their inputs (STARsolo filtered cells) and packages are confirmed present — they will be exercised in CP1.

## Scope decisions — CONFIRMED (user, 2026-05-24)

1. **Upstream counting tool fixed = STARsolo `star_default`** ✓ (kb_count available as optional sensitivity check).
2. **Ordering design = O1 (correct→QC) + O2 (QC→correct) across all 3 tools** ✓; pre-cell-calling orderings excluded for SoupX/DecontX (structurally impossible), folded into CellBender's native raw mode.
3. **Deliverable C reuses the existing Audit 3 C3 uncorrected baseline** ✓ as the no-correction reference.

## Recommendation

**CP0 complete — ready for CP1 (Deliverable A) on user go-ahead.** All matrices present (no re-counting/substitution/acquisition), all three tools installed and verified (CellBender end-to-end on GPU with documented patch), GPU free, downstream C3/CP6 pipeline reusable, scope decisions locked.
