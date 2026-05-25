# CellBender environment — `ambient_cb`

_CP0, 2026-05-24. Audit: ambient RNA correction._

## Conda env
- Name: `ambient_cb` (python 3.10), created with `mamba`, isolated from `liver_scvi`.
- Key versions (CellBender-only env; downstream scanpy work runs in the `audit3_counting` env):

| pkg | version |
|-----|---------|
| cellbender | 0.3.2 |
| torch | 2.3.1+cu121 |
| pyro-ppl | 1.8.6 |
| numpy | 1.26.4 |
| anndata | 0.11.4 |
| scipy | 1.15.3 |

GPU: RTX 4070 Ti (cc 8.9), CUDA 12.1 runtime via torch wheel. `torch.cuda.is_available() == True`.

## Required source patch — `checkpoint.py` (the actual fix)

**Symptom:** `cellbender remove-background` trains fine on GPU, then crashes at
checkpoint save with `TypeError: cannot pickle 'weakref.ReferenceType' object`
(`checkpoint.py:115`, `torch.save(model_obj, ...)`), which cascades to the
posterior step failing (`Checkpoint file ckpt.tar.gz does not exist`). No output
is produced.

**Root cause (diagnosed by recursive pickle probe of `model_obj`):** pyro caches
`param.unconstrained = weakref.ref(...)` on **every constrained `nn.Parameter`**
in the encoder/decoder (e.g. `encoder.z._modules.loc_out._parameters.weight.unconstrained`).
`torch.save` pickles the full model object including these weakrefs and fails.

**Reproduced across all version axes — NOT a version-pin issue:**
torch 2.11 & 2.3.1; pyro 1.9.1 & 1.8.6; numpy 2.2.6 & 1.26.4; cellbender 0.3.0 & 0.3.2.
All fail identically. The version pins above are a conservative known-good
CellBender-era stack left in place after debugging; **the patch below is what
actually resolves it**, independent of versions.

**Fix:** strip the cached `.unconstrained` weakrefs off all encoder/decoder
parameters immediately before `torch.save(model_obj, ...)`. Pyro lazily recreates
them on next access, and the posterior step reloads the checkpoint fresh — so
stripping is safe. Verified: full end-to-end run on `10x_pbmc_1k_v3` produces
`_filtered.h5` (1288 cells × 63241 genes), `_posterior.h5`, `_metrics.csv`, `.pdf`,
loadable via `cellbender.remove_background.downstream.anndata_from_h5`.

Patch location: inserted in `save_checkpoint()` just before the first
`torch.save(model_obj, ...)`. Reapply idempotently with
`apply_checkpoint_patch.py` in this directory (run after any env rebuild or
cellbender reinstall).

## Known non-fatal issues (do NOT block the audit)
- **HTML report fails** (`nbconvert`: `tmp.report.nbconvert.html` not found) *after*
  "Completed remove-background." The PDF report and all matrices are written first.
  Cosmetic only — the audit consumes corrected matrices + `_metrics.csv`, not the HTML.
- **Generic `anndata.read_h5ad` rejects CellBender output** (`unexpected keyword
  argument 'droplet_latents'`). Use CellBender's `anndata_from_h5` reader instead.

## Input staging (per dataset)
STARsolo writes plain v3-named files (`features.tsv`); CellBender misreads unzipped
files as CellRanger v2 and demands `genes.tsv`. Stage a gzipped v3 layout first:
```
gzip -c Gene/raw/matrix.mtx   > raw_v3/matrix.mtx.gz
gzip -c Gene/raw/features.tsv > raw_v3/features.tsv.gz
gzip -c Gene/raw/barcodes.tsv > raw_v3/barcodes.tsv.gz
```
Then `--input raw_v3`.
