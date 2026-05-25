#!/usr/bin/env python
"""Idempotently apply the CellBender checkpoint weakref-strip patch.

WHY
    CellBender 0.3.2 crashes at checkpoint save with
        TypeError: cannot pickle 'weakref.ReferenceType' object
    because pyro caches `param.unconstrained = weakref.ref(...)` on every
    constrained nn.Parameter, and torch.save(model_obj) tries to pickle it.
    No output is produced. This inserts a strip of those cached weakrefs just
    before the save. See cellbender_env.md for the full diagnosis.

USAGE (run with the ambient_cb interpreter; works from any directory / a clean clone)
    ~/miniforge3/envs/ambient_cb/bin/python apply_checkpoint_patch.py          # apply
    ~/miniforge3/envs/ambient_cb/bin/python apply_checkpoint_patch.py --check  # report only, no write
    ~/miniforge3/envs/ambient_cb/bin/python apply_checkpoint_patch.py --force  # patch even on a non-0.3.2 version

EXIT CODES
    0  patch applied now, OR already present (idempotent)
    1  wrong interpreter / CellBender not importable
    2  CellBender version mismatch (use --force to override after review)
    3  anchor line not found — CellBender layout changed; patch manually
    4  post-write verification failed (marker missing or file won't compile)
"""
import argparse
import py_compile
import sys
import tempfile

EXPECTED_VERSION = "0.3.2"
MARKER = "PATCH (audit: ambient correction"
ANCHOR = "            torch.save(model_obj, filebase + '_model.torch')"

PATCH = '''            # --- PATCH (audit: ambient correction CP0, 2026-05-24) ---
            # pyro caches `param.unconstrained = weakref.ref(...)` on every
            # constrained nn.Parameter; torch.save can't pickle the weakref.
            # Strip before saving; pyro recreates lazily, posterior reloads fresh.
            def _strip_unconstrained_weakrefs(m):
                if m is None:
                    return
                for _p in getattr(m, "_parameters", {}).values():
                    if _p is not None and hasattr(_p, "unconstrained"):
                        try:
                            object.__delattr__(_p, "unconstrained")
                        except Exception:
                            try:
                                del _p.unconstrained
                            except Exception:
                                pass
                for _c in getattr(m, "_modules", {}).values():
                    _strip_unconstrained_weakrefs(_c)
            _enc = getattr(model_obj, "encoder", None)
            if isinstance(_enc, dict):
                for _v in _enc.values():
                    _strip_unconstrained_weakrefs(_v)
            else:
                _strip_unconstrained_weakrefs(_enc)
            _strip_unconstrained_weakrefs(getattr(model_obj, "decoder", None))
            # --- END PATCH ---
'''


def fail(code, msg):
    print(f"FAIL [{code}]: {msg}", file=sys.stderr)
    sys.exit(code)


def main():
    ap = argparse.ArgumentParser(description="Apply CellBender checkpoint weakref-strip patch.")
    ap.add_argument("--check", action="store_true", help="report status only; do not modify")
    ap.add_argument("--force", action="store_true", help="apply even if CellBender version != %s" % EXPECTED_VERSION)
    args = ap.parse_args()

    # 1. CellBender importable with the active interpreter?
    try:
        import cellbender
        import cellbender.remove_background.checkpoint as ckpt
    except Exception as e:
        fail(1, f"cannot import cellbender ({e}). Run with the ambient_cb interpreter: "
                f"~/miniforge3/envs/ambient_cb/bin/python {sys.argv[0]}")

    version = getattr(cellbender, "__version__", "unknown")
    path = ckpt.__file__
    print(f"CellBender {version}  ->  {path}")
    print(f"interpreter: {sys.executable}")

    # 2. version guard
    if version != EXPECTED_VERSION and not args.force:
        fail(2, f"version {version} != developed-against {EXPECTED_VERSION}. "
                f"Review checkpoint.py, then re-run with --force if the anchor/save logic is unchanged.")

    src = open(path, encoding="utf-8").read()

    # 3. idempotent: already patched?
    if MARKER in src:
        print("STATUS: already patched (no change).")
        sys.exit(0)

    # 4. anchor present?
    if ANCHOR not in src:
        fail(3, f"anchor line not found in {path}; CellBender checkpoint layout changed. "
                f"Patch manually: strip `param.unconstrained` weakrefs before torch.save(model_obj).")

    if args.check:
        print("STATUS: NOT patched; anchor present; patchable. (--check: no changes written.)")
        sys.exit(0)

    # 5. apply (single, anchored insertion)
    if src.count(ANCHOR) != 1:
        fail(3, f"expected exactly 1 anchor, found {src.count(ANCHOR)}; refusing ambiguous patch.")
    new_src = src.replace(ANCHOR, PATCH + ANCHOR, 1)

    # 6. verify the patched source compiles BEFORE overwriting the installed file
    with tempfile.NamedTemporaryFile("w", suffix=".py", delete=False, encoding="utf-8") as tf:
        tf.write(new_src)
        tmp_path = tf.name
    try:
        py_compile.compile(tmp_path, doraise=True)
    except py_compile.PyCompileError as e:
        fail(4, f"patched source failed to compile, NOT written: {e}")

    open(path, "w", encoding="utf-8").write(new_src)

    # 7. post-write verification
    check = open(path, encoding="utf-8").read()
    if MARKER not in check:
        fail(4, "post-write check: marker missing after write.")
    try:
        py_compile.compile(path, doraise=True)
    except py_compile.PyCompileError as e:
        fail(4, f"post-write check: installed file won't compile: {e}")

    print(f"STATUS: patched OK ({path}).")
    sys.exit(0)


if __name__ == "__main__":
    main()
