# Audit-repo integration notes

Cross-cutting concerns that surfaced during audit work but aren't
in-scope for any specific audit. Triage at audit closeout or before
the next major audit starts.

---

## Known issues

- **`phase2/audit1_main/databases/msigdb_c5_go_bp_human.tsv` is 57 MB in git history (committed 2026-05-17 with schema-v1.0.3 merge).** Annoying but not blocking — GitHub warns but allows up to 100 MB per file. Fix options when it matters: (a) `git filter-repo` to rewrite history (destructive, requires force-push to origin/main); (b) migrate to Git LFS (adds friction + cost). Neither worth doing reactively. Re-evaluate if the audit corpus accumulates more large cache files and the repo size becomes a real concern.
