# Reference literature

Reference literature for positioning audit findings against published consensus.

Each file here is a published best-practices source. The accompanying
`_index.md` file maps the source's specific recommendations to existing audit
coverage, allowing audit findings to be tagged against the consensus they
extend, confirm, or contradict.

This is NOT a reading list of papers to audit. These are reference documents
used to sharpen the framing of existing audit findings.

Temporal context matters. Bioinformatics best-practices papers age fast in
some subdomains (spatial transcriptomics, deep learning methods, multi-omics
integration) and slow in others (architectural recommendations like pseudobulk
DE, Leiden clustering, MAD-based QC). Each index notes temporal context per
recommendation.

When citing a reference in a rule's prior_audit_relationship or in findings.md:
- "extends Heumos 2023 recommendation by quantifying ..."
- "confirms Heumos 2023 recommendation across N datasets ..."
- "contradicts Heumos 2023 implicit assumption that ..."
- "addresses gap in Heumos 2023 where field had moved on since publication"

Do not add new reference documents here without an immediate use case in an
in-progress or completed audit. The reference library exists to sharpen
existing audit framing, not to accumulate papers for potential future use.

Current references:
- heumos_2023: Heumos et al. 2023, "Best practices for single-cell analysis
  across modalities" (Nat Rev Genet). Comprehensive scRNA-seq + adjacent
  modality review. 3 years old at time of indexing; fast-moving subdomains
  have moved on.
