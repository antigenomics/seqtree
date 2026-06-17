Roadmap
=======

The v1 core ships both engines, scope + budget search, BLOSUM62 / custom matrices with linear
gaps, top / all hits, on-demand alignment, parallel batch and batch-of-batches, and auto-selecting
pairwise search. Planned next:

1. **Alignment polish** — affine gaps (gap open vs extend), CIGAR-style output, batch alignment.
2. **Position-specific scoring matrices** — ``penalty(pos, a, b)``; consume per-position PWMs.
3. **E-value / significance.** TCR and epitope reference sets are highly redundant and clustered,
   so a naive BLAST-style E-value over-counts. Two directions under consideration:

   - *Control-set expected hits* — the caller supplies a background set matched to the relevant
     null (pgen-matched TCRs; anchor-masked random peptides over an MHC distribution; read
     fractions for nucleotide errors). Searching it yields expected hits at each score → an
     empirical E-value. Payload-agnostic; a thin layer over ``search_batch``.
   - *tf-idf redundancy weighting* — down-weight hits to common / clustered motifs intrinsically,
     accounting for redundancy without an external control set.

   The direction will be chosen once the core ships and null-score distributions can be measured
   on real VDJdb / IEDB data.

4. **Memory** — succinct (LOUDS) trie / packed reference strings if reference counts reach the
   10M+ range; memory-mapped frozen index for zero-copy load.
5. **Batch-of-batches & distribution** — streaming query batches and optional process-level fan-out.

Downstream consumers (built in *other* libraries, not here): UMI collapsing, CDR3 nucleotide error
correction, VDJdb CDR3 matching, IEDB epitope matching. seqtree stays payload-agnostic; those
libraries apply V-gene / MHC / count filters on top of ``ref_id`` results.
