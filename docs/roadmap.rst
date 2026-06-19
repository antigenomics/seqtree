Roadmap
=======

The core ships both engines, scope + budget search, BLOSUM62 / PAM50 / custom matrices
(Gram-distance penalty) with linear gaps, **position-aware scoring** (``PositionalMatrix``, for
anchor masking + PSSMs), top / all hits, on-demand alignment, parallel batch and batch-of-batches,
auto-selecting pairwise search, a C++ **k-mer seed index** (``KmerIndex``), **control-set E-values**
for significance (see :doc:`evalue`), and **pMHC epitope homology search** with presentation-aware
E-values and neoantigen mimic discovery (see :doc:`pmhc`). Planned next:

1. **Alignment polish** — affine gaps (gap open vs extend), CIGAR-style output, batch alignment.
2. **Graded / native pMHC scoring** — ``PositionalMatrix`` already supports per-position PSSMs and
   graded TCR-hotspot weights; remaining work is a native local/best-window engine mode and
   predictor-driven proteome-to-presented-peptide pipelines for full mimic scans.
3. **Significance, continued.** Control-set + presentation-aware E-values ship; a *tf-idf redundancy
   weighting* — down-weighting hits to common / clustered motifs without an external control set — is
   still under consideration.
4. **Memory** — succinct (LOUDS) trie / packed reference strings if reference counts reach the
   10M+ range; memory-mapped frozen index for zero-copy load.
5. **Batch-of-batches & distribution** — streaming query batches and optional process-level fan-out.

Downstream consumers (built in *other* libraries, not here): UMI collapsing, CDR3 nucleotide error
correction, VDJdb CDR3 matching, IEDB epitope matching. seqtree stays payload-agnostic; those
libraries apply V-gene / MHC / count filters on top of ``ref_id`` results.

Benchmark findings
------------------

The simple max-edit-3 benchmark (``bench/bench_gnuplot.py``) makes the cost structure clear:
throughput is dominated by **scope**, not reference-set size. ``seqtm`` substitution-only search
stays cheap; allowing insertions/deletions widens the branch-and-bound frontier sharply, and the
``seqtrie`` banded DP pays a larger constant per candidate. Crucially, enumeration cost depends on
**reference-set redundancy**: a low-redundancy (uniform-random) set has a bushy trie with no shared
prefixes, so indel-heavy scopes blow up; real, clustered repertoires share structure and the trie
collapses, keeping the same scope tractable. This shapes the per-domain direction below.

Per-domain search strategy (open questions)
-------------------------------------------

The three target domains have very different null structure and need different indexing:

- **CDR3 (amino acid).** Index from the **3' / J-proximal end** rather than the 5' end: the J segment
  is longer and more conserved than the V-side N-region, so suffix-anchored tries share more prefix
  structure and prune faster. Reverse the sequence before insertion (and the query before search).
  See also https://arxiv.org/abs/2604.26190 for an alternative formulation.
- **UMIs / nucleotide.** Largely a birthday-problem + Hamming-tree regime: short fixed-length tags,
  substitution-dominated errors, no indels in the common case. A plain Hamming trie with a small
  ``max_subs`` budget is sufficient; the interesting part is collision probability, not search.
- **Epitopes (MHC-presented peptides).** The hardest case — the null is not uniform over peptide
  space but shaped by presentation. Need to characterize the real structure of MHC-I and MHC-II
  ligand repertoires (from eluted-ligand / immunopeptidomics data) or adopt a probabilistic
  presentation model (https://journals.aps.org/prxlife/abstract/10.1103/fbct-vzwm) to define a
  meaningful background before E-values over epitope hits are trustworthy.
