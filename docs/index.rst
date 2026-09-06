seqtree
=======

.. raw:: html

   <div class="proj-intro">
     <div>
       <p class="proj-intro__eyebrow">FUZZY SEQUENCE SEARCH</p>
       <p class="proj-intro__lead">seqtree finds biological sequences (amino-acid or nucleotide)
       within a fixed edit scope or score budget. Build an immutable index once, then search
       single queries or millions of queries in parallel. C++ core, minimal Python binding.</p>
       <p class="proj-intro__links">
         <a href="getting-started.html">Getting started</a>
         <span>&middot;</span>
         <a href="engines.html">Engines</a>
         <span>&middot;</span>
         <a href="api.html">API</a>
       </p>
     </div>
   </div>

   <div class="proj-card-grid">
     <a class="proj-card" href="getting-started.html">
       <h3>Getting Started</h3>
       <p>Install, build an index, run your first search.</p>
     </a>
     <a class="proj-card" href="engines.html">
       <h3>Engines &amp; Concepts</h3>
       <p>seqtm vs seqtrie, scope vs budget, scoring.</p>
     </a>
     <a class="proj-card" href="api.html">
       <h3>API Reference</h3>
       <p>Index, matrices, gap-block scoring, E-values.</p>
     </a>
     <a class="proj-card" href="text-index.html">
       <h3>Text Search</h3>
       <p>k-mismatch search over a proteome, one index for every query length.</p>
     </a>
   </div>

   <div class="proj-feature-grid">
     <div class="proj-feature">
       <h3>seqtm &mdash; branch-and-bound</h3>
       <p>Exact per-type edit caps (subs / ins / dels), a fast Hamming-only path, and an exact
       edit-type breakdown per hit. The workhorse for small edit distances: UMI collapse,
       CDR3 error correction, CDR3/epitope matching.</p>
     </div>
     <div class="proj-feature">
       <h3>seqtrie &mdash; banded DP</h3>
       <p>Matrix-weighted score budgets (BLOSUM62 + gap costs) with cost independent of the
       edit count. Best for similarity-scored searches over a total-edit or penalty budget.</p>
     </div>
     <div class="proj-feature">
       <h3>TextIndex &mdash; search a proteome</h3>
       <p>Exact k-mismatch search over a concatenated <em>text</em>, where a trie over reference
       strings would need one build per query length. <strong>k belongs to the index</strong>, so
       one build answers every length and every <code>max_subs</code>: a 9-mer within 2
       substitutions of the human proteome in <strong>1.3 ms</strong> on one thread.</p>
     </div>
     <div class="proj-feature">
       <h3>Pairwise alignment, no BioPython</h3>
       <p>Needleman&ndash;Wunsch and Smith&ndash;Waterman with affine or linear gaps, verified
       against <code>Bio.Align.PairwiseAligner</code> as an oracle with zero disagreements &mdash;
       and 65&ndash;87&times; faster.</p>
     </div>
     <div class="proj-feature">
       <h3>Gap blocks for V(D)J junctions</h3>
       <p>One contiguous indel, its position set by a prior rather than by the score alone &mdash;
       exactly optimal against unrestricted affine alignment on <strong>98.8%</strong> of related
       pairs. <code>IslandProfile</code> adds a per-island PWM over the same scale.</p>
     </div>
     <div class="proj-feature">
       <h3>Calibrated cutoffs, not fixed ones</h3>
       <p>E-values counted against a background control, inverted into the score cutoff that
       achieves a target false-positive rate <strong>per query</strong>. A control repertoire is
       dense near germline and sparse among rare junctions, so one fixed threshold buys a common
       query far more chance neighbours than a rare one.</p>
     </div>
   </div>

Every number on this site is reproducible: see :doc:`benchmarks` for throughput, thread scaling,
recall against ground truth, and the scripts that produce each figure.

Results are payload-agnostic — ``(ref_id, score, n_subs, n_ins, n_dels)``. Downstream libraries
map ``ref_id`` back to their own payloads (V gene, MHC, read counts) and filter there.

.. toctree::
   :hidden:
   :caption: Start here
   :maxdepth: 2

   getting-started
   engines
   examples

.. toctree::
   :hidden:
   :caption: Searching
   :maxdepth: 2

   text-index
   pairwise
   distance
   gapblock

.. toctree::
   :hidden:
   :caption: Significance
   :maxdepth: 2

   evalue
   pmhc

.. toctree::
   :hidden:
   :caption: Reference
   :maxdepth: 2

   api
   benchmarks
   roadmap
