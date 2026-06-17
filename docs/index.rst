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
       <p>Index, SearchParams, Hit, Alignment.</p>
     </a>
     <a class="proj-card" href="benchmarks.html">
       <h3>Benchmarks</h3>
       <p>Throughput, scaling, alignment cost.</p>
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
   </div>

Results are payload-agnostic — ``(ref_id, score, n_subs, n_ins, n_dels)``. Downstream libraries
map ``ref_id`` back to their own payloads (V gene, MHC, read counts) and filter there.

.. toctree::
   :hidden:
   :maxdepth: 2

   getting-started
   engines
   api
   benchmarks
   roadmap
