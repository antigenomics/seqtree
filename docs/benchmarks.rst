Benchmarks
==========

Two harnesses ship with the repo. Both can bootstrap realistic TCR CDR3 sequences from OLGA
(if installed) and otherwise fall back to seeded random sequences.

C++ (raw throughput + scaling)
------------------------------

.. code-block:: console

   cmake -S . -B build -DSEQTREE_BENCH=ON && cmake --build build
   ./build/seqtree_bench 1000 10000 100000 1000000

Reports build time, peak RSS, single-query latency (median / p99), batch throughput, and
thread scaling, followed by a two-engine comparison and per-call alignment cost.

Python (methods + recall)
-------------------------

``bench/bench_methods.py`` compares ``seqtm`` vs ``seqtrie`` across reference sizes, edit
scopes, and budgets (edit count and BLOSUM62 score), plus alignment-fetch cost:

.. code-block:: fish

   python bench/bench_methods.py
   env RUN_BENCHMARK=1 python bench/bench_methods.py --sizes 100000 1000000

``bench/bench.py`` measures **recall** against ground truth on the AIRR VDJdb table (queries are
mutated references with known parents), with throughput and peak RSS:

.. code-block:: fish

   python bench/bench.py
   env RUN_BENCHMARK=1 python bench/bench.py --sizes 1000000 --queries 1000000 --threads 16

E-value benchmark
-----------------

``bench/bench_evalue.py`` is the **true E-value benchmark**. For a target repertoire (VDJdb,
antigen-selected) scored against the ``airr_control`` background, at each scope/budget it reports the
number of **neighbours** (distinct hits, **excluding exact/self matches** — the queries are members
of the target, so the self-match is dropped per the punctured-null lemma), the **exact** self-hits
removed, the number of **collisions** (references re-reached via a different edit path — non-zero only
for ``seqtm`` with indels), and the **fraction of neighbours called significant** both at fixed
E-value cutoffs and after a Benjamini–Hochberg FDR correction across the query family:

.. code-block:: fish

   python bench/bench_evalue.py
   env RUN_BENCHMARK=1 python bench/bench_evalue.py --target-size 200000 --control-size 2000000

The discriminating result is the contrast between query sets: antigen-selected VDJdb queries produce
orders of magnitude more neighbours and are largely significant (BH FDR < 0.05), whereas background
(control) queries produce almost none and survive no correction. The smallest resolvable E-value is
:math:`N/M`, so finer fixed cutoffs (``E < 0.01``) require a control much larger than the target (the
``RUN_BENCHMARK`` tier uses :math:`M = 2{,}000{,}000`); the BH correction is what makes the
fixed-cutoff fractions trustworthy at small control sizes.

TCR-beta benchmark (gnuplot figures)
------------------------------------

``bench/bench_gnuplot.py`` is the main benchmark. It measures **two reference families separately**,
so the effect of sequence structure is visible rather than averaged away:

* **olga** — OLGA-generated human TRB CDR3 (a generative model, **no antigen motif**); queried with
  1000 fresh OLGA TRB sequences.
* **vdjdb** — VDJdb CDR3 **mutated** (real antigen-specific receptors, **with shared motif**
  structure); queried with 1000 held-out VDJdb CDR3.

Both families are expanded by substitution-mutation to each target size. Timings are over the 1000
queries. Figures are **vertically stacked two-panel SVGs**; **seqtm is drawn with a long dash and
seqtrie with a dash-dot**, and the reference family is encoded by colour.

.. code-block:: fish

   python bench/bench_gnuplot.py                       # fast tier: 10k / 100k
   env RUN_BENCHMARK=1 python bench/bench_gnuplot.py   # full tier: 10k / 100k / 1M / 10M

Each figure is written to ``bench/figures/<key>.svg`` (+ per-panel ``.tsv``). Requires ``gnuplot``
and ``olga-generate_sequences`` on PATH (``pip install olga``). The *scaling*, *matrix* and *per-op*
figures span all reference sizes; the edit-budget sweeps (*scope*, *selectivity*, *collisions*) run
at one representative size to stay tractable.

A note on engine semantics: at an edit budget *e*, **seqtm** explores the **Hamming ball**
(``max_subs=e``, substitutions only — the dominant TCR diversity/error mode) while **seqtrie**
explores the **edit-distance ball** (``max_total_edits=e``, substitutions *and* indels). They answer
subtly different questions, which shows as a higher match count for seqtrie at the same *e*.

Scaling and parallelism
~~~~~~~~~~~~~~~~~~~~~~~~~

Throughput (queries per millisecond) versus reference-set size, for both engines at 1, 4, and 8
threads (fixed scope: 2 substitutions), with the olga family on top and vdjdb below. Batches
parallelize near-linearly to 8 cores (~6.5–7×):

.. image:: _static/bench/scaling.svg
   :alt: throughput vs reference size per engine and thread count, olga and vdjdb
   :width: 80%

Edit budget
~~~~~~~~~~~

Cost (top) and selectivity (bottom) as the edit budget grows from 1 to 5. Throughput is governed by
**scope** far more than by reference-set size, and the match count grows steeply — by *e* = 5 a query
already pulls hundreds (seqtm Hamming ball) to thousands (seqtrie edit ball) of neighbours, so loose
budgets are rarely useful:

.. image:: _static/bench/scope.svg
   :alt: throughput and matches per query vs edit budget 1..5
   :width: 80%

Matrix scoring (BLOSUM62 / PAM50 / custom)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

seqtm scores substitutions through a substitution matrix, reporting the best (minimum-penalty) score
across all alignments to each reference. The **time overhead of matrix scoring is small** (within
~5–10 % of unit cost — one table lookup replaces a character compare), for both families:

.. image:: _static/bench/matrix.svg
   :alt: seqtm throughput unit vs BLOSUM62 vs PAM50, olga and vdjdb
   :width: 80%

Besides the built-in ``BLOSUM62`` and ``PAM50``, a custom matrix is supplied via
``SubstitutionMatrix.from_similarity`` (row/column order from ``seqtree.amino_acids()``).

Selectivity and collisions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Top: matches per query versus the seqtrie ``max_penalty`` budget — PAM50 is stricter than BLOSUM62 at
equal budget. Bottom: **collisions** — how often seqtm's branch-and-bound re-reaches the *same*
reference via a *different* edit path (reported by ``Index.collisions_batch``). Substitution-only
search never collides; once indels are allowed, collisions rise with the edit budget, and the
motif-rich **vdjdb** family collides far more than **olga** because shared structure makes a
reference reachable by many distinct edit paths:

.. image:: _static/bench/selectivity_collisions.svg
   :alt: matches per query vs penalty budget, and seqtm collisions per query vs edit budget
   :width: 80%

Per-operation costs
~~~~~~~~~~~~~~~~~~~~~

Top: fetching a global-alignment CIGAR (the C++ Needleman–Wunsch in ``Index.align``) is on-demand and
about a microsecond per call, roughly flat in reference count. Bottom: peak resident memory after the
index build, which scales with the reference count (the trie is shared by both engines):

.. image:: _static/bench/perop.svg
   :alt: align CIGAR fetch cost and peak RSS vs reference size
   :width: 80%

Indicative numbers
~~~~~~~~~~~~~~~~~~~

Apple M3, OLGA TRB references, 1000 queries, 8 threads (``bench/bench_gnuplot.py``, full tier):

.. list-table::
   :header-rows: 1

   * - metric
     - 10k
     - 100k
     - 1M
     - 10M
   * - seqtm, 2 subs (q/ms)
     - ~271
     - ~48
     - ~8.6
     - ~2.4
   * - seqtrie, edits≤2 (q/ms)
     - ~55
     - ~11
     - ~1.2
     - ~0.28
   * - align CIGAR fetch (µs)
     - ~0.9
     - ~1.1
     - ~1.1
     - ~1.8
   * - peak RSS (MB)
     - ~91
     - ~183
     - ~990
     - ~3170

The 8-thread speed-up is ~6.5–7×; matrix scoring stays within ~5–10 % of unit cost across all sizes.

Takeaway
~~~~~~~~

Throughput is governed by **scope** (edit budget) far more than reference-set size, parallelizes
near-linearly to 8 cores, and matrix scoring is nearly free. Sequence structure matters: the
motif-rich vdjdb family is denser and collides more under indels — enumeration cost ultimately
depends on reference redundancy (see :doc:`roadmap`).
