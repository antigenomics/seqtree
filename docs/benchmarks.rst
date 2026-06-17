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

Indicative numbers
-------------------

On an Apple M3 (16 cores), 1M amino-acid references, single-substitution queries:

.. list-table::
   :header-rows: 1

   * - metric
     - seqtm
     - seqtrie
     - notes
   * - build (1M)
     - ~0.8 s
     - ~0.8 s
     - shared trie
   * - peak RSS
     - ~1.5 GB
     - ~1.5 GB
     - well under 32 GB
   * - throughput
     - ~0.5 M q/s
     - ~0.1 M q/s
     - 16 threads, k=1
   * - alignment
     - ~0.3 µs
     - n/a
     - per call, on demand

``seqtm`` is markedly faster at small edit distances; ``seqtrie`` is the choice for
matrix-weighted budgets. Numbers vary with sequence length, scope, and hit density.
