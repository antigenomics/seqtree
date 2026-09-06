Engines & Concepts
==================

One immutable trie, two search drivers. ``engine="auto"`` (the default) always means
``seqtm``; ask for ``seqtrie`` by name when you want it.

Scope vs budget
---------------

A query matches a reference when it satisfies the **scope** and/or the **budget**:

- **Scope** — per-type edit caps: ``max_subs``, ``max_ins``, ``max_dels``, and an optional
  combined ``max_total_edits``. A per-type cap of ``0`` means *zero of that type*.
- **Budget** — a score threshold ``max_penalty`` under a substitution matrix and gap costs.

``max_total_edits`` is an independent total cap (``0`` means "no total cap", falling back to the
per-type sum). It is not clamped by the per-type caps, so ``seqtrie`` can be driven by it alone.

seqtm — branch-and-bound
------------------------

Enumerates each edit (substitution / insertion / deletion) while descending the trie, tracking
the per-type counts, and prunes as soon as a cap or the budget is exceeded. Consequences:

- Per-type caps are enforced **exactly**, and every hit reports an exact
  ``(n_subs, n_ins, n_dels)`` breakdown.
- A dedicated **Hamming-only** path runs when ``max_ins == max_dels == 0``.
- Cost grows with the number of allowed edits, so it is fastest at small distances (k = 1–3) —
  which covers UMI collapse, CDR3 error correction, and CDR3 / epitope matching.

seqtrie — banded DP
-------------------

Carries an edit-distance DP row down the trie and prunes a subtree once its best cell exceeds the
budget. Consequences:

- Handles a matrix-weighted **score budget** (BLOSUM62 or a custom matrix) and indels naturally.
- Cost is independent of the edit count, so it scales better to large budgets / long range.
- It tracks a single cost, so it enforces ``max_total_edits`` + ``max_penalty`` but **not** the
  per-type caps; ``n_subs`` / ``n_ins`` / ``n_dels`` are reported as ``0`` (use
  :meth:`~seqtree.Index.align` to recover the breakdown).

Choosing between them
---------------------

``engine="auto"`` resolves to ``seqtm``, always. It does **not** inspect the query: ``seqtrie``
ignores the per-type caps, so routing a capped search there would silently widen the ball, and a
matrix without an explicit ``max_penalty`` would leave the budget unbounded and scan the whole
index. Neither failure is visible in the results, so ``auto`` never selects it.

The practical consequence: ``seqtrie`` **runs only when you name it.** If you set a matrix and a
``max_penalty`` and leave ``engine`` alone, you get ``seqtm`` — which still enforces ``max_subs``,
and that defaults to ``0``:

.. code-block:: python

   p = seqtree.SearchParams(matrix="blosum62", max_penalty=12)          # engine="auto" -> seqtm
   idx.search(query, p)          # 1 hit: max_subs is still 0, so only exact matches pass

   p = seqtree.SearchParams(matrix="blosum62", max_penalty=12, engine="seqtrie")
   idx.search(query, p)          # 3 hits: the budget is the whole specification

Use ``seqtrie`` when a score budget is the entire specification and you do not care how the
distance decomposes. Use ``seqtm`` — so, the default — every other time.

Scoring
-------

Scores are non-negative penalties where ``0`` is an exact match. Similarity matrices such as
BLOSUM62 are converted at load time via the Gram→squared-distance transform
``pen[a][b] = sim[a][a] + sim[b][b] - 2*sim[a][b]`` — i.e. ``‖φ(a) - φ(b)‖²`` if the score is read
as an inner product ``sim(a,b) = ⟨φ(a), φ(b)⟩``. It is symmetric, ``0`` on the identity, and
non-negative for BLOSUM/PAM (the diagonal is each row's maximum), so every edit adds a non-negative
cost and the budget prune stays valid.

The transform roughly doubles the scale, so gap costs have to move with it.
``matrix.scale()`` reports the median mismatch penalty — 14 for BLOSUM62 — and the rule is
``gap_open = 2 * matrix.scale()``, i.e. **28 for BLOSUM62**. The default of ``1`` is only right
for unit cost; leaving it there under a matrix makes a gap ~14x cheaper than a substitution, and
every alignment degenerates into gaps.

With no matrix the cost is unit: 1 per substitution and 1 per gap position. Search charges gaps
linearly (``gap_extend`` applies to :meth:`~seqtree.Index.align` and the
:mod:`seqtree.pairwise` aligners, not to the trie search).

Alphabets
---------

``"aa"`` (20 amino acids in BLOSUM62 order plus ``B Z X *``), ``"nt"`` (``ACGT``), and
``"iupac"`` (nucleotide ambiguity codes). Encoding is case-insensitive; a symbol outside the
alphabet raises an error at build or search time.
