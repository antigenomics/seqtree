Engines & Concepts
==================

One immutable trie, two search drivers. Pick one with ``engine=`` or let
``engine="auto"`` choose per query.

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

``engine="auto"`` routes substitution-only / small-k indel work to ``seqtm`` and matrix-weighted
budgets to ``seqtrie``.

Scoring
-------

Scores are non-negative penalties where ``0`` is an exact match. Similarity matrices such as
BLOSUM62 are converted at load time via the Gram→squared-distance transform
``pen[a][b] = sim[a][a] + sim[b][b] - 2*sim[a][b]`` — i.e. ``‖φ(a) - φ(b)‖²`` if the score is read
as an inner product ``sim(a,b) = ⟨φ(a), φ(b)⟩``. It is symmetric, ``0`` on the identity, and
non-negative for BLOSUM/PAM (the diagonal is each row's maximum), so every edit adds a non-negative
cost and the budget prune stays valid. Note penalties are roughly twice the scale of the raw score
gaps, so in matrix mode set ``gap_open`` / ``max_penalty`` accordingly. With no matrix, the cost is
unit (1 per substitution and per gap position; linear gaps in v1).

Alphabets
---------

``"aa"`` (20 amino acids in BLOSUM62 order plus ``B Z X *``), ``"nt"`` (``ACGT``), and
``"iupac"`` (nucleotide ambiguity codes). Encoding is case-insensitive; a symbol outside the
alphabet raises an error at build or search time.
