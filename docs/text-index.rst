Text search
===========

.. currentmodule:: seqtree

:class:`Index` answers a question about a *set of reference strings*. :class:`TextIndex`
answers a question about a **text** — a proteome, a genome, a set of transcripts — where the
query is short and the reference is long and has no length of its own:

    Given a reference text and a large batch of short queries **of many different lengths**,
    find every position where a query matches within ``m`` substitutions. No gaps, full
    length, exact answer.

This is peptide-to-proteome mapping. It is not homology search: no scoring matrix in the
predicate, no E-value, no local alignment.

Why not ``Index``
-----------------

To ask a text question with :class:`Index` you have to enumerate every length-``L`` window of
the text as a separate reference string — and do it again for every distinct query length. The
human proteome has 68,389,335 nine-mer windows, so one such index costs several gigabytes, and
a query set spanning 45 distinct lengths costs 45 of them.

``TextIndex`` keys on a k-mer seed table over the flat text instead. **``k`` is a property of
the index, not of the query**, so one build answers every length and every ``max_subs``::

    from seqtree import TextIndex

    ix = TextIndex.build(records, alphabet="aa", k=4)
    res = ix.search_batch(["SIINFEKL", "KTAYIAKQRQISFVKSHFSRQ"], max_subs=2, threads=0)

    for query, hits in zip(queries, res):
        for h in hits:
            h.ref_id, h.offset, h.n_subs, h.mismatches

How it answers exactly
----------------------

Two paths share the one seed table, chosen per query on ``s = L / (max_subs + 1)``:

``s >= k`` — **pigeonhole.** Split the query into ``max_subs + 1`` disjoint blocks. A match
within ``max_subs`` substitutions must leave at least one block untouched, so at least one
block's leading k-mer is exact. Probe all of them and verify the candidates.

``s < k`` — **ball.** Enumerate the ``<= max_subs`` neighbourhood of the query's *first* ``k``
residues and probe every variant. Lossless for the same reason: a match carrying ``<= m``
mismatches carries at most ``m`` of them in its first ``k`` residues. Under proper
substitutions the enumeration is duplicate-free by construction, so it needs no sort and no
hash set.

Neither path is a heuristic. Completeness is checked against a brute-force scan over a grid of
query length, ``max_subs`` and ``k``, asserting **set equality** of ``(ref_id, offset,
n_subs)`` and of the mismatch detail — not merely that hits were found.

Choosing ``k``
--------------

``k`` sets how specific one seed probe is, so it trades the seed path against the ball path.
Measured on the human proteome (UP000005640, 147,506 records / 69,578,135 residues), one
thread, Apple M-series, ``bench/bench_text_index.py``:

======  ==========  ===================  ===================
``L``   ``m``       ms/query, ``k = 4``  ms/query, ``k = 5``
======  ==========  ===================  ===================
8       2           29.5                 **4.24**
9       2           34.5                 **4.86**
10      2           32.8                 **4.49**
11      2           27.3                 **3.26**
12      2           **0.139**            3.92
======  ==========  ===================  ===================

Hit counts are identical in every cell — ``k`` never changes the answer, only the work. The
step at ``L = 12`` is the dispatch: ``12 / (2 + 1) = 4``, so ``k = 4`` takes the seed path and
``k = 5`` does not. Building is cheap either way (**0.59 s** at ``k = 4``, **0.67 s** at
``k = 5``, ~750 MB peak), so a mixed corpus can afford one index per ``k`` — which is still one
build per ``k``, not one per query length.

Rule of thumb: **pick the largest ``k`` with ``k <= L / (max_subs + 1)``** for the bulk of the
query set, and accept the ball path for the tail. At ``L >= 12, m <= 3`` a ``k = 4`` index
answers in **0.12–0.20 ms/query** on the human proteome; the short-query ball path is one to
two orders of magnitude more expensive because each probe of a 4-residue seed returns ~210
positions on a text that size, and verifying a candidate is one cache miss.

Threads scale nearly linearly — human proteome, ``L = 12``, ``max_subs = 2``: 0.138 ms/query on
one thread, 0.017 on eight (8.3x), 0.010 on all cores (13.9x).

Results are flat arrays
-----------------------

A 445,000-query run returns a handful of arrays, not millions of Python objects.
:class:`TextResult` holds parallel arrays in CSR form and materialises hit objects only for
the query you ask for::

    res.num_hits                 # total across the batch
    res[i]                       # TextHit list for query i, built on demand
    res.arrays()                 # every array as a zero-copy ArrayView (buffer protocol)
    res.to_numpy()               # the same, as numpy views; numpy imported on the call

Hits for query ``i`` are ``[query_begin[i], query_begin[i+1])``; mismatches for hit ``h`` are
``[mm_begin[h], mm_begin[h+1])``. Order is ``(n_subs, ref_id, offset)`` — stable across runs
and across thread counts, so a downstream digest is stable.

Beyond a position
-----------------

**Mismatch detail is a pair, not a position.** Each mismatch is ``(pos, query_aa, text_aa)``,
so a caller ranking by chemistry can tell ``L→I`` from ``L→D`` without re-fetching the window.

**Optional scoring.** Pass ``matrix=`` and each hit also carries a ``score``: the substitution
matrix's similarity summed over the mismatched positions, so ``0`` is an exact match and
less-negative is more conservative. It only *scores* hits the Hamming predicate has already
accepted — it never changes which hits come back.

**Nearest shell, whole.** ``best_only=True`` walks the distance upward and stops at the first
shell that has anything in it, returning **all** of it. A query resolved at ``m = 0`` never
enters the ``m = 1`` pass, and no second index is needed.

**Group folding.** ``build(..., group_ids=[...])`` labels each record with an arbitrary integer
— gene ids, species, cluster labels; seqtree does not know what they mean. Then
``group_by=True`` collapses hits onto them::

    res.groups(i)   # [(group_id, min_subs, n_hits), ...], sorted by group id

More than one row means the nearest parents disagree, which makes a tie a first-class output
rather than something every caller re-derives from the raw hits.

**Truncation is reported.** A degenerate query can match a repeat region a very large number of
times, so ``max_hits`` caps it — *after* the sort, so the best hits survive — and sets
``res.truncated[i]``. A cap that is invisible is a recall bug wearing a performance costume.

The alphabet, which is where a silent wrong answer would live
-------------------------------------------------------------

Symbols are compared through the alphabet's codec, not as raw bytes.
For ``alphabet="aa"`` that codec is the 24 symbols ``ARNDCQEGHILKMFPSTWYVBZX*``, so:

* ``B``, ``Z``, ``X`` and ``*`` are **real symbols, not wildcards** — ``X`` matches ``X`` and
  is one substitution from ``A``. The human proteome carries 8,417 ``X`` residues and the mouse
  proteome 4,263, so this is not hypothetical.
* ``U`` (selenocysteine) is **not in the codec**. It encodes to the invalid marker, so no hit
  can cross a ``U`` in the text, and a *query* containing one is refused by name::

      ValueError: queries[42] ('CASSUGQYF'): symbol 'U' is not in the alphabet

  The same marker separates records, which is what structurally prevents a hit — or a seed —
  from spanning a record boundary, rather than a check that could be forgotten.

Every query must be at least ``k`` residues long. Seeds are indexed only where a whole k-mer
fits inside a record, so a shorter query could match within ``k - L`` of a record end and be
missed; ``TextIndex`` raises rather than answering incompletely.

Persistence
-----------

The file is flat, with the arrays at known aligned offsets, so loading it is a map rather than
a parse::

    ix.save("human_k4.sti")
    ix = TextIndex.load("human_k4.sti")            # mmap=True by default

With ``mmap=True`` the pages are shared read-only across processes, so a fan-out of workers
costs one copy of the index rather than one per worker.

Limits
------

``TextIndex`` answers a Hamming question over a text. Gaps, scored alignment and edit distance
stay with :class:`Index` and :mod:`seqtree.gapblock` — indels change the length and break the
pigeonhole argument the exactness rests on. Postings are 32-bit, capping the text at 2\\
:sup:`32` symbols; ``k`` is capped so the direct-addressed table stays under 2\\ :sup:`27`
buckets (``k <= 5`` for amino acids).
