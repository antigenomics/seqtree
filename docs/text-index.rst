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

``TextIndex`` keys on a k-mer seed table over the flat text instead. ``k`` **is a property of
the index, not of the query**, so one build answers every length and every ``max_subs``::

    from seqtree import TextIndex

    ix = TextIndex.build(records, alphabet="aa", k=4)
    res = ix.search_batch(["SIINFEKL", "KTAYIAKQRQISFVKSHFSRQ"], max_subs=2, threads=0)

    for query, hits in zip(queries, res):
        for h in hits:
            h.ref_id, h.offset, h.n_subs, h.mismatches

How it answers exactly
----------------------

One rule, applied per query. Split the query into ``b`` disjoint blocks of width
``bw = L / b >= k``, and probe block *j*'s leading k-mer at radius ``c_j`` — every k-mer within
``c_j`` substitutions of it.

    **Lossless if and only if** ``sum(c_j) >= max_subs - b + 1``.

    The scheme can only miss an occurrence whose errors satisfy ``e_j > c_j`` for *every* block.
    The lightest such error vector is ``e_j = c_j + 1``, of total weight ``sum(c_j) + b``, so no
    occurrence within ``max_subs`` can hide once ``sum(c_j) + b > max_subs``. Errors falling
    outside the blocks only lower the total, so an uncovered tail is free.

Probing one block at radius ``c`` costs ``N(c) = sum(C(k,i)·(A-1)^i, i <= c)``, whose increments
climb steeply in ``c`` and are the same for every block — so the cheapest legal scheme takes as
many blocks as fit and spreads the budget as thinly as possible::

    b = min(max_subs + 1, L // k)
    r = max(0, max_subs - b + 1)          # budget units to place
    c_j = r // b, plus one for the last (r % b) blocks

At ``b = max_subs + 1`` every ``c_j`` is 0 and this is plain pigeonhole — one exact seed per
error, at least one of which must survive. At ``b = 1`` it is a single substitution ball over
``q[0:k]``. **The useful schemes are in between**, and they are where the short queries live: at
``L = 8, max_subs = 2, k = 4`` two disjoint 4-mers fit, so ``c = (0, 1)`` probes **94** variants
where a one-block ball probes **3,267**.

Two details are not incidental. The spare units go on the **last** blocks, worth a further 1.4×:
a candidate already matches its own block's k-mer, so left-to-right verification meets the
*unconstrained* prefix first and rejects after a residue or two — and the high-budget block is
the one contributing nearly all the candidates. And a start that two blocks both reach is
returned **once**, emitted by the lowest-indexed block that could have produced it — a test read
straight off the mismatch positions verification has already computed, so deduplication needs no
sort and no hash set.

None of this is a heuristic. Completeness is checked against a brute-force scan over a grid of
query length, ``max_subs`` and ``k``, asserting **set equality** of ``(ref_id, offset,
n_subs)`` and of the mismatch detail — not merely that hits were found — and separately that the
answer is *identical* across ``k``, so a dispatch bug cannot pass by agreeing with itself.

Choosing ``k``
--------------

``k`` never changes the answer, only the work: a bigger ``k`` makes one probe more specific
(each 4-mer bucket on the human proteome holds ~208 positions, each 5-mer bucket ~9) but needs
``L >= b·k`` to fit the same number of blocks. Measured on the human proteome (UP000005640,
147,506 records / 69,578,135 residues), one thread, Apple M-series,
``bench/bench_text_index.py``:

======  ==========  ===================  ===================
``L``   ``m``       ms/query, ``k = 4``  ms/query, ``k = 5``
======  ==========  ===================  ===================
8       2           **1.468**            4.694
9       2           **1.334**            4.526
10      2           1.298                **0.129**
11      2           1.138                **0.108**
12      2           **0.097**            0.126
15      2           0.086                **0.014**
10      3           3.595                **0.347**
12      3           1.567                **0.394**
15      3           1.514                **0.138**
======  ==========  ===================  ===================

Hit counts are identical in every cell — the benchmark asserts it rather than reporting it, so a
dispatch that got ``b`` or a budget wrong would fail there on real text and not only in a unit
test.

The crossover is ``L = 2k``: below it only one block fits, so the query pays a full radius-``m``
ball over its leading k-mer and the bigger ``k`` is much worse. **Pick the largest** ``k``
**with** ``2k <= L`` for the bulk of the query set — ``k = 4`` down to length 8, ``k = 5`` from 10 up.
At ``k = 5`` every ``L >= 10`` cell above is **at or under 0.4 ms/query**. Building is cheap
either way (0.55 s at ``k = 4``, 0.59 s at ``k = 5``), so a corpus spanning both can simply hold
both indexes — still one build per ``k``, never one per query length.

Threads scale well — human proteome, ``L = 12``, ``max_subs = 2``, 50,000 queries: 0.103
ms/query on one thread, 0.013 on eight (7.8x), 0.008 on all cores (12.6x). The ratio is capped
by the serial CSR flatten at the end of a batch rather than by the search, which is why it
falls short of linear on a workload the search scheme made this fast.

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
costs one copy of the index rather than one per worker — and, since there is nothing to parse,
the load does not scale with the file.

============  =====  ============  ===========  =========  ================  ===============
text          ``k``  residues      file         × text     save              load
============  =====  ============  ===========  =========  ================  ===============
mouse         4      23,131,234    116.6 MB     5.0        0.14 s            0.10 ms (mmap)
human         4      69,578,135    348.2 MB     5.0        0.46 s            0.10 ms (mmap)
human         5      69,578,135    378.1 MB     5.4        0.54 s            0.10 ms (mmap)
============  =====  ============  ===========  =========  ================  ===============

Reading the file instead (``mmap=False``) costs 0.12 s on the human index. **Roughly 5× the
text** is what a *positional* seed index costs: 79 % of the file is ``post_ids``, one ``uint32``
per in-record k-mer start, and every one of those 69.1 M positions has to be addressable for
the answer to be exact. Small texts pay a different bill — the bucket table is directly
addressed, so it has a floor of ``A^k × 4`` bytes (1.3 MB at ``k = 4``, 31.8 MB at ``k = 5``)
regardless of how little text there is.

A ``.sti`` is treated as untrusted input: the header's record, bucket and posting counts are
cross-checked against each other and against the file length before any of them is used to size
an allocation or to index into the mapping.

Insertions and deletions
------------------------

``max_indels`` allows gaps. The two caps are independent, so ``max_subs=2, max_indels=1``
accepts two substitutions **and** one gap — not three edits of any kind:

.. code-block:: python

   res = ix.search_batch(peptides, max_subs=1, max_indels=1, threads=0)
   for h in res[0]:
       print(h.ref_id, h.offset, h.length, h.n_subs, h.n_ins, h.n_dels)
       print(ix.ref_seq(h.ref_id)[h.offset : h.offset + h.length])

**The seed table does not change, and neither does the scheme.** With
``b = max_subs + max_indels + 1`` disjoint blocks the total edit count cannot reach the block
count, so the same pigeonhole argument still puts a **zero-error** block somewhere — and a
zero-error block matches the text exactly, gaps or not, so the existing exact lookup finds it.
Nothing about the index, the build, or the probe enumeration is different.

What changes is verification. A seed no longer pins the start: the prefix before the block may
have gained or lost up to ``max_indels`` residues, so every start in
``[pos - base - max_indels, pos - base + max_indels]`` is aligned, by a banded DP over the
states ``(query position, net offset, indels spent)``.

Two consequences worth knowing before you call it:

* **A match is no longer** ``len(query)`` **residues wide.** Read
  :attr:`~seqtree.TextHit.length`; it is ``len(query) + n_dels - n_ins``.
* **Every seed block must be exact**, so a query must be at least
  ``(max_subs + max_indels + 1) * k`` long. At ``k = 4``, ``max_subs=1, max_indels=1`` needs 12
  residues. A shorter query raises, naming the bound — it is not answered partially. Build a
  second index at a smaller ``k`` if your queries are shorter than that.

An alignment must **begin and end on an aligned pair**. A gap at either edge would stretch the
reported interval over a residue that matches nothing, and a leading gap is just the same
occurrence starting one residue over — allowing it would report one match two or three times.

``mismatches`` is empty on this path: recovering which columns were substituted needs an
alignment traceback the verifier does not keep. The counts are exact; the per-column detail is
only available at ``max_indels = 0``. ``matrix=`` is likewise ignored.

Cost, measured on a synthetic 5,000,000-residue corpus, ``k = 4``, one thread, against the same
query set with ``max_indels = 0``:

.. list-table::
   :header-rows: 1

   * - query length
     - ``max_subs``
     - ms/query, no indels
     - ms/query, ``max_indels=1``
   * - 16
     - 1
     - 0.001
     - 0.028
   * - 16
     - 2
     - 0.001
     - 0.043
   * - 20
     - 3
     - 0.002
     - 0.065
   * - 24
     - 3
     - 0.002
     - 0.067

Gapped search costs 25–37× an ungapped one at the same length and ``max_subs``: the banded DP
replaces a byte comparison that early-exits after a residue or two, and it runs once per start
in the window. It is still well under 0.1 ms/query. Completeness is pinned the same way the
substitution path is — **set equality against an independent brute force**, over
``k`` ∈ {3, 4, 5} × ``max_subs`` 0–3 × ``max_indels`` 1–2, 6,732 occurrences, zero missing and
zero extra.

Limits
------

``TextIndex`` answers a Hamming question over a text. Gaps, scored alignment and edit distance
stay with :class:`Index` and :mod:`seqtree.gapblock` — indels change the length and break the
pigeonhole argument the exactness rests on. Postings are 32-bit, capping the text at 2\\
:sup:`32` symbols; ``k`` is capped so the direct-addressed table stays under 2\\ :sup:`27`
buckets (``k <= 5`` for amino acids).
