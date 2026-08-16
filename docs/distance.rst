Plain edit distances: Hamming and Levenshtein
==============================================

:mod:`seqtree.distance` is the *unweighted* corner of the library: Hamming and Levenshtein
distances on raw characters, unit costs, no substitution matrix and no alphabet. When all you
need is "how many edits apart are these two strings", you should not have to build a
:class:`~seqtree.SubstitutionMatrix` or add ``python-Levenshtein`` / ``rapidfuzz`` as a
dependency. seqtree still needs nothing at runtime. The same module also *enumerates* a Hamming
ball — the members, not the distance — and takes the deduplicated union over many centres.

.. contents::
   :local:
   :depth: 1

Two distances
-------------

.. code-block:: python

   from seqtree.distance import hamming, levenshtein

   hamming("CASSLGQYF", "CASSPGQYF")          # 1  -- differing positions, equal length only
   levenshtein("kitten", "sitting")           # 3  -- insertions + deletions + substitutions
   levenshtein("CASSLGQAYEQYF", "CASSPGQAYEQF")  # 2

* **Hamming** counts the positions at which two **equal-length** sequences differ. On a length
  mismatch it raises :class:`ValueError` — the distance is simply undefined there, and silently
  returning something would hide a bug in the caller.
* **Levenshtein** is the classic edit distance: the fewest single-character insertions, deletions
  and substitutions that turn one string into the other, each costing 1. Any lengths.

Comparison is **case-sensitive**, byte for byte. This is the one place the library does *not* fold
case — the search engines treat ``a`` and ``A`` as equal, but a generic string distance should
report the difference it is asked about.

Matrices, in parallel
----------------------

Every ``a`` against every ``b`` in one GIL-released, multi-threaded C++ call, returned as a
:class:`~seqtree.ScoreMatrix` that ``numpy.asarray`` wraps without copying:

.. code-block:: python

   import numpy as np
   from seqtree.distance import hamming_matrix, levenshtein_matrix

   umis = ["ACGTACGT", "ACGTACGA", "TTTTACGT"]
   d = np.asarray(hamming_matrix(umis, umis, threads=0))   # 0 = all cores
   d.shape        # (3, 3)
   d.diagonal()   # all zero -- a sequence's distance to itself

   d = np.asarray(levenshtein_matrix(cdr3s, prototypes, threads=0))

:func:`~seqtree.distance.hamming_matrix` raises :class:`ValueError` if any pair has mismatched
lengths, so it is the right tool for a set of fixed-length tags (UMIs, barcodes, one-length CDR3s);
:func:`~seqtree.distance.levenshtein_matrix` places no such constraint.

Enumerating a Hamming ball
--------------------------

The functions above *score* a pair you already hold. :func:`~seqtree.distance.neighbourhood`
*generates* one instead — every sequence within ``r`` substitutions of a centre:

.. code-block:: python

   from seqtree.distance import neighbourhood, neighbourhood_union, union_size

   neighbourhood("CASSLGQYF")                     # 172 sequences = 19*9 + 1
   neighbourhood("CASSLGQYF", include_self=False) # 171
   neighbourhood("CASSLGQYF", 2, shell=True)      # 12,996 = 19^2 * C(9,2), distance exactly 2
   neighbourhood("A", 1, alphabet="ACGT")         # ['A', 'C', 'G', 'T']

**Substitution only**, so every member has the length of the centre — Hamming distance is
undefined across lengths, and an indel is a different question (:doc:`gapblock`). Over a
``k``-letter alphabet the closed ball holds ``sum((k-1)**d * comb(L, d) for d in range(r+1))``;
``alphabet=None`` means the 20 standard residues, i.e. :func:`~seqtree.amino_acids` **minus** the
ambiguity codes ``B``/``Z``/``X`` and the stop ``*``, which that function does include.

:func:`~seqtree.distance.neighbourhood_union` is the one that earns its keep. The union of many
balls is not their concatenation, and near-duplicate centres — co-specific TCR junctions, an error
family around one UMI — overlap heavily:

.. code-block:: python

   junctions = [...]                              # a specificity group
   union_size(junctions)                          # how big before you commit to it
   for seq in neighbourhood_union(junctions):     # each distinct sequence exactly once
       ...

Each sequence is emitted once, deduplicated by a multi-source breadth-first walk *during*
generation — the ``sum(19*L_i)`` multiset is never materialised. ``shell=True`` returns the
members whose distance to the **nearest** centre is exactly ``r``; ``include_self=False`` drops
the whole ``r = 0`` shell, which for a union is exactly ``set(seqs)`` — including a centre that
happens to be one substitution from another.

How much the dedup buys depends on how tight the group is. 200 junctions of length 14, all
drawn within distance ``d`` of a common centre, ``r = 1``:

.. list-table::
   :header-rows: 1

   * - spread ``d``
     - sum of the balls
     - union
     - double-counted
   * - 1
     - 53,400
     - 31,122
     - 41.7%
   * - 2
     - 53,400
     - 50,971
     - 4.5%
   * - 3
     - 53,400
     - 53,312
     - 0.2%

At the scale a precursor-frequency calculation asks for — 300 junctions of length 14 at ``r = 1``,
80,100 distinct sequences — the whole union takes **23 ms** on one M3 core (``union_size`` alone,
which skips the sort and the result list, 11 ms). ``r = 2`` over the same 300 is 9.9 M sequences,
6.8 s and ~1.8 GB: there the dedup is a memory question, not a style one, and the walk is the only
thing keeping it to one copy.

When to use which
-----------------

* Fixed-length tags, substitution-only errors (UMI collapse, barcode demultiplexing) → **Hamming**.
* Mixed lengths, indels in play → **Levenshtein**.
* A **weighted** alignment — a substitution matrix, affine gaps, local mode → :doc:`pairwise`, not
  this module. These two are deliberately unweighted; ``hamming`` with an identity matrix is not
  what BLOSUM-scored search means.
* A large fuzzy search under an edit *budget* rather than a full distance matrix →
  :doc:`the search engines <engines>`, which prune instead of scoring every pair.
* The ball's *members*, not a distance — every variant you must then score with something that
  cannot be indexed (a ``P_gen`` model, an external predictor) →
  :func:`~seqtree.distance.neighbourhood_union`. If the thing you are searching is already in an
  :class:`~seqtree.Index`, search it instead: enumerating 19·L candidates to look each one up is
  what the trie exists to avoid.

See also
--------

* :doc:`pairwise` — weighted Needleman–Wunsch / Smith–Waterman on a substitution matrix.
* :doc:`engines` — indexed fuzzy search under an edit-scope or score budget.
