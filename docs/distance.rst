Plain edit distances: Hamming and Levenshtein
==============================================

:mod:`seqtree.distance` is the *unweighted* corner of the library: Hamming and Levenshtein
distances on raw characters, unit costs, no substitution matrix and no alphabet. When all you
need is "how many edits apart are these two strings", you should not have to build a
:class:`~seqtree.SubstitutionMatrix` or add ``python-Levenshtein`` / ``rapidfuzz`` as a
dependency. seqtree still needs nothing at runtime.

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

When to use which
-----------------

* Fixed-length tags, substitution-only errors (UMI collapse, barcode demultiplexing) → **Hamming**.
* Mixed lengths, indels in play → **Levenshtein**.
* A **weighted** alignment — a substitution matrix, affine gaps, local mode → :doc:`pairwise`, not
  this module. These two are deliberately unweighted; ``hamming`` with an identity matrix is not
  what BLOSUM-scored search means.
* A large fuzzy search under an edit *budget* rather than a full distance matrix →
  :doc:`the search engines <engines>`, which prune instead of scoring every pair.

See also
--------

* :doc:`pairwise` — weighted Needleman–Wunsch / Smith–Waterman on a substitution matrix.
* :doc:`engines` — indexed fuzzy search under an edit-scope or score budget.
