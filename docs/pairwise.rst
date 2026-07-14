Pairwise alignment without BioPython
=====================================

:mod:`seqtree.pairwise` is Needleman–Wunsch and Smith–Waterman on the raw log-odds scale — an
ordinary protein aligner, so that reaching for BioPython is no longer necessary just to score two
sequences against each other.

.. contents::
   :local:
   :depth: 1

Two scales, one matrix
----------------------

Everything else in seqtree **minimises a non-negative penalty**. That is what a search ball and an
E-value need: a distance-like cost that is zero when two sequences are identical.

An alignment score is the opposite. It **maximises a signed similarity**, the way BLAST and
BioPython do, because a positive score is what "these residues are more alike than chance" means.

Both live on the same :class:`~seqtree.SubstitutionMatrix`:

.. code-block:: python

   mat = seqtree.SubstitutionMatrix.blosum62()
   mat.penalty("A", "K")      # 6  -- non-negative, zero on the diagonal
   mat.similarity("A", "K")   # -1 -- the raw log-odds, signed

.. important::

   The penalty is the Gram transform of the similarity, ``pen = s(a,a) + s(b,b) − 2·s(a,b)``, and
   that transform is **lossy**: it forces the diagonal to zero, destroying ``s(a,a)``. BLOSUM62 has
   ``s(C,C) = 9`` and ``s(W,W) = 11``, but both become ``penalty == 0``. A similarity therefore
   *cannot* be recovered from a penalty, which is why the raw grid is stored rather than
   reconstructed.

Aligning
--------

.. code-block:: python

   from seqtree.pairwise import score, align, score_matrix, dist_matrix

   score("CASSLGQAYEQYF", "CASSPGQAYEQF", mat)                  # 45  (global, BLAST defaults)
   score("WWWAAAWWW", "KKKAAAKKK", mat, mode="local")           # 12  (only the AAA core)
   score("AAA", "AAAAA", mat, gap_open=5, gap_extend=5)         # linear gaps

   aln = align("CASSLGQAYEQYF", "CASSPGQAYEQF", mat)
   aln.aligned_query, aln.aligned_ref   # ('CASSLGQAYEQYF', 'CASSPGQAYEQ-F')

Four conventions, each one a place a reimplementation goes quietly wrong:

* a gap run of length ``L`` costs ``gap_open + (L-1)·gap_extend`` — ``gap_open`` is the cost of the
  *first* gap column, not a surcharge on top of it;
* **``gap_open == gap_extend`` gives linear gaps.** There is no separate mode, and none is needed;
* ``mode="global"`` charges end gaps like any other — true Needleman–Wunsch, not semi-global;
* ``mode="local"`` never lets the score fall below zero and takes the best cell anywhere —
  Smith–Waterman.

Gap costs are **positive magnitudes** and are subtracted. BLAST's protein defaults are
``gap_open=11, gap_extend=1``; BioPython's ``PairwiseAligner("blastp")`` preset uses ``12, 1``.

It is a drop-in
---------------

The point of the module is to *be* BioPython for this job, so it is tested against it as an oracle:
three matrices (BLOSUM45/62/80) × ten gap-and-mode settings (affine, linear, free, extreme; global
and local) × sixty sequence shapes (equal, unequal, identical, one-indel, wildly unequal).

**Zero disagreements** — including on real germline V genes, where the quantity of interest is the
distance ``d = s(a,a) + s(b,b) − 2·s(a,b)``.

BioPython is a **test-only** dependency. seqtree itself has zero required runtime dependencies and
never imports it.

And it is faster
----------------

There is no Python in the per-pair loop — a C++ Gotoh with ``O(min(m,n))`` memory, and the batch
paths release the GIL. All-against-all, BLOSUM62, global, gap 11/1, on an M3:

==============================  ==================  ===================  ==============  =========
sequence length                 seqtree, 1 thread   seqtree, 16 threads  BioPython       speedup
==============================  ==================  ===================  ==============  =========
15 (a CDR3 junction)            1.7 M pairs/s       **20.1 M pairs/s**   0.31 M pairs/s  **65×**
90 (a germline V gene)          72 k pairs/s        **893 k pairs/s**    10 k pairs/s    **87×**
==============================  ==================  ===================  ==============  =========

Distances, directly
-------------------

The quantity a prototype-distance embedding actually consumes is not the score but

.. math::

   d(a, b) = s(a,a) + s(b,b) - 2\\,s(a,b)

— the Gram transform applied at the *sequence* level. It is non-negative, symmetric, and zero on
the diagonal, so it is a distance. :func:`~seqtree.pairwise.dist_matrix` computes it without a
Python loop, taking each self-score once per sequence rather than once per pair:

.. code-block:: python

   import numpy as np
   d = np.asarray(dist_matrix(v_genes, v_genes, mat, gap_open=12, gap_extend=1, threads=0))
   d.shape        # (179, 179)
   d.diagonal()   # all zero

Which matrix?
-------------

``BLOSUM45`` (shallow, for remote homologs), ``BLOSUM62`` (the default), and ``BLOSUM80`` (deep, for
close ones) all ship, along with ``PAM250``, ``PAM100`` and ``structural``. They are accepted
wherever a matrix name is, including in :class:`~seqtree.SearchParams`.

See also
--------

* :doc:`gapblock` — the *other* aligner: one contiguous indel, minimising a penalty, for junctions.
* ``bench/bench_pairwise.py`` — reproduces the throughput table above.
