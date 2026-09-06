API Reference
=============

.. currentmodule:: seqtree

Grouped by the job you are doing. For a starting point: :class:`Index` searches a *set of
sequences*, :class:`TextIndex` searches a *long text*, and :mod:`seqtree.pairwise` aligns *two
sequences*.

Searching a set of sequences
----------------------------

The core loop: build an :class:`Index` from your references, describe the search with
:class:`SearchParams`, get back :class:`Hit` objects carrying a ``ref_id`` you map to your own
payload.

.. autoclass:: Index
   :members:
   :undoc-members:

.. autoclass:: SearchParams
   :members:
   :undoc-members:

.. autoclass:: Hit
   :members:
   :undoc-members:

.. autoclass:: Alignment
   :members:
   :undoc-members:

Searching a long text
---------------------

For a short query against a proteome or genome, where enumerating every window as its own
reference is not affordable. One index answers every query length -- see :doc:`text-index`.

.. autoclass:: TextIndex
   :members:
   :undoc-members:

.. autoclass:: TextResult
   :members:
   :undoc-members:

.. autoclass:: TextHit
   :members:
   :undoc-members:

.. autoclass:: ArrayView
   :members:
   :undoc-members:

Scoring
-------

How a mismatch is priced. Penalties are non-negative and zero on a match, so a score is a
distance: lower is better, everywhere in seqtree.

.. autoclass:: SubstitutionMatrix
   :members:
   :undoc-members:

.. autoclass:: PositionalMatrix
   :members:
   :undoc-members:

Dense matrices
--------------

Every query against every reference in one GIL-released call, for when nothing can be pruned.
Returned by :func:`seqtree.pairwise.score_matrix`, :func:`seqtree.distance.hamming_matrix`, and
:func:`seqtree.gapblock.score_matrix`.

.. autoclass:: ScoreMatrix
   :members:
   :undoc-members:

Seed-and-extend
---------------

Candidate generation at million scale: match query k-mers, merge the posting lists, rank. Used
by :mod:`seqtree.pmhc`.

.. autoclass:: KmerIndex
   :members:
   :undoc-members:

.. autoclass:: Candidate
   :members:
   :undoc-members:

Alphabets and batch helpers
---------------------------

.. autofunction:: alphabet_symbols

.. autofunction:: amino_acids

.. autofunction:: pairwise_batch

Significance -- is this hit real?
---------------------------------

A score alone is not evidence: a germline-adjacent query collects neighbours by chance. These
count the same ball in a background control and report an E-value, or invert it for the score
cutoff that achieves a target false-positive rate. See :doc:`evalue`.

.. autofunction:: evalues

.. autofunction:: load_control

.. autofunction:: threshold_for_evalue

.. autofunction:: thetas_from_scores

.. autofunction:: seqtree.evalue.evalue_result

Pairwise alignment (Needleman-Wunsch / Smith-Waterman)
------------------------------------------------------

.. automodule:: seqtree.pairwise
   :members:
   :undoc-members:
   :show-inheritance:

Plain edit distances (Hamming / Levenshtein)
--------------------------------------------

.. automodule:: seqtree.distance
   :members:
   :undoc-members:
   :show-inheritance:

Gap-block alignment
-------------------

.. automodule:: seqtree.gapblock
   :members:
   :undoc-members:
   :show-inheritance:
   :exclude-members: ScoreMatrix

Seed E-values
-------------

.. automodule:: seqtree.seeds
   :members:
   :undoc-members:
   :show-inheritance:

Layout and anchors
------------------

.. automodule:: seqtree.layout
   :members:
   :undoc-members:
   :show-inheritance:

Epitope (pMHC) search
---------------------

.. automodule:: seqtree.pmhc
   :members:
   :undoc-members:
   :show-inheritance:
