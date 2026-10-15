API Reference
=============

.. currentmodule:: seqtree

Batch CPU budgets
-----------------

``threads=1`` runs inline. ``threads=0`` selects the detected available CPU budget:
Linux affinity, root cgroup v2 quota and SLURM allocation; macOS physical cores;
Windows process affinity. Legacy nonpositive values also select automatic mode.
Positive counts bound workers; a batch launches at most its query count. A request
above the detected allocation raises ``ValueError``; explicit counts are never clamped.
The caller performs one share of the work: a budget of B uses B compute threads
in total, creating only B-1 additional OS threads when there is enough work.
Pools are joined before the call returns; use one batch rather than an outer pool
of internally threaded searches. Cgroup layouts outside the detected root quota
and Windows CPU-rate job limits require an explicit thread budget.

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

For edit-count summaries, ``Index.edit_histogram_batch(queries, params, threads=0, exclude_exact=False, max_indels=-1)``
returns one list of exact counts per query, indexed by total substitutions,
insertions and deletions. Set ``params.max_total_edits`` explicitly to a positive
value; bins run from zero through that cap. Reference duplicates count separately.
``exclude_exact=True`` retains bin zero with count zero. Nonpositive ``threads``
uses the detected available CPU budget, matching the other batch APIs. This requires ``seqtm`` or ``auto``
and all-hit mode; the search uses one reusable native hit buffer per worker rather
than materializing all hits in Python. Counts are not cumulative; sum the bins up
to a radius to obtain its ball count.
The optional final ``max_indels`` argument limits total insertions plus deletions;
``-1`` preserves the directional-only limits. Bounded unit-cost profiles reuse bit-vector rejection
and exact constrained dynamic programming; other profiles use the generic search
worker. Existing eligible count paths are retained when the bound is redundant.
Explicit ``seqtrie`` remains unsupported.

Floating alignment target reductions
------------------------------------

``Index.alignment_max_batch`` filters unit-cost search candidates by requested
numeric target IDs, aligns each relevant reference sequence once, and reduces
annotation scores without exporting the query-by-reference pair table. It releases
the GIL and uses the same batch CPU budget as search. Build the index with unique
sequences. Bound output memory with contiguous query chunks against that index.
This is a numeric reduction API; callers supply their own group mappings and counts.

The positional arguments are ``queries, params, similarity, similarity_symbols,
max_indels``. ``similarity`` is a square finite float64 grid in the explicit
``similarity_symbols`` order; asymmetric grids are allowed. Every input residue
must be represented. Candidate ``params`` must describe uncapped unit-cost all-hit
``seqtm``/``auto`` search. ``max_indels`` independently caps total insertion and
deletion residues. Candidate existence does not constrain the selected alignment.

Global Needleman--Wunsch maximizes ``U - lambda*g``, where ``U`` sums nongap
similarities, ``g`` counts gap residues including terminal gaps, and
``lambda = max(abs(similarity)) > 0``. Strict predecessor ties use M/X/Y for the
match and insertion states, M/Y/X for the deletion state, and M/X/Y at the terminal
cell. The normalized feature is ``U - max(H_query, H_reference)``, where ``H`` is
the sum of diagonal self similarities. No integer scaling or Gram transform is used.

Required keyword arrays describe the annotations:

* ``query_target_offsets`` (Q+1) and ``query_target_ids`` (M) are requested-target
  CSR; repeated targets within a query are rejected.
* ``reference_target_offsets`` (R+1) and ``reference_target_ids`` (L) map each
  reference to annotation lanes. ``target_class_ids`` maps each target to a class.
* ``reference_class_offsets`` (R+1) and ``reference_class_ids`` enumerate each
  reference's unique class roles, independently of target annotations.
* ``query_group_ids`` (C by Q), ``reference_group_ids`` (C by L), and
  ``group_offsets`` (C small rectangular matrices) describe additive components.
  Each matrix is indexed by query and reference group ID; ``-1`` is unsupported.
* ``reference_other_counts`` (L) supplies distinct neighboring references with
  another class, and ``reference_multiplicities`` (L) supplies positive observation
  counts. These counts are supplied and checked, never estimated by this reducer.
* ``coefficients`` contains ``[intercept, normalized_coefficient, gap_coefficient]``.

The final predictor is the intercept plus group offsets plus the fitted normalized
and gap terms. Its complementary-log-log transform is multiplied by
``-log10((k+1)/N)``, where ``N`` is the number of unique reference/class roles and
``k`` is the supplied other-class count for that annotation. With
``exclude_exact=True``, every role and annotation of the exact query sequence is
removed; its contribution to an other-class count is subtracted once. This requires
symmetric directional indel caps. Inconsistent counters raise instead of clipping.
These supplied-count semantics do not make the output a calibrated probability.

The result's ``arrays()`` method returns read-only owned zero-copy buffers. Per
requested target: ``scores``, ``raw_similarity``, ``normalized_similarity``, ``pfit``,
``informativeness``, ``gap_counts``, ``best_ref_ids``, ``best_lane_ids``, ``counts``,
``modeled_counts``, ``aa_counts``, and ``exact_counts``. Counts distinguish observation
multiplicity, supported multiplicity, distinct reference sequences, and raw exact
multiplicity before exclusion. IDs are -1 if unavailable. Per query:
``candidate_counts`` counts all full-index candidates before filtering,
``nw_pairs`` counts actual alignments, and ``query_supported`` records group support.
Supported no-hit scores are zero; unsupported scores are NaN while counts remain.
Tied maxima choose the smallest annotation lane ID. Invalid CSR, IDs, duplicate
references or roles, nonfinite coefficients, and overflow raise explicitly.


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

The same question for substring hits in a text, where counts are per distinct record
rather than per offset and the ball is a mismatch radius rather than a penalty. The
module docstring of :mod:`seqtree.text_evalue` explains both choices.

.. autofunction:: text_evalues

.. autofunction:: text_threshold_for_evalue

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
