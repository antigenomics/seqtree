Gap-block alignment and calibrated cutoffs
==========================================

Two questions come up whenever you search a set of V(D)J junctions: *where does the gap go*, and
*how close is close enough*. seqtree answers the first with a restricted alignment plus a prior,
and the second by inverting an E-value. Both answers are measured, not assumed.

.. contents::
   :local:
   :depth: 1

Why one gap block
-----------------

A junction's length varies because of V/J trimming and N-addition — **one** contiguous indel event,
not a scatter of them. So restrict the alignment to a single gap block of length ``d = |m - n|`` and
let it sit anywhere:

.. code-block:: python

   from seqtree import SubstitutionMatrix
   from seqtree.gapblock import gapblock_score

   mat = SubstitutionMatrix.blosum62()
   score, block_pos = gapblock_score("CASSLGQAYEQYF", "CASSGQAYEQYF", mat, gap_open=2 * mat.scale())

``gapblock_score`` enumerates every one of the ``L + 1`` block positions and returns the best, in
``O(min(m, n))`` rather than the ``O(mn)`` of a full DP. It is exact: zero disagreements against
brute-force layout enumeration over 55,727 random pairs.

The restriction costs nothing where it applies. Against a **model-independent** structural oracle
(iterative superposition plus an unrestricted affine DP, free to open any number of blocks) over
3,049 crystal junction pairs from 199 unique sequences, the true residue correspondence is a single
contiguous block in **95.2–100 %** of cases, for both chains and every ``d`` from 1 to 4.

And where it does *not* apply, the restriction protects you. At ``gap_open = 2 * scale``:

=============================  ===========================  ===============
pairs                          gap-block == affine exactly  median excess
=============================  ===========================  ===============
related (1 indel + 0–2 subs)   **98.8 %**                   2
unrelated (independent)        0.0 %                        **106**
=============================  ===========================  ===============

On unrelated sequences unrestricted affine alignment always undercuts the block score — by inventing
an alignment that does not exist. Extra gap freedom buys manufactured similarity.

.. warning::

   ``gap_open`` defaults to **1**, which is wrong for every real matrix. The Gram transform puts a
   typical BLOSUM62 mismatch at ``SubstitutionMatrix.scale() == 14``, so a ``gap_open`` of 1 makes
   gaps fourteen times cheaper than substitutions and every alignment degenerates to gaps. Pass
   ``gap_open = 2 * matrix.scale()``.

Choosing the block position
---------------------------

A sequence score alone cannot place the block. Measured against structure, minimum-BLOSUM62 agrees
with the structurally correct position about as often as picking at random, and a hard central pin
agrees with the score-only choice on only **10.6 %** of pairs. So the score gets a prior.

.. code-block:: python

   from seqtree.gapblock import central_prior, profile_prior, frame_prior

   lam = int(1.5 * mat.scale())            # 21 for BLOSUM62
   prior = central_prior(lam)
   score, pos = gapblock_score(q, r, mat, gap_open=2 * mat.scale(), gap_prior=prior)

A prior is any callable ``prior(block_start, block_length, longer_length) -> int``. It must satisfy
**exactly two** invariants:

1. ``prior >= 0`` — otherwise trie pruning stops being admissible.
2. ``prior == 0`` when ``d == 0`` — otherwise ``s(q, q) != 0``, the score stops defining a ball, and
   the whole E-value construction (:doc:`evalue`) collapses.

Monotonicity in ``d`` is *not* required, and is false for :func:`~seqtree.gapblock.central_prior`:
growing a leading block drags its midpoint toward the centre, so the penalty falls.

Three priors ship:

``central_prior(lam)``
    ``lam * |block_midpoint - m/2|``. The block sits at the loop apex — measured at Cys-offset 6 for
    both TRA and TRB, and it does **not** drift with ``d``.

``profile_prior(lam, w)``
    ``lam * sum(w(j, m) for j in block)``, charging per unit of positional weight the block deletes.

``frame_prior(lam, c)``
    ``lam * |i - c|``, pinning the block to a fixed column. The only transitive rule — see below.

A rule is a column frame
------------------------

Pairwise-optimal gap placement is **not transitive**. Align *A* to *B* and *B* to *C* independently
and the two column assignments do not compose, so a set of unequal-length sequences has no
consistent column index — and therefore no position weight matrix.

A frame rule fixes that, but only if its block start does not depend on ``d``. Under
:func:`~seqtree.gapblock.central_prior` the start drifts, and two shorter members end up related by
*two* blocks:

.. code-block:: text

   len 14 (d=0)   C A S S L G Q A Y E Q Y F F
   len 13 (d=1)   C A S S L G - Q A Y E Q Y F     block start 6
   len 11 (d=3)   C A S S L - - - G Q A Y E Q     block start 5  <- it moved

   rows 2 and 3: the residues at columns 5 and 7 of the longer are both unmatched,
   on opposite sides of the shorter's gap. Two blocks, not one.

Pin the block instead, and embedding reproduces the pairwise alignment exactly:

.. code-block:: python

   >>> from seqtree.gapblock import embed_in_frame
   >>> for s in ("CASSLGQGAYEQYF", "CASSLGQAYEQYF", "CASSGQAYEQYF"):
   ...     print(embed_in_frame(s, 14, 4))
   CASSLGQGAYEQYF
   CASS-LGQAYEQYF
   CASS--GQAYEQYF

Columns ``0..c-1`` are left-anchored on the conserved Cys, columns ``c+d..W-1`` right-anchored on the
Phe. Now every column means the same thing in every member, and a PWM is well defined.

Profiling an island
-------------------

:class:`~seqtree.gapblock.IslandProfile` is that PWM. Its column penalty is measured against the
column's **own consensus**, ``pen(j, a) = round(lam * log(p_max_j / p_j(a)))`` — a textbook log-odds
score is signed, and a signed score is not a ball. This one is ``>= 0`` and zero on the consensus,
so it flows through :func:`~seqtree.thetas_from_scores` unchanged.

.. code-block:: python

   from seqtree.gapblock import IslandProfile

   profile = IslandProfile.fit(island_members)     # c defaults to the entropy-optimal column
   control_scores = profile.score_batch(control_seqs)
   theta = seqtree.thetas_from_scores([control_scores], n_target=N, m_control=M,
                                      e_target=0.05, theta_max=1 << 20)[0]

**Whether it beats scoring against every member depends entirely on how strict your cutoff is**, and
the cutoff moves with ``N``. The E-value's ``k = floor(e_target · M / N)`` is the number of control
neighbours the cutoff may admit, so the false-positive rate is ``k / M``:

* Building islands *within* one epitope group puts ``N`` at the group size (median 88 in VDJdb), so
  ``k`` has median 142 out of ``M = 250,000`` — an FPR of ``5.7e-4``.
* Annotating a whole repertoire against known islands puts ``N`` at ≈ 20,000. Then
  ``e_target = 0.05`` yields ``k = 0``, which :func:`~seqtree.thetas_from_scores` reports as ``-1``:
  the rule of three certifies no ``E`` below ``3N/M = 0.236``. At that smallest certifiable ``E``,
  ``k = 3`` — an FPR of ``1.2e-5``.

Recall on held-out members of 108 calibrated VDJdb islands of ≥ 10 (human TRB, three splits each,
paired bootstrap over islands, 250,000 control junctions as negatives):

===================  ==========  ==================  =================  =======================
regime               FPR         min-over-members    ``IslandProfile``  difference [95% CI]
===================  ==========  ==================  =================  =======================
loose reference      1 %         **99.5 %**          99.1 %             −0.40 [−1.09, +0.14]
per-epitope islands  0.0568 %    88.3 %              **89.3 %**         +0.93 [−0.80, +2.79]
repertoire           0.0012 %    37.6 %              **48.5 %**         +10.90 [+7.69, +14.21]
===================  ==========  ==================  =================  =======================

So there is **no significant difference while you are building the islands**, and a large one when
you use them to annotate a repertoire. On islands of ≥ 50 members the repertoire-regime gap is
9.8 % against 22.6 %.

.. warning::

   A profile does **not** generalise beyond its island. Junctions specific to the same epitope that
   landed in a *different* island are recovered 3.5 % of the time by the profile and 3.7 % by
   min-over-members at a 1 % FPR, and by neither at either operating point. Distinct islands of one
   epitope share no motif that similarity can find. Fit a profile to recognise *this* island's
   members more sharply, not to discover new ones.

   Nor is it a compression. 14 columns × 21 symbols × 4 B is 1,176 B against 182 B of member
   strings; an island needs 84 members before the profile is the smaller of the two, which 3.7 % of
   real islands reach.

The entropy-optimal frame column is modal at ``c = 6`` across real islands — the same place the
crystal structures put the block, arrived at from sequence alone.

Searching a reference set
-------------------------

:class:`~seqtree.gapblock.GapBlockIndex` reuses the ordinary Hamming engine over deletion variants —
no new C++, and no separate aligner.

.. code-block:: python

   from seqtree.gapblock import GapBlockIndex

   gbi = GapBlockIndex(cdr3s, "aa", d_max=2)
   for ref_id, score, block_len, block_pos in gbi.search(
           query, max_penalty=40, matrix=mat,
           gap_open=2 * mat.scale(), gap_prior=central_prior(lam)):
       ...

``d_max`` bounds the length difference. On VDJdb same-epitope pairs, ``d <= 1`` covers 49.2 % of
them and ``d <= 3`` covers 85.5 %, so the default of 1 is deliberately conservative — raise it when
you care about length-different neighbours, and pay for it in build memory.

Scoring every pair
------------------

A search prunes; a **distance-vector embedding cannot**. Scoring *n* clonotypes against a few
thousand fixed prototypes needs all *n × K* cells, so there is nothing for a trie to skip.
:func:`~seqtree.gapblock.score_matrix` is that dense path, in C++ with the GIL released:

.. code-block:: python

   from seqtree.gapblock import score_matrix, central_prior

   sm = score_matrix(clonotypes, prototypes, mat,
                     gap_open=2 * mat.scale(),
                     gap_prior=central_prior(lam),
                     threads=0)                # 0 = one per core
   d = numpy.asarray(sm)                       # (n, K) int32, no copy

Measured on an M3, 3,000 prototypes, human TRB junctions:

=================================  ===========  ===============
rung                               M pairs/s    vs pure Python
=================================  ===========  ===============
``gapblock_score`` (Python)        0.41         1×
``score_matrix``, 1 thread         **51.3**     125×
``score_matrix``, 16 threads       **532.7**    1294×
=================================  ===========  ===============

The prior is free — it is flattened once into an ``[m][d][i]`` lookup cube, so the kernel never
re-enters Python. The result is ``int32`` and carries the CPython buffer protocol, so
``numpy.asarray`` wraps it without copying and seqtree keeps its zero runtime dependencies. Budget
``4 * n * K`` bytes (1.2 GB at 100k × 3000) and chunk the queries if that does not fit.

.. note::

   The scores are already the distance. There is no ``d = s(a,a) + s(b,b) - 2·s(a,b)`` step to do:
   the Gram transform is applied per *residue* when the matrix is built, so the alignment score is
   non-negative, zero on the diagonal, and symmetric by construction.

One knob exists only for interoperability. :func:`~seqtree.gapblock.positions_prior` restricts the
block to a fixed list of starts (negative values counting from the end), reproducing the
``gap_positions=(3, 4, -4, -3)`` convention that other junction aligners hardcode:

.. code-block:: python

   from seqtree.gapblock import positions_prior
   sm = score_matrix(clonotypes, prototypes, mat, gap_open=28,
                     gap_prior=positions_prior((3, 4, -4, -3)))

Reach for it to reproduce someone else's numbers, not to improve your own. On human TRB retrieval at
a matched false-positive rate, candidate starts ``(3, 4, mid)`` reached precision **0.156** against
**0.414** for a single hard-pinned centre: scoring several placements and keeping the best is *more*
freedom, and freedom is what manufactures similarity. A fixed set of four starts is also a TRB-shaped
assumption — an IGH junction runs to 50 residues, and a block pinned within four of an anchor there
means something quite different.

Calibrated cutoffs
------------------

.. important::

   **A fixed score cutoff is not a calibrated cutoff.** A control repertoire is dense near germline
   and sparse among rare junctions, so one threshold buys a common query many more chance neighbours
   than a rare one. Building a neighbour graph on human TRB at ``gapblock_score <= 60``, **31.7%** of
   size-matched *random control* junctions land in a component of ≥ 5. That structure is invented by
   the threshold; per-query cutoffs remove all of it.

:func:`~seqtree.threshold_for_evalue` inverts ``E = (N/M) * n_control`` into the score cutoff that
achieves a target E — **per query**. Because scores are integers the inversion is exact rather than
a root-find: sort a query's control-hit scores and the answer sits just below the ``(k+1)``-th
smallest, with ``k = floor(e_target * M / N)``.

.. code-block:: python

   import seqtree

   control = seqtree.load_control("human_trb_aa")          # 250k bundled
   target = seqtree.Index.build(vdjdb_cdr3s, "aa")

   ceiling = seqtree.SearchParams(max_subs=14, max_penalty=50,
                                  matrix="BLOSUM62", engine="seqtm")
   thetas = seqtree.threshold_for_evalue(target, control, queries, ceiling,
                                         e_target=0.05, exclude_exact=True)

   for q, theta in zip(queries, thetas):
       if theta < 0:
           continue        # e_target < 3N/M: this control cannot certify the bar
       hits = [h for h in target.search(q, ceiling) if h.score <= theta]

One control scan supplies every query's cutoff. A ``-1`` is not a failure to compute — it is the
rule of three telling you the control is too small for the E you asked for. Enlarge it rather than
lowering the bar.

Re-run the same neighbour graph on E-value edges (an edge needs *mutual* significance,
``score <= min(theta_a, theta_b)``) and the picture inverts:

===========  ==============  ================  ================
arm          edges per node  in island ≥ 5     in island ≥ 20
===========  ==============  ================  ================
real         **2.334**       **0.233**         **0.152**
control      0.021           0.000             0.000
===========  ==============  ================  ================

The control forms 19,248 singletons, 223 pairs, three components of 3–4, and **nothing larger**. Its
realised edge rate lands on ``E* = 0.05`` as designed, which is the check that the calibration is
honest.

See also
--------

* :doc:`evalue` — the theory the cutoff inverts.
* :doc:`api` — full signatures for :mod:`seqtree.gapblock`, :mod:`seqtree.seeds` and
  :mod:`seqtree.layout`.
* ``bench/bench_gapblock.py`` — reproduces the gap-freedom ladder table above.
