Examples
========

Four runnable scripts in ``examples/``. Each is self-contained: the only data is the 250k control
repertoire bundled with the package, and the only import is ``seqtree`` plus the standard library.

.. code-block:: fish

   python examples/01_gapped_search.py
   python examples/02_sequence_dendrogram.py
   python examples/03_indel_positions.py
   python examples/04_island_profile.py

.. contents::
   :local:
   :depth: 1

1. Gapped search, with a cutoff you can defend
----------------------------------------------

``examples/01_gapped_search.py`` plants twelve relatives of ``CASSLGQAYEQYF`` in a slice of the
control -- some substituted, some a residue shorter, some a residue longer -- and looks for them.

A Hamming ball can only ever reach references of the query's own length:

.. code-block:: text

   Hamming ball (max_penalty 40): 10 hits, 0 of them length-different
   gap-block ball  (max_penalty 60): 44 hits, 12 of them length-different

An indel costs ``gap_open = 28``, so gapped hits never top a score-sorted list; the script surfaces
them separately. Then it derives the cutoff instead of choosing one:

.. code-block:: text

   N = 2,012 references, M = 250,000 control. Rule of three: no E* below 3N/M = 0.024.
     E* = 1.0    theta per query:  CASSLG..=32  CSARDI..=50
     E* = 0.05   theta per query:  CASSLG..=16  CSARDI..=50
     E* = 0.001  theta per query:  CASSLG..=unreachable  CSARDI..=unreachable

        E*  theta   hits  planted  background  precision
       1.0     32      5        4           1       0.80
      0.05     16      2        2           0       1.00

Two things to take from this. The two queries get **different** cutoffs: ``CASSLGQAYEQYF`` sits in a
dense, near-germline neighbourhood where the same score buys far more chance neighbours, so it is
held to a stricter bar (32 against 50 at ``E* = 1.0``). And ``E* = 0.001`` comes back ``-1`` -- with
250k control sequences no cutoff can certify a bar below ``3N/M``. That is the rule of three, not a
bug.

2. A dendrogram over an island
-------------------------------

``examples/02_sequence_dendrogram.py`` builds a family with internal structure -- a substitution
branch and a branch defined by one deletion -- and clusters it by average linkage on the gap-block
distance. Thirty lines of stdlib; seqtree has no runtime dependencies and the example keeps it that
way. Output is Newick plus an ASCII rendering.

.. code-block:: text

   merge                                           height
   CASSLGAYEQYF        + CASSVGAYEQYF               6.0  -----
   CASSWGQAYEQYF       + CASSWGQSYEQYF              6.0  -----
   ...
   (...)               + (...)                     52.3  ------------------------------------

The root sits at the gap cost, 28 units, because that is what separates the deletion branch from the
rest. Then the tree is cut where the control says to cut it:

.. code-block:: text

   N = 20,012 (the annotation set), M = 250,000 (the control)
     E* = 1.0   mean theta 24.4   8/11 merges survive the cut   -> 4 clusters
     E* = 0.1   mean theta 15.1   6/11 merges survive the cut   -> 6 clusters

.. warning::

   Do not build one of these over a whole repertoire: a dense condensed distance matrix over the
   ~90k CDR3 of a VDJdb shortlist is **32 GB**. Cluster first, draw trees inside islands.

   And do not read anything into merges *above* the calibrated cutoff. Distinct sequence islands for
   one epitope share no motif that similarity can find -- only ~0.5% of cross-island co-specific
   pairs share even a rare central 4-mer, and none share a 6-mer. Above ``theta`` the tree is
   arithmetic, not biology.

3. Where the indel goes
------------------------

``examples/03_indel_positions.py`` plants a deletion at a known position, adds ``k`` substitutions,
and asks each rule to find it again.

.. code-block:: text

   === deletion planted in the CORE (positions 4 .. L-5): exact recovery ===
     rule                         k=0 subs   k=1 subs   k=2 subs   k=3 subs
     fixed centre                    27.3%      25.0%      23.4%      24.4%
     central prior lam=21            28.3%      25.8%      24.0%      25.0%
     flat (score alone)             100.0%      92.4%      81.8%      76.5%

This is easy to misread. When two sequences **really are** one deletion plus a little noise, the
score finds the gap. So why constrain it at all? Because that ground truth is a *construction*.
Against TCR-pMHC crystal structures -- where "correct" means "the residues that actually superpose"
-- the score agrees about a tenth of the time and a central prior 30–42%. In retrieval at a matched
false-positive rate, letting the score choose cuts precision from **0.414** to **0.176**.

The prior is an assumption, and the script prices it. Move the deletion to the anchors, where V(D)J
recombination does not put it, and the prior is confidently wrong while the score is unbothered:

.. code-block:: text

   === deletion planted at the EDGE (within 2 of an anchor), k = 1 sub ===
     rule                          exact   within 1
     fixed centre                   0.0%       0.0%
     central prior lam=21           0.0%       0.0%
     flat (score alone)            89.7%      98.6%

Finally, the marginal: where does each rule *land* on ordinary, unrelated pairs?

.. code-block:: text

   rule                        i/L mean     sd   histogram of i/L (0 = Cys end, 1 = Phe end)
   fixed centre                   0.496  0.026   ....##....
   central prior lam=21           0.496  0.027   ....##....
   flat (score alone)             0.515  0.260   #+########

Left free, the score puts the block within two residues of an anchor on **25.9%** of unrelated
pairs. That is where a spurious low score hides: shift the whole sequence by one and the conserved
``CASS`` / ``EQYF`` frame still lines up. It is the same mechanism by which unrestricted affine
alignment undercuts the block score by a median of 106 penalty units on unrelated pairs. It does not
find a better alignment. It manufactures one.

4. An island profile, and the cutoff that decides whether to build one
----------------------------------------------------------------------

``examples/04_island_profile.py`` fits an :class:`~seqtree.gapblock.IslandProfile` to a planted
island and scores twenty held-out members against all 250,000 control junctions — twice, once per
scorer. The point is not which wins on one island. It is that *the E-value decides which question
you are asking*:

.. code-block:: text

   regime                         N     k        FPR   note
   per-epitope island            90   138    0.0552%   E*=0.05
   repertoire annotation     20,000     3    0.0012%   E*=0.05 unreachable (3N/M=0.240); using E*=0.240

``k = floor(e_target · M / N)`` is how many control neighbours the calibrated cutoff may admit, so
the false-positive rate is ``k / M`` — and it moves by a factor of fifty depending on whether you are
building islands inside one epitope group or annotating a whole repertoire against them. At the first
cutoff the two scorers are indistinguishable; at the second, over 108 real VDJdb islands, the profile
recovers 48.5 % of held-out members against 37.6 %.

The script says so out loud rather than reporting its own twenty-member table as a result. Note also
that ``E* = 0.05`` is simply *unreachable* at ``N = 20,000``: the rule of three certifies nothing
below ``3N/M``, and :func:`~seqtree.thetas_from_scores` returns ``-1`` instead of a cutoff it cannot
defend.

It also demonstrates :func:`~seqtree.gapblock.score_matrix`: min-over-members is a row minimum of the
250,000 × 40 gap-block matrix, which the C++ kernel produces in milliseconds.

See also
--------

* :doc:`gapblock` — the model behind all four.
* :doc:`evalue` — where the cutoff comes from.
