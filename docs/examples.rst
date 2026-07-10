Examples
========

Three runnable scripts in ``examples/``. Each is self-contained: the only data is the 250k control
repertoire bundled with the package, and the only import is ``seqtree`` plus the standard library.

.. code-block:: fish

   python examples/01_gapped_search.py
   python examples/02_sequence_dendrogram.py
   python examples/03_indel_positions.py

.. contents::
   :local:
   :depth: 1

1. Gapped search, with a cutoff you can defend
----------------------------------------------

``examples/01_gapped_search.py`` plants twelve relatives of ``CASSLGQAYEQYF`` in a slice of the
control -- some substituted, some a residue shorter, some a residue longer -- and looks for them.

A Hamming ball can only ever reach references of the query's own length:

.. code-block:: text

   Hamming ball (max_penalty 40): 32 hits, 0 of them length-different
   gap-block ball  (max_penalty 60): 170 hits, 75 of them length-different

An indel costs ``gap_open = 28``, so gapped hits never top a score-sorted list; the script surfaces
them separately. Then it derives the cutoff instead of choosing one:

.. code-block:: text

   N = 2,012 references, M = 250,000 control. Rule of three: no E* below 3N/M = 0.024.
     E* = 1.0    theta per query:  CASSLG..=17  CASSLS..=26
     E* = 0.05   theta per query:  CASSLG..=7   CASSLS..=12
     E* = 0.001  theta per query:  CASSLG..=unreachable  CASSLS..=unreachable

        E*  theta   hits  planted  background  precision
       1.0     17      8        5           3       0.62
      0.05      7      2        2           0       1.00

Two things to take from this. The two queries get **different** cutoffs: ``CASSLGQAYEQYF`` sits in a
dense, near-germline neighbourhood where the same score buys far more chance neighbours, so it is
held to a stricter bar. And ``E* = 0.001`` comes back ``-1`` -- with 250k control sequences no
cutoff can certify a bar below ``3N/M``. That is the rule of three, not a bug.

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
     E* = 1.0   mean theta 11.2   4/11 merges survive the cut   -> 8 clusters
     E* = 0.1   mean theta  6.5   2/11 merges survive the cut   -> 10 clusters

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
     fixed centre                    25.5%      27.6%      26.7%      25.4%
     central prior lam=21            25.9%      28.1%      27.1%      25.8%
     flat (score alone)             100.0%      91.2%      83.7%      76.2%

This is easy to misread. When two sequences **really are** one deletion plus a little noise, the
score finds the gap. So why constrain it at all? Because that ground truth is a *construction*.
Against TCR-pMHC crystal structures -- where "correct" means "the residues that actually superpose"
-- the score agrees about a tenth of the time and a central prior 30–42%. In retrieval at a matched
false-positive rate, letting the score choose **halves** precision.

The prior is an assumption, and the script prices it. Move the deletion to the anchors, where V(D)J
recombination does not put it, and the prior is confidently wrong while the score is unbothered:

.. code-block:: text

   === deletion planted at the EDGE (within 2 of an anchor), k = 1 sub ===
     rule                          exact   within 1
     fixed centre                   0.0%       0.0%
     central prior lam=21           0.0%       0.0%
     flat (score alone)            90.5%      98.7%

Finally, the marginal: where does each rule *land* on ordinary, unrelated pairs?

.. code-block:: text

   rule                        i/L mean     sd   histogram of i/L (0 = Cys end, 1 = Phe end)
   fixed centre                   0.496  0.028   ....##....
   central prior lam=21           0.495  0.028   ....##....
   flat (score alone)             0.524  0.246   ++######+#

Left free, the score puts the block within two residues of an anchor on **22.4%** of unrelated
pairs. That is where a spurious low score hides: shift the whole sequence by one and the conserved
``CASS`` / ``EQYF`` frame still lines up. It is the same mechanism by which unrestricted affine
alignment undercuts the block score by a median of 106 penalty units on unrelated pairs. It does not
find a better alignment. It manufactures one.

See also
--------

* :doc:`gapblock` — the model behind all three.
* :doc:`evalue` — where the cutoff comes from.
