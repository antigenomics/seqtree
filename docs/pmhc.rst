Epitope (pMHC) search
=====================

TCR cross-reactivity is driven by a shared **central, TCR-facing motif**, not the MHC anchors
(Dolton et al., Cell 2023: one HLA-A\*02:01 TCR recognizes EAAGIGILTV / LLLGIGILVL / NLSALGIFST via
``x-x-x-A/G-I/L-G-I-x-x-x``). ``seqtree.pmhc`` therefore **masks anchor positions** and searches
anchor-masked TCR-facing **k-mers** over the C++ :class:`~seqtree.KmerIndex` (seed-and-gather).

- **MHC class I** (global, length 8–11): anchors at P2 and the C-terminus (PΩ) are masked; the
  central bulge + P1 drive homology across lengths.
- **MHC class II** (local, trim/shift): the 9-mer core register is handled register-agnostically —
  shared central k-mers match regardless of flanks or register.

Anchor positions are **parametrized** (presets per class, overridable) and masking is realized by a
zero-cost :class:`~seqtree.PositionalMatrix` (weight 0 = masked/free; up-weight a hotspot for a
graded TCR-facing distance).

Homology search
---------------

.. code-block:: python

   from seqtree.pmhc import PMHCStore

   store = PMHCStore.from_pmhc("pmhc_full.tsv.gz")        # IEDB-style epitope table
   # or PMHCStore.from_records([{ "epitope": ..., "mhc": ..., "mhc_class": "MHCI"}, ...])

   for h in store.search_homologs("EAAGIGILTV", "mhc1", mhc="HLA-A*02:01", max_subs=2):
       print(h.epitope, h.mhc, h.shared_kmers, h.score)

``mhc=`` restricts to a single allele (cross-allele groove similarity is out of scope). The Dolton
trio is a built-in positive control: each epitope recovers the other two as mutual homologs.

Neoantigen mimic discovery + E-values
-------------------------------------

Significance is **presentation-aware** — calibrated per allele against an anchor-masked presented
background, so anchor sharing (presentation, not recognition) does not inflate hits:

.. code-block:: python

   from seqtree.pmhc import find_mimics

   res = find_mimics(neoantigen, self_set=human_presented, bacterial_sets={"ecoli": ecoli_presented})
   for source, r in res.items():
       print(source, "E=%.3g" % r["E"], "p=%.3g" % r["p_enrichment"],
             [h.epitope for h in r["hits"][:5]])

``self_set`` / ``bacterial_sets`` are presented peptides for the relevant allele (from ``pmhc_data``
or any caller-supplied proteome slice). The math — per-allele null, the ``E=(N/M)·n_control``
estimator, the cross-reactivity distance, and the limitation that distinct alleles are distinct
nulls — is in ``appendix/evalue.tex`` §"Epitopes".

Reverse problem: peptide → allele
---------------------------------

Flipping the weights to score the **anchors** ranks the likely presenting alleles (a lightweight
presentation prior, not a trained predictor):

.. code-block:: python

   store.assign_allele("KLEEEEEEV", "mhc1")   # -> [(allele, score, n_match, n_allele), ...]
