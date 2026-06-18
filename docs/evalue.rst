E-values for TCR hits
=====================

Fuzzy search tells you *which* references are near a query; an **E-value** tells you whether that
proximity is *surprising*. For TCR repertoires the hard part is biological redundancy: convergent
V(D)J recombination, public clones, and clonal expansion mean a query in a common region of
sequence space has many neighbours for purely generative reasons. A naive BLAST E-value (an
i.i.d.-letter null) would call these wildly significant.

.. admonition:: Full derivation

   This page summarises the method; the complete derivation — theorems, finite-sample bounds, the
   selection Q-factor, multiple-testing correction, and the epitope detection-complexity formalism —
   is in the technical appendix: :download:`A control-calibrated E-value for fuzzy TCR sequence
   search (PDF) <_static/evalue-appendix.pdf>` (LaTeX source in ``appendix/evalue.tex``).

seqtree calibrates against a **background control** instead. With a target index of ``N`` unique
clonotypes (e.g. VDJdb) and a control index of ``M`` unique clonotypes from an unselected
repertoire, for a query ``q`` at a fixed scope/budget:

.. math::

   E(q) = \frac{N}{M}\, n_{\mathrm{control}}(q), \qquad
   p_{\mathrm{any}} = 1 - e^{-E}, \qquad
   p_{\mathrm{enrich}} = \Pr\!\big(\mathrm{Poisson}(E) \ge n_{\mathrm{target}}(q)\big).

Redundancy explained by the background process inflates ``n_control`` and hence ``E``, so such hits
are **not** significant; antigen-driven convergence shows up as ``n_target`` exceeding ``E``. This is
the **TCRNET** approach — counting sequence neighbours against a real-world control repertoire —
put on a rigorous, finite-sample footing; it reduces to the classical Karlin–Altschul E-value when
the background is an i.i.d. product measure and alignments are ungapped.

An empirical control already **is** the post-selection law: thymic and peripheral selection are baked
into a real repertoire, so no extra correction is applied when calibrating against it. Where a
*generative* null is needed instead (rare queries with too few control neighbours), the V(D)J
generation probability :math:`P_\mathrm{gen}` is reshaped by the per-sequence Elhanati selection
factor :math:`Q` (:math:`\langle Q\rangle_{P_\mathrm{gen}} = 1`), giving the post-selection law
:math:`P_0 \propto Q\,P_\mathrm{gen}`. A single global thymic-acceptance fraction only rescales
absolute frequencies and is a crude last-resort fallback. See the appendix remark on the selection
factor for the full treatment.

Usage
-----

.. code-block:: python

   import seqtree

   control = seqtree.load_control("human_trb_aa", size=1_000_000)   # cached after first build
   target = seqtree.Index.build(vdjdb_cdr3s, alphabet="aa")          # unique clonotypes

   p = seqtree.SearchParams(max_subs=1, engine="seqtm")
   for q, r in zip(queries, seqtree.evalues(target, control, queries, p)):
       if r["p_enrichment"] < 1e-3:
           print(q, r["E"], r["n_target"], r["n_control"])

``load_control`` ships a small bundled subset for quick use and downloads larger controls from the
``isalgo/airr_control`` dataset on demand (needs ``huggingface_hub``); both are deduplicated to
unique clonotypes. Meaningful E-values need a control at least as large as the target — see the
precision bound below.

Theory
------

The full derivation — Poisson approximation with an explicit Chen–Stein / Le Cam error bound,
the self-match / punctured-null lemma (benchmark-only exact-hit exclusion), clonotype-collapsing for
over-dispersion, the tf-idf = self-information equivalence, multiple-testing control (Bonferroni and
Benjamini–Hochberg), the control-size requirement :math:`M \gtrsim N/(\rho^2 E^\ast)`, the closest-hit
Gumbel law, **epitope detection complexity** from the degree distribution (worked NLV vs GIL example),
the Karlin–Altschul reduction, and the epitope-presentation limitation — is in the technical appendix
linked at the top of this page (LaTeX source ``appendix/evalue.tex``; rebuild with ``make -C appendix``).
