"""Presentation-aware E-values for peptide-MHC epitope homology.

The null is allele-conditional and anchor-masked: significance is computed over the
TCR-facing readout against a background of peptides presented by the same MHC, so
shared anchors (presentation, not recognition) do not inflate hits. The arithmetic is
:func:`seqtree.evalue.evalue_result`; only the choice of background differs. See
``appendix/evalue.tex``.
"""
from .evalue import evalue_result


def homolog_evalue(n_target, n_control, n_ref, m_control):
    """E-value of a query's homolog neighbourhood against a per-allele background.

    Args:
        n_target: homologs found in the searched (target) set.
        n_control: homologs found in the per-allele presented background.
        n_ref: size of the target set, ``N``.
        m_control: size of the background, ``M``.

    Returns:
        dict with ``n_target``, ``n_control``, ``E``, ``p_any``, ``p_enrichment``,
        ``rule_of_three``.
    """
    return evalue_result(n_target, n_control, n_ref, m_control)
