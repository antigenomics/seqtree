"""Presentation-aware E-values for peptide-MHC epitope homology.

The null is allele-conditional and anchor-masked: significance is computed over the
TCR-facing readout against a background of peptides presented by the same MHC, so
shared anchors (presentation, not recognition) do not inflate hits. Reuses the
Poisson tail of :mod:`seqtree.evalue`; see ``appendix/evalue.tex``.
"""
import math

from .evalue import _poisson_sf


def homolog_evalue(n_target, n_control, n_ref, m_control):
    """E-value of a query's homolog neighbourhood.

    n_target  : homologs found in the searched set (target)
    n_control : homologs found in the per-allele presented background
    n_ref     : size of the target set (N)
    m_control : size of the background (M)
    """
    if m_control <= 0:
        return {"n_target": n_target, "n_control": n_control, "E": float("inf"),
                "p_any": 1.0, "p_enrichment": 1.0, "rule_of_three": False}
    rule3 = n_control == 0
    E = (3.0 if rule3 else float(n_control)) * n_ref / m_control
    return {
        "n_target": n_target,
        "n_control": n_control,
        "E": E,
        "p_any": 1.0 - math.exp(-E) if E < 700 else 1.0,
        "p_enrichment": _poisson_sf(n_target, E),
        "rule_of_three": rule3,
    }
