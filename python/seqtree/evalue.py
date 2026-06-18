"""Control-calibrated, BLAST-style E-values for sequence hits.

For a query ``q`` at a fixed scope/budget, with a target index of ``N`` unique
clonotypes and a background control index of ``M`` unique clonotypes:

    E(q)          = (N / M) * n_control(q)          # per-query Poisson intensity / E-value
    p_any(q)      = 1 - exp(-E)                      # P(>= 1 background hit this close)
    p_enrichment  = P(Poisson(E) >= n_target(q))     # excess neighbors over background

Redundancy explained by background V(D)J biology inflates ``n_control`` and hence
``E``, so such hits are *not* significant; antigen-driven convergence shows up as
``n_target`` exceeding ``E``. See ``appendix/evalue.tex`` for the derivation. Both
indices must be deduplicated to unique clonotypes (``load_control`` does this).
"""
import math
import warnings


def _poisson_sf(k, lam):
    """P(Poisson(lam) >= k) via the complementary CDF (stable for small k)."""
    if k <= 0:
        return 1.0
    if lam <= 0.0:
        return 0.0
    term = math.exp(-lam)
    cdf = term
    for i in range(1, k):
        term *= lam / i
        cdf += term
    return min(1.0, max(0.0, 1.0 - cdf))


def evalues(target, control, queries, params, threads=0):
    """Compute control-calibrated E-values for each query.

    Args:
        target: :class:`Index` to score hits in (e.g. VDJdb), unique clonotypes.
        control: background :class:`Index` (e.g. healthy-donor control).
        queries: list of query strings.
        params: :class:`SearchParams` defining the scope/budget (the ball).
        threads: worker threads for the batch searches (0 = all cores).

    Returns:
        One dict per query with ``n_target, n_control, E, p_any, p_enrichment,
        rule_of_three``.
    """
    N, M = len(target), len(control)
    if M == 0:
        raise ValueError("control index is empty")
    if M < N:
        warnings.warn(f"control ({M}) smaller than target ({N}); E-values may be imprecise")

    n_t = [len(h) for h in target.search_batch(queries, params, threads)]
    n_c = [len(h) for h in control.search_batch(queries, params, threads)]

    out = []
    for nt, nc in zip(n_t, n_c):
        rule3 = nc == 0
        # rule-of-three upper bound on the background rate when the control ball is empty
        E = (3.0 if rule3 else float(nc)) * N / M
        out.append({
            "n_target": nt,
            "n_control": nc,
            "E": E,
            "p_any": 1.0 - math.exp(-E),
            "p_enrichment": _poisson_sf(nt, E),
            "rule_of_three": rule3,
        })
    return out
