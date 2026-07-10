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

:func:`threshold_for_evalue` inverts this: it turns a target E into the score cutoff that
achieves it. A *fixed* score cutoff is not calibrated -- the control is far denser near
germline than in the rare-junction region, so the same cutoff buys many more chance
neighbours for a common query than for a rare one. The cutoff has to be per query.
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


def _counts(hitlists, exclude_exact):
    """Per-query hit counts, optionally dropping exact (distance-0) self/identity hits."""
    if exclude_exact:
        return [sum(1 for h in hl if h.score > 0) for hl in hitlists]
    return [len(hl) for hl in hitlists]


def evalues(target, control, queries, params, threads=0, exclude_exact=False):
    """Compute control-calibrated E-values for each query.

    Args:
        target: :class:`Index` to score hits in (e.g. VDJdb), unique clonotypes.
        control: background :class:`Index` (e.g. healthy-donor control).
        queries: list of query strings.
        params: :class:`SearchParams` defining the scope/budget (the ball).
        threads: worker threads for the batch searches (0 = all cores).
        exclude_exact: drop distance-0 (exact / self) hits from both target and control
            counts. Set this when queries may themselves be members of the target or
            control (e.g. a VDJdb-vs-VDJdb scan) so the trivial self-match is not counted.

    Returns:
        One dict per query with ``n_target, n_control, E, p_any, p_enrichment,
        rule_of_three``.
    """
    N, M = len(target), len(control)
    if M == 0:
        raise ValueError("control index is empty")
    if M < N:
        warnings.warn(f"control ({M}) smaller than target ({N}); E-values may be imprecise")

    n_t = _counts(target.search_batch(queries, params, threads), exclude_exact)
    n_c = _counts(control.search_batch(queries, params, threads), exclude_exact)

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


def _theta(scores, k, c_max, theta_max):
    """Largest cutoff in ``[0, theta_max]`` whose control count stays within ``k``."""
    if k < 1:
        # One control hit already overshoots, and an empty ball still reports the
        # rule-of-three bound 3N/M, which exceeds e_target whenever c_max < 3.
        return -1
    theta = theta_max if len(scores) <= k else min(theta_max, scores[k] - 1)
    if theta < 0:
        return -1
    n_control = sum(1 for s in scores if s <= theta)
    if n_control == 0 and c_max < 3.0:
        return -1
    return theta


def thetas_from_scores(control_scores, n_target, m_control, e_target, theta_max,
                       *, exclude_exact=False):
    """Invert ``E = (N/M) * n_control`` for the score cutoff, one cutoff per query.

    Scores are integers, so the inversion is exact rather than a root-find: sort a query's
    control-hit scores and the answer is the value just below the ``(k+1)``-th smallest,
    where ``k = floor(e_target * M / N)`` is the largest control count the target E allows.

    Args:
        control_scores: Per query, the scores of its control hits found at ``theta_max``.
            Hits above ``theta_max`` are irrelevant and may be omitted.
        n_target: ``N``, size of the target index.
        m_control: ``M``, size of the control index.
        e_target: Desired E-value, e.g. ``0.05``.
        theta_max: Score ceiling the control was searched at. Returned cutoffs never exceed it.
        exclude_exact: Drop distance-0 hits, matching :func:`evalues`.

    Returns:
        One integer cutoff per query, or ``-1`` where no cutoff achieves ``e_target``. That
        happens when ``e_target < 3N/M``: with only ``M`` control sequences, even an empty
        ball certifies no better than the rule-of-three bound. Enlarge the control.

    Example:
        >>> # E(4) = (100/100)*1 = 1.0 is allowed; E(5) = 2.0 is not.
        >>> thetas_from_scores([[2, 5, 9]], n_target=100, m_control=100, e_target=1.0,
        ...                    theta_max=10)
        [4]
    """
    if n_target <= 0 or m_control <= 0:
        raise ValueError("n_target and m_control must be > 0")
    if e_target <= 0.0:
        raise ValueError("e_target must be > 0")
    if theta_max < 0:
        raise ValueError("theta_max must be >= 0")

    c_max = e_target * m_control / n_target
    k = int(c_max)
    out = []
    for scores in control_scores:
        s = sorted(x for x in scores if not (exclude_exact and x == 0))
        out.append(_theta(s, k, c_max, theta_max))
    return out


def threshold_for_evalue(target, control, queries, params, e_target, threads=0,
                         exclude_exact=False):
    """Per-query score cutoff achieving ``e_target`` against ``control``.

    One control search at ``params.max_penalty`` supplies every cutoff, so calling this and
    then filtering a target search at the same ceiling costs two scans, not two per query.

    Args:
        target: :class:`Index` the E-value is expressed against (only its size is used).
        control: background :class:`Index`.
        queries: list of query strings.
        params: :class:`SearchParams`; ``max_penalty`` is the ceiling ``theta_max`` and must
            be positive.
        threads: worker threads for the batch search (0 = all cores).
        e_target: Desired E-value.
        exclude_exact: Drop distance-0 hits, matching :func:`evalues`.

    Returns:
        One integer cutoff per query; ``-1`` where ``e_target`` is unreachable at this
        control size.

    Raises:
        ValueError: If ``params.max_penalty`` is not positive.
    """
    theta_max = int(params.max_penalty)
    if theta_max <= 0:
        raise ValueError(
            "threshold_for_evalue needs params.max_penalty > 0 as the score ceiling to "
            "search the control at; cutoffs are then found within [0, max_penalty]"
        )
    scores = [[h.score for h in hl] for hl in control.search_batch(queries, params, threads)]
    return thetas_from_scores(scores, len(target), len(control), e_target, theta_max,
                              exclude_exact=exclude_exact)
