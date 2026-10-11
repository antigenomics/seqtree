"""Control-calibrated E-values for substring hits in a text.

:func:`seqtree.evalues` answers "is this whole-sequence neighbourhood bigger than
background?" for an :class:`~seqtree.Index`. This module asks the same question of a
:class:`~seqtree.TextIndex`, where a short query is matched *inside* long records --
a peptide against a proteome, a read against a transcript set. The arithmetic is
:func:`seqtree.evalue.evalue_result`; only the ball and the counting unit differ.

    E(q)          = (N / M) * n_control(q)
    p_any(q)      = 1 - exp(-E)
    p_enrichment  = P(Poisson(E) >= n_target(q))

Three things are specific to text search, and each is a decision rather than a detail.

**The ball is ``max_subs``, never ``score``.** ``TextResult.score`` is a similarity --
``0`` for an exact match, more negative for a worse one -- the opposite sign to the
non-negative penalty :func:`seqtree.thetas_from_scores` expects. ``n_subs`` is a
non-negative integer, zero on an exact match, and monotone in the scope, so it is the
quantity that satisfies the nesting property the derivation rests on.

**Counts are per distinct record, not per offset.** One query usually hits the same
record at several offsets -- a repeat, a low-complexity stretch, a tandem domain -- and
those events are strongly dependent. Counting them separately gives an over-dispersed,
compound-Poisson count and an E-value that is too large; collapsing to distinct records
avoids counting dependent offsets repeatedly. The Poisson approximation still
requires independent records and rare background hits. This is the same reason
``appendix/evalue.tex`` counts unique clonotypes rather than reads, and the same reason
it says to count the union of a dependent family rather than correct for its size.
``N`` and ``M`` are therefore ``num_refs``, the number of records in each text.

**The control has to be a matched text.** The estimator assumes target and control are
drawn from one background: same residue composition, same record-length distribution,
same low-complexity content. A shuffled or reversed copy of the target text, or an
unrelated proteome of similar composition, is the usual choice. There is no
``load_control`` for texts -- that choice is the caller's, because it is the assumption
the whole calibration rests on, and nothing here can check it.
"""
from .evalue import evalue_result, thetas_from_scores

__all__ = ["text_evalues", "text_threshold_for_evalue"]


def _search(index, queries, max_subs, threads):
    """One native batch search and record reduction; refuse incomplete counts."""
    return index.search_batch(queries, max_subs=max_subs, threads=threads).record_distances()


def _counts(distances, exclude_exact):
    """Distinct records per query, optionally dropping exact-match records."""
    return [sum(d > 0 for d in row) if exclude_exact else len(row) for row in distances]


def text_evalues(target, control, queries, max_subs=0, threads=0, exclude_exact=False):
    """Control-calibrated E-values for each query's substring neighbourhood.

    Args:
        target: :class:`~seqtree.TextIndex` to score hits in, e.g. a proteome.
        control: background :class:`~seqtree.TextIndex`, matched to ``target`` in
            composition and record-length distribution. See the module note -- this
            choice is the assumption the calibration rests on.
        queries: query strings, each at least ``target.k`` long.
        max_subs: substitution radius defining the ball, searched identically in both
            texts.
        threads: worker threads for the batch searches (0 = all cores).
        exclude_exact: drop 0-mismatch hits from both counts. Set this when a query may
            occur verbatim in either text. A record with any exact occurrence is
            excluded, even if it also contains inexact occurrences.

    Returns:
        One dict per query with ``n_target``, ``n_control``, ``E``, ``p_any``,
        ``p_enrichment`` and ``rule_of_three``. Counts are distinct records, not offsets.

    Raises:
        ValueError: If the control text is empty, or if a search was truncated, which
            would make a count a floor rather than a count.

    Example:
        >>> import seqtree
        >>> proteome = seqtree.TextIndex.build(["MKVLAAGIVGLNLGGK"], "aa", k=4)
        >>> shuffled = seqtree.TextIndex.build(["GKLNVGALMGIVKGLA"], "aa", k=4)
        >>> seqtree.text_evalues(proteome, shuffled, ["AAGIVGL"])[0]["n_target"]
        1
    """
    queries = list(queries)
    if target.alphabet != control.alphabet:
        raise ValueError("target and control alphabets must match")
    n_ref, m_control = target.num_refs, control.num_refs
    if m_control == 0:
        raise ValueError("control text is empty")

    target_hits = _counts(
        _search(target, queries, max_subs, threads), exclude_exact
    )
    control_hits = _counts(
        _search(control, queries, max_subs, threads), exclude_exact
    )
    return [
        evalue_result(n_t, n_c, n_ref, m_control)
        for n_t, n_c in zip(target_hits, control_hits)
    ]


def text_threshold_for_evalue(target, control, queries, max_subs, e_target, threads=0,
                              exclude_exact=False):
    """Per-query mismatch cutoff achieving ``e_target`` against ``control``.

    A fixed radius is not calibrated: at ``max_subs=2`` a low-complexity query picks up
    thousands of chance occurrences while a compositionally unusual one picks up none,
    so the radius that buys a given E-value has to be found per query. One control
    search at ``max_subs`` supplies every cutoff.

    Args:
        target: :class:`~seqtree.TextIndex` the E-value is expressed against. Only
            ``num_refs`` is used.
        control: background :class:`~seqtree.TextIndex`.
        queries: query strings.
        max_subs: radius the control is searched at, and the ceiling for the returned
            cutoffs.
        e_target: desired E-value, e.g. ``0.05``.
        threads: worker threads for the batch search (0 = all cores).
        exclude_exact: drop 0-mismatch hits, matching :func:`text_evalues`.

    Returns:
        One integer cutoff per query in ``[0, max_subs]``, or ``-1`` where no radius
        achieves ``e_target``. That happens when ``e_target < 3N/M``: with only ``M``
        control records, even an empty ball certifies no better than the rule-of-three
        bound. Enlarge the control text.

    Raises:
        ValueError: If ``max_subs`` is negative, or a search was truncated.
    """
    if max_subs < 0:
        raise ValueError("max_subs must be >= 0")
    if target.alphabet != control.alphabet:
        raise ValueError("target and control alphabets must match")
    if control.num_refs == 0:
        raise ValueError("control text is empty")
    distances = _search(control, list(queries), max_subs, threads)
    return thetas_from_scores(
        distances, target.num_refs, control.num_refs, e_target, max_subs,
        exclude_exact=exclude_exact,
    )
