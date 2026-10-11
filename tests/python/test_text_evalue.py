"""E-values for substring hits.

The C++ suite owns the completeness of the search itself. These tests own the
calibration layer: that a planted motif is significant while background is not, that
counts collapse dependent offsets onto records, that the ball nests, and that a count
which is really a floor is refused rather than reported.
"""
import random

import pytest

import seqtree

ALPHABET = "ACDEFGHIKLMNPQRSTVWY"
MOTIF = "WYWMHCWY"


def _records(n, length, seed):
    rng = random.Random(seed)
    return ["".join(rng.choice(ALPHABET) for _ in range(length)) for _ in range(n)]


@pytest.fixture(scope="module")
def planted():
    """A target where 40 of 400 records carry a motif, and a matched control."""
    target = _records(400, 120, seed=1)
    for i in range(0, 40):
        target[i] = target[i][:50] + MOTIF + target[i][50 + len(MOTIF):]
    control = _records(400, 120, seed=2)
    return (seqtree.TextIndex.build(target, "aa", k=4),
            seqtree.TextIndex.build(control, "aa", k=4))


def test_planted_motif_is_significant_and_background_is_not(planted):
    target, control = planted
    background = _records(1, 8, seed=99)[0]
    results = seqtree.text_evalues(target, control, [MOTIF, background], max_subs=1)

    planted_result, background_result = results
    assert planted_result["n_target"] >= 40
    assert planted_result["p_enrichment"] < 1e-6
    assert background_result["p_enrichment"] > 0.01


def test_counts_are_records_not_offsets():
    """A motif repeated inside one record counts once, not once per occurrence.

    Overlapping occurrences in a repeat are strongly dependent; counting them
    separately is what makes the Poisson estimator over-disperse.
    """
    repeat = "AAAA" + MOTIF + "AAAA" + MOTIF + "AAAA" + MOTIF + "AAAA"
    target = seqtree.TextIndex.build([repeat], "aa", k=4)
    control = seqtree.TextIndex.build(_records(50, 60, seed=3), "aa", k=4)

    hits = target.search_batch([MOTIF], max_subs=0)
    assert hits.num_hits == 3

    result = seqtree.text_evalues(target, control, [MOTIF], max_subs=0)[0]
    assert result["n_target"] == 1


def test_evalue_is_monotone_in_the_radius(planted):
    target, control = planted
    previous = -1.0
    for radius in (0, 1, 2):
        result = seqtree.text_evalues(target, control, [MOTIF], max_subs=radius)[0]
        assert result["E"] >= previous
        previous = result["E"]


def test_exclude_exact_drops_only_zero_mismatch_records(planted):
    target, control = planted
    with_exact = seqtree.text_evalues(target, control, [MOTIF], max_subs=1)[0]
    without = seqtree.text_evalues(
        target, control, [MOTIF], max_subs=1, exclude_exact=True
    )[0]
    assert without["n_target"] < with_exact["n_target"]


def test_truncated_search_is_refused(planted):
    """A capped count is a floor; an E-value computed from it would be wrong."""
    target, control = planted
    capped = target.search_batch([MOTIF], max_subs=1, max_hits=1)
    assert any(memoryview(capped.arrays()["truncated"]).tolist())
    with pytest.raises(ValueError, match="truncated"):
        capped.record_distances()


def test_empty_control_is_rejected(planted):
    target, _ = planted
    empty = seqtree.TextIndex.build([], "aa", k=4)
    with pytest.raises(ValueError, match="control text is empty"):
        seqtree.text_evalues(target, empty, [MOTIF], max_subs=0)


def test_threshold_round_trips_to_the_requested_evalue(planted):
    target, control = planted
    queries = [MOTIF, _records(1, 10, seed=7)[0]]
    cutoffs = seqtree.text_threshold_for_evalue(
        target, control, queries, max_subs=2, e_target=1.0
    )
    assert len(cutoffs) == len(queries)
    for query, cutoff in zip(queries, cutoffs):
        if cutoff < 0:
            continue
        achieved = seqtree.text_evalues(target, control, [query], max_subs=cutoff)[0]
        assert achieved["E"] <= 1.0 + 1e-9


def test_threshold_is_unreachable_when_the_control_is_too_small(planted):
    """e_target below 3N/M cannot be certified however tight the ball."""
    target, control = planted
    cutoffs = seqtree.text_threshold_for_evalue(
        target, control, [MOTIF], max_subs=2, e_target=1e-9
    )
    assert cutoffs == [-1]


def test_negative_radius_is_rejected(planted):
    target, control = planted
    with pytest.raises(ValueError, match="max_subs must be >= 0"):
        seqtree.text_threshold_for_evalue(
            target, control, [MOTIF], max_subs=-1, e_target=1.0
        )


def test_record_distances_and_exact_record_exclusion():
    # Build-time groups do not change the counting unit; a record containing an
    # exact AND an inexact occurrence is wholly excluded by exclude_exact.
    target = seqtree.TextIndex.build(["AAAAACAA", "AAAC", "AAAC"], k=2, group_ids=[0, 0, 0])
    control = seqtree.TextIndex.build(["CCCC", "CCCC", "CCCC"], k=2)
    hits = target.search_batch(["AAAA"], max_subs=1, group_by=True)
    assert hits.record_distances() == [[0, 1, 1]]
    assert seqtree.text_evalues(target, control, iter(["AAAA"]), max_subs=1)[0]["n_target"] == 3
    assert seqtree.text_evalues(target, control, ["AAAA"], max_subs=1, exclude_exact=True)[0]["n_target"] == 2
    assert target.search_batch([], max_subs=1).record_distances() == []


def test_calibration_requires_matching_alphabets_and_nonempty_control():
    target = seqtree.TextIndex.build(["AAAA"], k=2)
    control = seqtree.TextIndex.build(["AAAA"], alphabet="nt", k=2)
    for call in (seqtree.text_evalues, seqtree.text_threshold_for_evalue):
        opts = dict(max_subs=1)
        if call is seqtree.text_threshold_for_evalue:
            opts["e_target"] = 1.0
        with pytest.raises(ValueError, match="alphabets"):
            call(target, control, ["AAAA"], **opts)
    with pytest.raises(ValueError, match="control text is empty"):
        seqtree.text_threshold_for_evalue(target, seqtree.TextIndex.build([], k=2), ["AAAA"], 1, 1.0)
