import random

import pytest

import seqtree
from seqtree.evalue import _poisson_sf, evalues, thetas_from_scores, threshold_for_evalue

AA = "ACDEFGHIKLMNPQRSTVWY"


def _rand(n, rng, length=14):
    return ["".join(rng.choice(AA) for _ in range(length)) for _ in range(n)]


def _mutate(seq, k, rng):
    s = list(seq)
    for pos in rng.sample(range(len(seq)), k):
        s[pos] = rng.choice(AA)
    return "".join(s)


def _clustered(seeds, per_seed, rng, max_subs=4):
    """A control with the density structure of a real repertoire: dense around germline."""
    return [_mutate(s, rng.randint(1, max_subs), rng) for s in seeds for _ in range(per_seed)]


def _hamming_params(theta, qlen):
    """Unit cost, so a score is a Hamming distance and cutoffs are small integers."""
    return seqtree.SearchParams(max_subs=qlen, max_penalty=theta, engine="seqtm")


def test_poisson_sf_basics():
    assert _poisson_sf(0, 5.0) == 1.0
    assert _poisson_sf(1, 0.0) == 0.0
    # P(Poisson(1) >= 1) = 1 - e^-1
    assert abs(_poisson_sf(1, 1.0) - (1 - 2.718281828 ** -1)) < 1e-6


def test_planted_cluster_is_significant_background_is_not():
    rng = random.Random(0)
    control = seqtree.Index.build(_rand(20000, rng), alphabet="aa")

    q = "".join(rng.choice(AA) for _ in range(14))
    cluster = []
    for _ in range(50):
        s = list(q)
        s[rng.randrange(14)] = rng.choice(AA)
        cluster.append("".join(s))
    target = seqtree.Index.build(_rand(20000, rng) + cluster, alphabet="aa")

    p = seqtree.SearchParams(max_subs=1, engine="seqtm")
    r = evalues(target, control, [q], p)[0]
    assert r["n_target"] >= 40           # the planted neighbors are found
    assert r["p_enrichment"] < 1e-3      # far more neighbors than the background predicts

    qb = "".join(rng.choice(AA) for _ in range(14))  # a random background query
    rb = evalues(target, control, [qb], p)[0]
    assert rb["p_enrichment"] > 0.01     # not significant


def test_exclude_exact_drops_self_hit():
    rng = random.Random(7)
    q = "CASSPGTEAFF"
    mut = []
    for _ in range(2):
        s = list(q)
        s[rng.randrange(len(q))] = rng.choice(AA)
        mut.append("".join(s))
    target = seqtree.Index.build([q] + mut, alphabet="aa")  # q is a member of the target
    control = seqtree.Index.build(_rand(5000, rng), alphabet="aa")
    p = seqtree.SearchParams(max_subs=1, engine="seqtm")

    full = evalues(target, control, [q], p, exclude_exact=False)[0]
    punc = evalues(target, control, [q], p, exclude_exact=True)[0]
    assert punc["n_target"] == full["n_target"] - 1  # the exact self-match is removed


def test_evalue_monotone_in_scope():
    rng = random.Random(1)
    db = _rand(10000, rng)
    target = seqtree.Index.build(db, alphabet="aa")
    control = seqtree.Index.build(_rand(10000, rng), alphabet="aa")
    q = [db[0]]
    e1 = evalues(target, control, q, seqtree.SearchParams(max_subs=1, engine="seqtm"))[0]["E"]
    e2 = evalues(target, control, q, seqtree.SearchParams(max_subs=3, engine="seqtm"))[0]["E"]
    assert e2 >= e1


# ---------------------------------------------------------------- threshold inversion

def test_thetas_from_scores_is_the_exact_inverse():
    # c_max = e_target * M / N = 1.0 * 100/100 = 1, so at most one control hit is allowed.
    assert thetas_from_scores([[2, 5, 9]], 100, 100, 1.0, 10) == [4]
    # Two allowed: cut just below the third-smallest.
    assert thetas_from_scores([[2, 5, 9]], 100, 100, 2.0, 10) == [8]
    # Ten allowed but only three exist: the ceiling stands.
    assert thetas_from_scores([[2, 5, 9]], 100, 100, 10.0, 10) == [10]


def test_theta_is_monotone_in_e_target():
    scores = [[1, 3, 3, 7, 12, 20]]
    prev = None
    for e in (8.0, 4.0, 2.0, 1.0, 0.5):
        (th,) = thetas_from_scores(scores, 10, 100, e, 30)
        if prev is not None:
            assert th <= prev, f"theta must shrink as e_target shrinks (e={e})"
        prev = th


def test_theta_reports_unreachable_when_the_control_is_too_small():
    # e_target < 3N/M: even an empty ball only certifies the rule-of-three bound.
    assert 3 * 100 / 1000 == 0.3
    assert thetas_from_scores([[]], n_target=100, m_control=1000, e_target=0.2, theta_max=9) == [-1]
    # Loosen past 3N/M and the empty ball becomes admissible.
    assert thetas_from_scores([[]], n_target=100, m_control=1000, e_target=0.3, theta_max=9) == [9]


def test_theta_honours_the_rule_of_three_when_ties_empty_the_ball():
    # Three hits all at score 5, k = 2 -> cut at 4, which leaves the ball EMPTY.
    # E is then 3N/M = 0.6, above e_target = 0.2, so no cutoff works.
    assert thetas_from_scores([[5, 5, 5]], n_target=100, m_control=500, e_target=0.2, theta_max=9) == [-1]


def test_thetas_from_scores_rejects_bad_arguments():
    with pytest.raises(ValueError, match="e_target"):
        thetas_from_scores([[1]], 10, 10, 0.0, 5)
    with pytest.raises(ValueError, match="n_target and m_control"):
        thetas_from_scores([[1]], 0, 10, 1.0, 5)
    with pytest.raises(ValueError, match="theta_max"):
        thetas_from_scores([[1]], 10, 10, 1.0, -1)


def test_theta_is_unreachable_when_even_an_exact_ball_is_too_crowded():
    """k = 1 but the two smallest control scores are both 0: the cutoff would be -1."""
    assert thetas_from_scores([[0, 0, 4]], n_target=100, m_control=100, e_target=1.0,
                              theta_max=9) == [-1]
    # c_max < 1 -> not even one control hit is affordable, and an empty ball reports 3N/M.
    assert thetas_from_scores([[7]], n_target=100, m_control=100, e_target=0.5,
                              theta_max=9) == [-1]


def test_evalues_rejects_an_empty_control():
    rng = random.Random(13)
    target = seqtree.Index.build(_rand(10, rng), "aa")
    empty = seqtree.Index.build([], "aa")
    with pytest.raises(ValueError, match="control index is empty"):
        evalues(target, empty, ["CASSPGTEAFF"], seqtree.SearchParams(max_subs=1))


def test_threshold_for_evalue_needs_a_score_ceiling():
    rng = random.Random(11)
    idx = seqtree.Index.build(_rand(50, rng), "aa")
    with pytest.raises(ValueError, match="max_penalty"):
        threshold_for_evalue(idx, idx, ["CASSPGTEAFF"], seqtree.SearchParams(max_subs=1), 1.0)


def test_threshold_round_trips_through_evalues():
    """The returned theta is the largest cutoff with E <= e_target: theta+1 must overshoot."""
    rng = random.Random(3)
    seeds = _rand(40, rng)
    control = seqtree.Index.build(_clustered(seeds, 100, rng), "aa")   # M = 4000
    target = seqtree.Index.build(_clustered(seeds, 10, rng), "aa")     # N = 400
    queries = [_mutate(s, 1, rng) for s in seeds[:12]]

    theta_max, e_target = 8, 1.0
    thetas = threshold_for_evalue(target, control, queries,
                                  _hamming_params(theta_max, 14), e_target)

    checked = 0
    for q, th in zip(queries, thetas):
        if th < 1:
            continue
        at = evalues(target, control, [q], _hamming_params(th, len(q)))[0]["E"]
        assert at <= e_target, f"E({th}) = {at} overshoots {e_target}"
        if th < theta_max:
            above = evalues(target, control, [q], _hamming_params(th + 1, len(q)))[0]["E"]
            assert above > e_target, f"theta={th} was not maximal: E({th + 1}) = {above}"
        checked += 1
    assert checked >= 5, f"only {checked} queries produced a usable cutoff"


def test_calibrated_cutoff_varies_across_queries():
    """The whole point: a common query needs a tighter cutoff than a rare one."""
    rng = random.Random(5)
    hub = "CASSLGQAYEQYF"
    dense = [_mutate(hub, rng.randint(1, 2), rng) for _ in range(2000)]
    control = seqtree.Index.build(dense + _rand(2000, rng, len(hub)), "aa")
    target = seqtree.Index.build(_rand(500, rng, len(hub)), "aa")
    rare = _rand(1, rng, len(hub))[0]

    th_hub, th_rare = threshold_for_evalue(
        target, control, [hub, rare], _hamming_params(6, len(hub)), e_target=1.0
    )
    assert th_hub < th_rare, f"hub cutoff {th_hub} should be tighter than rare cutoff {th_rare}"
