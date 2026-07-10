"""Single-gap-block alignment: differential tests against brute force and against Gotoh.

The closed form in gapblock.py is easy to get subtly wrong -- an earlier draft silently
dropped the trailing block position and was out of bounds when the ref was longer than the
query. These tests enumerate every layout independently and compare.
"""
import random

import pytest

import seqtree as st
from seqtree.gapblock import (
    GapBlockIndex,
    central_prior,
    deletion_variants,
    gap_cost,
    gapblock_score,
)

AA = "ACDEFGHIKLMNPQRSTVWY"
INF = 1 << 30
BLOSUM = st.SubstitutionMatrix.blosum62()
PEN = {(a, b): BLOSUM.penalty(a, b) for a in AA for b in AA}

# gap_open < gap_extend and gap_open == 0 both exercise the G(0) guard.
GAP_COSTS = [(1, 1), (4, 1), (28, 1), (0, 3), (2, 5), (5, 5), (0, 0)]


def brute_block(q, r, go, ge, prior=None):
    """Enumerate every single contiguous gap-block layout. The reference implementation.

    The gap prior applies only when there is a block to place (d > 0); with equal lengths
    it must not fire, or s(q, q) would stop being zero.
    """
    m, n = len(q), len(r)
    d = abs(m - n)
    shorter = min(m, n)
    prior = prior if (prior is not None and d > 0) else (lambda i, length: 0)
    best = INF
    for i in range(shorter + 1):
        head = sum(PEN[q[j], r[j]] for j in range(i))
        if m >= n:
            tail = sum(PEN[q[j + d], r[j]] for j in range(i, n))
        else:
            tail = sum(PEN[q[j], r[j + d]] for j in range(i, m))
        best = min(best, head + tail + prior(i, shorter))
    return best + gap_cost(d, go, ge)


def gotoh(q, r, go, ge):
    """Unrestricted affine alignment: any number of gap blocks."""
    m, n = len(q), len(r)
    M = [[INF] * (n + 1) for _ in range(m + 1)]
    X = [[INF] * (n + 1) for _ in range(m + 1)]
    Y = [[INF] * (n + 1) for _ in range(m + 1)]
    M[0][0] = 0
    for i in range(1, m + 1):
        X[i][0] = go + (i - 1) * ge
    for j in range(1, n + 1):
        Y[0][j] = go + (j - 1) * ge
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            M[i][j] = PEN[q[i - 1], r[j - 1]] + min(M[i - 1][j - 1], X[i - 1][j - 1], Y[i - 1][j - 1])
            X[i][j] = min(M[i - 1][j] + go, X[i - 1][j] + ge, Y[i - 1][j] + go)
            Y[i][j] = min(M[i][j - 1] + go, Y[i][j - 1] + ge, X[i][j - 1] + go)
    return min(M[m][n], X[m][n], Y[m][n])


def rand_pair(rng, d_choices=(0, 0, 1, 1, 2, 3, 4)):
    m = rng.randint(6, 22)
    d = rng.choice(d_choices)
    n = m - d if rng.random() < 0.5 else m + d
    n = max(1, n)
    return ("".join(rng.choice(AA) for _ in range(m)),
            "".join(rng.choice(AA) for _ in range(n)))


def test_gap_cost_zero_guard():
    # gap_open + (d-1)*gap_extend would be NEGATIVE at d == 0 when gap_open < gap_extend.
    assert gap_cost(0, 1, 9) == 0
    assert gap_cost(1, 11, 1) == 11
    assert gap_cost(3, 11, 1) == 13
    with pytest.raises(ValueError):
        gap_cost(-1, 1, 1)


def test_deletion_variants_enumerate_every_layout():
    # d == 0 has no block: one identity variant, not len(q)+1 copies of it.
    assert deletion_variants("CAST", 0) == [(0, "CAST")]
    assert deletion_variants("CAST", 1) == [(0, "AST"), (1, "CST"), (2, "CAT"), (3, "CAS")]
    assert len(deletion_variants("CASSLGQAYEQYF", 2)) == 12
    with pytest.raises(ValueError):
        deletion_variants("CAST", 5)


def test_gap_prior_does_not_fire_without_a_gap():
    """s(q, q) == 0 is required by the ball definition; an unconditional prior broke it."""
    prior = central_prior(20)
    for q in ("CASSLGQAYEQYF", "CASSLGQAYEQYFF"):  # odd and even lengths
        assert gapblock_score(q, q, BLOSUM, 28, 1, prior)[0] == 0
    # equal-length, non-identical: still no positional cost
    a, b = "CASSLGQAYEQYF", "CASSIRSSYEQYF"
    assert (gapblock_score(a, b, BLOSUM, 28, 1, prior)[0]
            == gapblock_score(a, b, BLOSUM, 28, 1, None)[0])


def test_identity_scores_zero():
    for q in ("CASSLGQAYEQYF", "A", "GG"):
        assert gapblock_score(q, q, BLOSUM)[0] == 0
        assert gapblock_score(q, q, None)[0] == 0


def test_matches_brute_force_over_every_layout():
    rng = random.Random(7)
    for _ in range(4000):
        q, r = rand_pair(rng)
        go, ge = rng.choice(GAP_COSTS)
        assert gapblock_score(q, r, BLOSUM, go, ge)[0] == brute_block(q, r, go, ge)


def test_matches_brute_force_with_a_gap_prior():
    rng = random.Random(11)
    prior = central_prior(20)
    for _ in range(2000):
        q, r = rand_pair(rng, d_choices=(1, 2, 3))
        go, ge = rng.choice(GAP_COSTS)
        assert gapblock_score(q, r, BLOSUM, go, ge, prior)[0] == brute_block(q, r, go, ge, prior)


def test_symmetric_in_its_arguments():
    rng = random.Random(13)
    for _ in range(1500):
        q, r = rand_pair(rng)
        go, ge = rng.choice(GAP_COSTS)
        assert gapblock_score(q, r, BLOSUM, go, ge)[0] == gapblock_score(r, q, BLOSUM, go, ge)[0]


def test_reported_block_position_reproduces_the_score():
    rng = random.Random(17)
    for _ in range(500):
        q, r = rand_pair(rng, d_choices=(1, 2))
        score, i = gapblock_score(q, r, BLOSUM, 28, 1)
        m, n, d = len(q), len(r), abs(len(q) - len(r))
        head = sum(PEN[q[j], r[j]] for j in range(i))
        tail = (sum(PEN[q[j + d], r[j]] for j in range(i, n)) if m >= n
                else sum(PEN[q[j], r[j + d]] for j in range(i, m)))
        assert head + tail + gap_cost(d, 28, 1) == score


def test_is_a_strict_restriction_of_affine():
    """A single block can never beat unrestricted affine, and must tie on a pure length change."""
    rng = random.Random(19)
    ties = 0
    for _ in range(400):
        q, r = rand_pair(rng, d_choices=(1, 2))
        go, ge = 28, 1
        block = gapblock_score(q, r, BLOSUM, go, ge)[0]
        assert block >= gotoh(q, r, go, ge)

    # A pure length change (ref is the query with one contiguous run removed): the block
    # model is exactly optimal, so it must equal Gotoh every time.
    for _ in range(300):
        q = "".join(rng.choice(AA) for _ in range(rng.randint(10, 18)))
        i, d = rng.randrange(1, len(q) - 3), rng.choice((1, 2))
        r = q[:i] + q[i + d:]
        assert gapblock_score(q, r, BLOSUM, 28, 1)[0] == gotoh(q, r, 28, 1)
        ties += 1
    assert ties == 300


def test_negative_gap_costs_rejected():
    with pytest.raises(ValueError):
        gapblock_score("CAST", "CAT", BLOSUM, gap_open=-1)
    with pytest.raises(ValueError):
        central_prior(-1)


# ---------------------------------------------------------------- GapBlockIndex

REFS = ["CASSLGQAYEQYF", "CASSLGQAYEQYFF", "CASSLGQAYEQY", "CASRQGAWDTQYF",
        "CAWSVSGGGTDTQYF", "CASSIRSSYEQYF", "CSARDRTGNGYTF", "CASSLAPGATNEKLF"]


@pytest.mark.parametrize("gap_prior", [None, central_prior(20)])
def test_index_search_matches_the_scorer_on_every_ref(gap_prior):
    """The Hamming engine over deletion variants must reproduce the closed form exactly."""
    gbi = GapBlockIndex(REFS, "aa", d_max=1)
    budget = 120
    for q in REFS + ["CASSLGQAYEQYW", "CASSLGQYEQYF", "CASSLGQAYEQYFFF"]:
        got = {rid: score for rid, score, _, _ in
               gbi.search(q, budget, BLOSUM, gap_open=28, gap_prior=gap_prior)}
        want = {}
        for rid, r in enumerate(REFS):
            if abs(len(q) - len(r)) > 1:
                continue  # outside d_max
            s, _ = gapblock_score(q, r, BLOSUM, 28, 1, gap_prior)
            if s <= budget:
                want[rid] = s
        assert got == want, f"query {q!r}"


def test_index_finds_refs_both_shorter_and_longer():
    gbi = GapBlockIndex(REFS, "aa", d_max=1)
    hits = {gbi.refs[rid]: (score, d) for rid, score, d, _ in
            gbi.search("CASSLGQAYEQYF", 200, BLOSUM, gap_open=28)}
    assert hits["CASSLGQAYEQYF"] == (0, 0)      # exact
    assert hits["CASSLGQAYEQYFF"][1] == 1       # ref longer by one
    assert hits["CASSLGQAYEQY"][1] == 1         # ref shorter by one


def test_index_zero_budget_returns_only_exact_matches():
    gbi = GapBlockIndex(REFS, "aa", d_max=1)
    hits = gbi.search("CASSLGQAYEQYF", 0, BLOSUM, gap_open=28)
    assert [gbi.refs[rid] for rid, *_ in hits] == ["CASSLGQAYEQYF"]


def test_index_dmax_zero_is_hamming():
    gbi = GapBlockIndex(REFS, "aa", d_max=0)
    hits = gbi.search("CASSLGQAYEQYF", 200, BLOSUM, gap_open=28)
    assert all(len(gbi.refs[rid]) == 13 for rid, *_ in hits)


def test_score_is_a_valid_ball_for_evalues():
    """s >= 0 with s(q, q) == 0 is what appendix/evalue.tex needs; check it holds with a prior."""
    rng = random.Random(23)
    prior = central_prior(20)
    for _ in range(500):
        q, r = rand_pair(rng)
        s, _ = gapblock_score(q, r, BLOSUM, 28, 1, prior)
        assert s >= 0
        assert gapblock_score(q, q, BLOSUM, 28, 1, prior)[0] == 0
