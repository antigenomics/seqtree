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
    embed_in_frame,
    frame_prior,
    gap_cost,
    gapblock_score,
    profile_prior,
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
    shorter, longer = min(m, n), max(m, n)
    prior = prior if (prior is not None and d > 0) else (lambda i, d_, m_: 0)
    best = INF
    for i in range(shorter + 1):
        head = sum(PEN[q[j], r[j]] for j in range(i))
        if m >= n:
            tail = sum(PEN[q[j + d], r[j]] for j in range(i, n))
        else:
            tail = sum(PEN[q[j], r[j + d]] for j in range(i, m))
        best = min(best, head + tail + prior(i, d, longer))
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
    for factory in (central_prior, lambda x: profile_prior(x, lambda j, m: 1.0),
                    lambda x: frame_prior(x, 4)):
        with pytest.raises(ValueError):
            factory(-1)


# ---------------------------------------------------------------- the prior protocol

def test_central_prior_is_bit_identical_after_the_reparametrisation():
    """It used to take (i, shorter). shorter == m - d, so abs(2i - shorter) == abs(2i + d - m)."""
    prior = central_prior(21)
    for m in range(6, 20):
        for d in range(0, 5):
            for i in range(0, m - d + 1):
                old = 0 if d == 0 else int(21 * abs(2 * i - (m - d)) // 2)
                assert prior(i, d, m) == old, (i, d, m)


def test_priors_obey_the_two_required_invariants():
    """Zero at d == 0 (so s(q,q) == 0) and non-negative (so trie pruning stays admissible)."""
    w = [0.77, 0.76, 0.71, 0.52, 0.22, 0.05, 0.07, 0.21, 0.44, 0.73, 0.83, 0.86, 0.93, 0.93]
    priors = [central_prior(21), profile_prior(30, lambda j, m: w[j]), frame_prior(21, 5)]
    for prior in priors:
        for m in range(6, 15):
            for i in range(m + 1):
                assert prior(i, 0, m) == 0
            for d in range(1, min(5, m)):
                for i in range(m - d + 1):
                    assert prior(i, d, m) >= 0


def test_central_prior_is_not_monotone_in_d_and_need_not_be():
    """Growing a leading block drags its midpoint toward the centre, so the penalty falls.
    Pruning holds (i, d) fixed and the ball needs only s >= 0 with s(q,q) == 0, so this is
    sound -- but nothing may assume monotonicity of a general prior."""
    prior = central_prior(21)
    assert prior(0, 1, 6) == 52 and prior(0, 2, 6) == 42

    # profile_prior, being a sum of non-negative weights, IS monotone.
    w = [0.9, 0.8, 0.1, 0.0, 0.1, 0.8, 0.9]
    p = profile_prior(100, lambda j, m: w[j])
    for i in range(4):
        for d in range(1, 7 - i - 1):
            assert p(i, d + 1, 7) >= p(i, d, 7)


def test_profile_prior_sums_the_weights_it_deletes():
    w = [0.9, 0.8, 0.1, 0.0, 0.1, 0.8, 0.9]
    prior = profile_prior(100, lambda j, m: w[j])
    assert prior(2, 3, 7) == int(100 * (0.1 + 0.0 + 0.1))    # cheapest window, in the valley
    assert prior(0, 3, 7) == int(100 * (0.9 + 0.8 + 0.1))    # eats the templated flank
    assert prior(0, 0, 7) == 0                                # empty block, empty sum
    # a flat zero weight must reproduce the no-prior score exactly
    q, r = "CASSLGQAYEQYF", "CASSLGQAYEQY"
    free = profile_prior(1000, lambda j, m: 0.0)
    assert gapblock_score(q, r, BLOSUM, 28, 1, free) == gapblock_score(q, r, BLOSUM, 28, 1, None)


def test_profile_prior_accepts_a_fixed_weight_vector():
    w = [0.0, 1.0, 0.0]
    assert profile_prior(10, w)(1, 1, 3) == 10
    assert profile_prior(10, w)(0, 1, 3) == 0


def test_profile_prior_matches_brute_force():
    rng = random.Random(29)
    w = lambda j, m: abs(j - (m - 1) / 2) / m  # noqa: E731  -- a U-shaped weight
    prior = profile_prior(40, w)
    for _ in range(2000):
        q, r = rand_pair(rng, d_choices=(1, 2, 3))
        go, ge = rng.choice(GAP_COSTS)
        assert gapblock_score(q, r, BLOSUM, go, ge, prior)[0] == brute_block(q, r, go, ge, prior)


def test_frame_prior_pins_the_block_and_keeps_the_ball():
    hard = frame_prior(10**6, 4)
    for q, r in (("CASSLGQAYEQYF", "CASSGQAYEQYF"), ("CASSLGQAYEQYF", "CASSLGQAYEQYFFF")):
        assert gapblock_score(q, r, BLOSUM, 28, 1, hard)[1] == 4
    assert gapblock_score("CASSLGQAYEQYF", "CASSLGQAYEQYF", BLOSUM, 28, 1, hard)[0] == 0


# ---------------------------------------------------------------- the frame

def _induced_blocks(ga, gb):
    """Maximal runs of columns where exactly one of the two gapped strings has a gap.

    One run == the pair is related by a single contiguous gap block. Two or more == the frame
    did not induce a single-block correspondence.
    """
    runs, cur = [], None
    for col, (x, y) in enumerate(zip(ga, gb)):
        state = (x == "-") - (y == "-")   # +1: only a gapped, -1: only b gapped, 0: neither
        if state != 0 and state == cur:
            runs[-1][2] = col             # extend the open run
        elif state != 0:
            runs.append([state, col, col])
        cur = state
    return runs


def test_embed_in_frame_round_trips():
    assert embed_in_frame("CASSLF", 6, 3) == "CASSLF"
    assert embed_in_frame("CASLF", 6, 3) == "CAS-LF"
    for s, w, c in (("CASSLF", 5, 3), ("CASSLF", 8, 9)):
        with pytest.raises(ValueError):
            embed_in_frame(s, w, c)


def test_a_constant_c_frame_is_transitive():
    """Embedding two members reproduces their pairwise single-block alignment. This is the
    property that makes a column index -- and therefore a PWM -- well defined."""
    rng = random.Random(31)
    c = 4
    hard = frame_prior(10**6, c)
    for _ in range(300):
        width = rng.randint(10, 18)
        a = "".join(rng.choice(AA) for _ in range(rng.randint(c + 1, width)))
        b = "".join(rng.choice(AA) for _ in range(rng.randint(c + 1, width)))
        ga, gb = embed_in_frame(a, width, c), embed_in_frame(b, width, c)

        runs = _induced_blocks(ga, gb)
        assert len(runs) <= 1, f"frame induced {len(runs)} blocks for {a!r} vs {b!r}"
        if len(a) != len(b):
            assert len(runs) == 1
            # The block starts at residue index c of the longer sequence, which is frame
            # column c + (the shared gap columns both members already carry).
            shared = min(width - len(a), width - len(b))
            assert runs[0][1] == c + shared
            assert gapblock_score(a, b, BLOSUM, 28, 1, hard)[1] == c


def test_a_length_dependent_frame_is_not_transitive():
    """The central rule puts the block midpoint at the frame centre, so its start drifts with
    d. Two shorter members are then related by TWO blocks, not one -- no consistent columns."""
    width = 14
    # central start = (width - d) // 2
    def central_embed(seq):
        d = width - len(seq)
        return embed_in_frame(seq, width, (width - d) // 2)

    a, b = central_embed("CASSLGQAYEQYF"), central_embed("CASSLGQAYEQ")   # d = 1 and d = 3
    assert a == "CASSLG-QAYEQYF"     # block start 6
    assert b == "CASSL---GQAYEQ"     # block start 5 -- it drifted
    runs = _induced_blocks(a, b)
    assert len(runs) == 2, f"expected the correspondence to split; got {runs}"
    assert [r[1] for r in runs] == [5, 7]   # 'G' and 'Q' of the longer, on opposite sides

    # the constant-c rule on the same two sequences yields exactly one block
    assert len(_induced_blocks(embed_in_frame("CASSLGQAYEQYF", width, 4),
                               embed_in_frame("CASSLGQAYEQ", width, 4))) == 1


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


def test_index_rejects_bad_arguments():
    with pytest.raises(ValueError, match="d_max"):
        GapBlockIndex(REFS, "aa", d_max=-1)
    with pytest.raises(ValueError, match="max_penalty"):
        GapBlockIndex(REFS, "aa").search("CASSLGQAYEQYF", -1, BLOSUM)
    with pytest.raises(ValueError, match="c must be"):
        frame_prior(1, -1)


def test_index_len_and_unit_cost_default():
    gbi = GapBlockIndex(REFS, "aa", d_max=1)
    assert len(gbi) == len(REFS)
    # matrix=None must fall back to unit cost with gap_open 1, not crash on scale()
    hits = gbi.search("CASSLGQAYEQYF", 3, None)
    assert any(gbi.refs[rid] == "CASSLGQAYEQYF" and s == 0 for rid, s, _, _ in hits)


def test_index_skips_refs_shorter_than_the_block():
    """A ref of length 1 has no 2-deletion variant; it must be skipped, not crash."""
    gbi = GapBlockIndex(["A", "CASSLGQAYEQYF"], "aa", d_max=2)
    assert len(gbi) == 2
    hits = gbi.search("CASSLGQAYEQYF", 60, BLOSUM, gap_open=28)
    assert [gbi.refs[rid] for rid, *_ in hits] == ["CASSLGQAYEQYF"]


def test_index_handles_an_empty_auxiliary_level():
    """Every ref shorter than d -> that level's index is empty and must be tolerated.

    The budget has to cover gap_cost(3) or the level is skipped before it is ever consulted.
    """
    gbi = GapBlockIndex(["AC", "AG"], "aa", d_max=3)
    assert gbi._var[3][0] is None
    hits = gbi.search("AC", 5, None)          # unit cost: gap_cost(3, 1, 1) == 3 <= 5
    assert hits[0][1] == 0 and gbi.refs[hits[0][0]] == "AC"


def test_score_is_a_valid_ball_for_evalues():
    """s >= 0 with s(q, q) == 0 is what appendix/evalue.tex needs; check it holds with a prior."""
    rng = random.Random(23)
    prior = central_prior(20)
    for _ in range(500):
        q, r = rand_pair(rng)
        s, _ = gapblock_score(q, r, BLOSUM, 28, 1, prior)
        assert s >= 0
        assert gapblock_score(q, q, BLOSUM, 28, 1, prior)[0] == 0
