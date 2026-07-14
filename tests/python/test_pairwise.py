"""Needleman-Wunsch / Smith-Waterman, pinned against BioPython as an oracle.

The point of this module is to *be* BioPython for ordinary protein alignment, so the tests are
differential wherever BioPython is installed. It is an optional test dependency, never a runtime
one -- the invariant tests below run regardless, and the oracle tests skip.

Every convention here is somewhere a reimplementation goes quietly wrong:
  * a gap of length L costs open + (L-1)*extend, NOT open + L*extend
  * open == extend is linear gaps
  * global charges end gaps (true NW, not semi-global)
  * local floors at zero and takes the best cell anywhere (SW)
"""
import random

import pytest

import seqtree
from seqtree.pairwise import align, dist_matrix, score, score_matrix

AA = "ACDEFGHIKLMNPQRSTVWY"
MATRICES = {
    "BLOSUM62": seqtree.SubstitutionMatrix.blosum62(),
    "BLOSUM45": seqtree.SubstitutionMatrix.blosum45(),
    "BLOSUM80": seqtree.SubstitutionMatrix.blosum80(),
}
# mode, gap_open, gap_extend -- affine, linear (open == extend), free, and extreme
SETTINGS = [
    ("global", 11, 1), ("global", 12, 1), ("global", 5, 5), ("global", 1, 1),
    ("global", 0, 0), ("global", 100, 5),
    # gap_extend > gap_open is legal and makes closing-and-reopening cheaper than extending.
    # It is the regime that exposed a wrong re-scorer, and it must stay in the grid.
    ("global", 1, 3), ("global", 0, 5), ("local", 1, 10),
    # huge costs: these used to overflow int32 and come back POSITIVE.
    ("global", 10 ** 8, 1), ("global", 2 ** 30, 1),
    ("local", 11, 1), ("local", 5, 5), ("local", 1, 1), ("local", 12, 1),
]


def rand_seq(rng, lo, hi):
    return "".join(rng.choice(AA) for _ in range(rng.randint(lo, hi)))


def sample_pairs(rng, n=60):
    """A spread of shapes, because the bugs live at the edges."""
    pairs = []
    for _ in range(n // 2):
        pairs.append((rand_seq(rng, 5, 25), rand_seq(rng, 5, 25)))
    for _ in range(n // 6):
        s = rand_seq(rng, 10, 20)
        pairs.append((s, s))                                    # identical
    for _ in range(n // 6):
        s = rand_seq(rng, 12, 20)
        i, d = rng.randrange(1, len(s) - 1), rng.randint(1, 4)
        pairs.append((s, s[:i] + s[i + d:]))                    # one indel
    for _ in range(n // 6):
        pairs.append((rand_seq(rng, 1, 3), rand_seq(rng, 15, 25)))   # wildly unequal
    return pairs


# ---------------------------------------------------------------------------------------------
# the oracle: BioPython, if it is installed
# ---------------------------------------------------------------------------------------------

def _oracle(mode, matrix_name, gap_open, gap_extend):
    Align = pytest.importorskip("Bio.Align", reason="BioPython is an optional test-only oracle")
    from Bio.Align import substitution_matrices

    a = Align.PairwiseAligner()
    a.mode = mode
    a.substitution_matrix = substitution_matrices.load(matrix_name)
    a.open_gap_score = -gap_open
    a.extend_gap_score = -gap_extend
    return a


@pytest.mark.parametrize("matrix_name", sorted(MATRICES))
@pytest.mark.parametrize("mode,gap_open,gap_extend", SETTINGS)
def test_matches_biopython_exactly(matrix_name, mode, gap_open, gap_extend):
    """The whole point: drop-in for Bio.Align.PairwiseAligner on protein sequences."""
    oracle = _oracle(mode, matrix_name, gap_open, gap_extend)
    mat = MATRICES[matrix_name]
    rng = random.Random(hash((matrix_name, mode, gap_open, gap_extend)) & 0xFFFF)

    for q, r in sample_pairs(rng):
        expected = int(oracle.align(q, r).score)
        got = score(q, r, mat, mode=mode, gap_open=gap_open, gap_extend=gap_extend)
        assert got == expected, f"{matrix_name} {mode} open={gap_open} ext={gap_extend}: {q} vs {r}"


def test_matches_biopython_on_real_germline_distances():
    """The use case that motivated this: d = s(a,a) + s(b,b) - 2 s(a,b) over V genes."""
    Align = pytest.importorskip("Bio.Align", reason="BioPython is an optional test-only oracle")
    import numpy as np

    oracle = Align.PairwiseAligner("blastp")          # global, BLOSUM62, open 12, extend 1
    mat = MATRICES["BLOSUM62"]
    rng = random.Random(3)
    seqs = [rand_seq(rng, 60, 95) for _ in range(12)]

    D = np.asarray(dist_matrix(seqs, seqs, mat, gap_open=12, gap_extend=1))
    self_ = [oracle.align(s, s).score for s in seqs]
    for i, a in enumerate(seqs):
        for j, b in enumerate(seqs):
            expected = self_[i] + self_[j] - 2 * oracle.align(a, b).score
            assert D[i, j] == int(expected)


# ---------------------------------------------------------------------------------------------
# invariants -- these hold with or without the oracle
# ---------------------------------------------------------------------------------------------

def test_a_gap_of_length_L_costs_open_plus_L_minus_1_extends():
    """Not open + L*extend. Getting this backwards is the classic off-by-one."""
    m = MATRICES["BLOSUM62"]
    s_aa = m.similarity("A", "A")
    assert s_aa == 4
    for L in (1, 2, 3, 5):
        got = score("AAA", "AAA" + "A" * L, m, gap_open=11, gap_extend=1)
        assert got == 3 * s_aa - (11 + (L - 1) * 1)


def test_linear_gaps_are_just_open_equal_to_extend():
    m = MATRICES["BLOSUM62"]
    for L in (1, 2, 3):
        assert score("AAA", "AAA" + "A" * L, m, gap_open=5, gap_extend=5) == 12 - 5 * L


def test_global_charges_end_gaps():
    """True Needleman-Wunsch. Semi-global/overlap alignment would leave them free."""
    m = MATRICES["BLOSUM62"]
    got = score("AAA", "KKAAA", m, gap_open=11, gap_extend=1)
    assert got == 12 - (11 + 1)          # 3 A/A matches, one leading gap of length 2
    assert got != 12                     # a free end gap would give this


def test_local_never_goes_below_zero_and_finds_the_core():
    m = MATRICES["BLOSUM62"]
    assert score("WWW", "KKK", m, mode="local") == 0
    assert score("WWWAAAWWW", "KKKAAAKKK", m, mode="local") == 12   # the AAA core only
    assert score("WWWAAAWWW", "KKKAAAKKK", m, mode="global") < 0


def test_local_and_global_agree_when_the_whole_alignment_is_positive():
    m = MATRICES["BLOSUM62"]
    s = "CASSLGQAYEQYF"
    assert score(s, s, m, mode="local") == score(s, s, m, mode="global")


def test_identity_maximises_the_score():
    rng = random.Random(11)
    m = MATRICES["BLOSUM62"]
    for _ in range(50):
        a, b = rand_seq(rng, 8, 20), rand_seq(rng, 8, 20)
        assert score(a, a, m) >= score(a, b, m)


def test_score_is_symmetric():
    rng = random.Random(13)
    m = MATRICES["BLOSUM62"]
    for _ in range(50):
        a, b = rand_seq(rng, 5, 20), rand_seq(rng, 5, 20)
        for mode in ("global", "local"):
            assert score(a, b, m, mode=mode) == score(b, a, m, mode=mode)


def rescore(aln, matrix, gap_open, gap_extend):
    """Re-score an alignment from its columns, the honest way.

    A gap run only continues if it is the **same kind** of gap. An X-gap (a query residue against
    ``-``) immediately followed by a Y-gap (``-`` against a ref residue) is *two* opens, not one
    open and one extend -- they are different states in the recurrence.

    That distinction is not academic: with ``gap_extend > gap_open`` it is cheaper to close and
    reopen than to extend, so adjacent gaps of opposite type really do occur. An earlier version
    of this helper ignored it, which made this test pass vacuously at 11/1 and report 178 phantom
    failures at 1/3 -- against a traceback that was correct all along.
    """
    total, prev = 0, None            # prev in {None, "X", "Y"}
    for x, y in zip(aln.aligned_query, aln.aligned_ref):
        if y == "-":                 # query residue against a gap
            total -= gap_extend if prev == "X" else gap_open
            prev = "X"
        elif x == "-":               # gap against a ref residue
            total -= gap_extend if prev == "Y" else gap_open
            prev = "Y"
        else:
            total += matrix.similarity(x, y)
            prev = None
    return total


@pytest.mark.parametrize("mode", ["global", "local"])
@pytest.mark.parametrize("go,ge", [(11, 1), (5, 5), (1, 3), (0, 0), (1, 10), (2, 100)])
def test_align_reproduces_its_own_score(mode, go, ge):
    """The traceback must be an alignment that actually scores what was reported.

    ``gap_extend > gap_open`` is in the grid on purpose: it is legal, it makes closing and
    reopening cheaper than extending, and it is exactly the regime a naive re-scorer gets wrong.
    """
    rng = random.Random(hash((mode, go, ge)) & 0xFFFF)
    m = MATRICES["BLOSUM62"]
    for _ in range(40):
        a, b = rand_seq(rng, 4, 20), rand_seq(rng, 4, 20)
        aln = align(a, b, m, mode=mode, gap_open=go, gap_extend=ge)
        assert len(aln.aligned_query) == len(aln.aligned_ref) == len(aln.ops)
        assert rescore(aln, m, go, ge) == aln.score
        assert aln.score == score(a, b, m, mode=mode, gap_open=go, gap_extend=ge)
        if mode == "global":         # local returns the matched sub-sequences, not the whole input
            assert aln.aligned_query.replace("-", "") == a
            assert aln.aligned_ref.replace("-", "") == b


def test_huge_gap_costs_do_not_overflow_into_a_positive_score():
    """The DP used to run in int32 and wrap.

    At ``gap_open = 2**30`` the sentinel (INT32_MIN/4) minus the gap fell below INT32_MIN, so an
    alignment that should have scored about **-10^9** came back as a large **positive** number --
    a silently wrong answer, the worst kind. The DP is int64 now, so these are simply correct.
    """
    m = MATRICES["BLOSUM62"]
    q, r = "CASSLGQAYEQYF", "CASSPGQAYEQF"
    assert score(q, r, m, gap_open=10 ** 8) == -99_999_944
    assert score(q, r, m, gap_open=2 ** 30) == -1_073_741_768      # was +2_147_483_639
    assert score(q, r, m, gap_open=2 ** 31 - 1) == -2_147_483_591
    for go in (2 ** 30, 2 ** 31 - 1):
        assert score(q, r, m, gap_open=go) < 0, "a huge gap cost must not score positive"


def test_a_score_that_cannot_fit_int32_is_refused_not_truncated():
    """``Alignment.score`` and ``ScoreMatrix`` are int32. Where a result genuinely needs more,
    say so rather than hand back a wrapped number."""
    m = MATRICES["BLOSUM62"]
    with pytest.raises(OverflowError, match="32 bits"):
        score("A" * 400, "W" * 600, m, gap_open=2 ** 31 - 1, gap_extend=2 ** 31 - 1)


# ---------------------------------------------------------------------------------------------
# distances
# ---------------------------------------------------------------------------------------------

def test_dist_matrix_is_a_distance():
    rng = random.Random(19)
    m = MATRICES["BLOSUM62"]
    seqs = [rand_seq(rng, 10, 20) for _ in range(15)]
    D = dist_matrix(seqs, seqs, m)
    for i in range(len(seqs)):
        assert D[i, i] == 0                       # zero on the diagonal
        for j in range(len(seqs)):
            assert D[i, j] >= 0                   # non-negative
            assert D[i, j] == D[j, i]             # symmetric


def test_dist_matrix_equals_the_hand_rolled_gram_transform():
    """d = s(a,a) + s(b,b) - 2 s(a,b), the thing every BioPython user writes by hand."""
    rng = random.Random(23)
    m = MATRICES["BLOSUM62"]
    q = [rand_seq(rng, 8, 16) for _ in range(6)]
    r = [rand_seq(rng, 8, 16) for _ in range(5)]
    D = dist_matrix(q, r, m)
    for i, a in enumerate(q):
        for j, b in enumerate(r):
            expected = score(a, a, m) + score(b, b, m) - 2 * score(a, b, m)
            assert D[i, j] == expected


def test_score_matrix_matches_the_scalar_scorer():
    rng = random.Random(29)
    m = MATRICES["BLOSUM80"]
    q = [rand_seq(rng, 6, 18) for _ in range(8)]
    r = [rand_seq(rng, 6, 18) for _ in range(7)]
    for mode in ("global", "local"):
        sm = score_matrix(q, r, m, mode=mode, threads=2)
        assert sm.shape == (len(q), len(r))
        for i, a in enumerate(q):
            for j, b in enumerate(r):
                assert sm[i, j] == score(a, b, m, mode=mode)


def test_score_matrix_is_thread_invariant():
    rng = random.Random(31)
    m = MATRICES["BLOSUM62"]
    q = [rand_seq(rng, 8, 20) for _ in range(24)]
    r = [rand_seq(rng, 8, 20) for _ in range(12)]
    ref = [score_matrix(q, r, m, threads=1).row(i) for i in range(len(q))]
    for threads in (2, 4, 0):
        got = score_matrix(q, r, m, threads=threads)
        assert [got.row(i) for i in range(len(q))] == ref


# ---------------------------------------------------------------------------------------------
# matrices
# ---------------------------------------------------------------------------------------------

def test_blosum45_and_80_carry_the_documented_anchors():
    b45, b62, b80 = (MATRICES["BLOSUM45"], MATRICES["BLOSUM62"], MATRICES["BLOSUM80"])
    assert (b45.similarity("A", "A"), b45.similarity("W", "W")) == (5, 15)
    assert (b62.similarity("A", "A"), b62.similarity("W", "W")) == (4, 11)
    assert (b80.similarity("A", "A"), b80.similarity("W", "W")) == (7, 16)


def test_the_similarity_view_is_signed_and_the_penalty_view_is_not():
    """Both views of one matrix. The Gram transform is lossy, so both must be stored."""
    m = MATRICES["BLOSUM62"]
    assert m.similarity("A", "K") < 0            # signed log-odds
    assert m.penalty("A", "K") > 0               # non-negative cost
    assert m.penalty("W", "W") == m.penalty("C", "C") == 0     # diagonal destroyed by the Gram
    assert m.similarity("W", "W") != m.similarity("C", "C")    # but recoverable from sim


def test_named_matrices_reach_the_search_path_too():
    for name in ("BLOSUM45", "BLOSUM80"):
        p = seqtree.SearchParams(max_subs=2, max_penalty=40, matrix=name, engine="seqtm")
        idx = seqtree.Index.build(["CASSLGQAYEQYF"], "aa")
        assert any(h.ref_id == 0 for h in idx.search("CASSLGQAYEQYF", p))


# ---------------------------------------------------------------------------------------------
# input handling
# ---------------------------------------------------------------------------------------------

def test_empty_sequences():
    m = MATRICES["BLOSUM62"]
    assert score("", "", m) == 0
    assert score("", "AAA", m) == -(11 + 2 * 1)     # global charges the gap
    assert score("", "AAA", m, mode="local") == 0   # local just takes nothing
    assert score_matrix([], ["A"], m).shape == (0, 1)


def test_bad_arguments_raise():
    m = MATRICES["BLOSUM62"]
    with pytest.raises(ValueError, match="gap_open"):
        score("AA", "AA", m, gap_open=-1)
    with pytest.raises(ValueError, match="mode"):
        score("AA", "AA", m, mode="semiglobal")
    with pytest.raises(ValueError, match="alphabet|symbol"):
        score("AA#", "AA", m)
