"""Hamming and Levenshtein distances, checked against pure-Python oracles.

Unit costs, no matrix, no alphabet -- these are plain string distances. The oracles here are
trivial reference implementations, so the tests are differential over random data and need no
external dependency (unlike test_pairwise.py, which oracles against BioPython).
"""
import random

import numpy as np
import pytest

import seqtree
from seqtree.distance import (
    hamming,
    hamming_matrix,
    levenshtein,
    levenshtein_matrix,
    neighbourhood,
    neighbourhood_union,
    union_size,
)

AA = "ACDEFGHIKLMNPQRSTVWY"


def ref_hamming(a, b):
    assert len(a) == len(b)
    return sum(x != y for x, y in zip(a, b))


def ref_levenshtein(a, b):
    m, n = len(a), len(b)
    prev = list(range(n + 1))
    for i in range(1, m + 1):
        cur = [i] + [0] * n
        for j in range(1, n + 1):
            cur[j] = min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] != b[j - 1]))
        prev = cur
    return prev[n]


def rand_seq(rng, alphabet, lo, hi):
    return "".join(rng.choice(alphabet) for _ in range(rng.randint(lo, hi)))


# --- scalar, known values ---------------------------------------------------------------

@pytest.mark.parametrize("a,b,d", [
    ("AAAA", "AAAA", 0),
    ("AAAA", "ATAT", 2),
    ("CASSLGQYF", "CASSPGQYF", 1),
    ("", "", 0),
])
def test_hamming_known(a, b, d):
    assert hamming(a, b) == d


@pytest.mark.parametrize("a,b,d", [
    ("kitten", "sitting", 3),
    ("flaw", "lawn", 2),
    ("CASSLGQAYEQYF", "CASSPGQAYEQF", 2),
    ("", "", 0),
    ("abc", "", 3),
    ("", "abc", 3),
    ("abc", "abc", 0),
])
def test_levenshtein_known(a, b, d):
    assert levenshtein(a, b) == d


def test_hamming_requires_equal_length():
    with pytest.raises(ValueError, match="equal-length"):
        hamming("AAAA", "AAA")


def test_case_sensitive():
    # Unlike the search Codec, these do NOT fold case.
    assert hamming("ABC", "abc") == 3
    assert levenshtein("ABC", "abc") == 3


# --- scalar, differential vs the oracle -------------------------------------------------

def test_hamming_matches_oracle_random():
    rng = random.Random(1)
    for _ in range(2000):
        L = rng.randint(0, 20)
        a, b = rand_seq(rng, AA, L, L), rand_seq(rng, AA, L, L)
        assert hamming(a, b) == ref_hamming(a, b)


def test_levenshtein_matches_oracle_random():
    rng = random.Random(2)
    for _ in range(2000):
        a = rand_seq(rng, "ACGT", 0, 18)
        b = rand_seq(rng, "ACGT", 0, 18)
        assert levenshtein(a, b) == ref_levenshtein(a, b)


def test_levenshtein_is_symmetric():
    rng = random.Random(3)
    for _ in range(500):
        a, b = rand_seq(rng, AA, 0, 15), rand_seq(rng, AA, 0, 15)
        assert levenshtein(a, b) == levenshtein(b, a)


def test_levenshtein_bounds():
    # 0 <= lev(a,b) <= max(len); and |len(a)-len(b)| <= lev <= max(len).
    rng = random.Random(4)
    for _ in range(500):
        a, b = rand_seq(rng, AA, 0, 15), rand_seq(rng, AA, 0, 15)
        d = levenshtein(a, b)
        assert abs(len(a) - len(b)) <= d <= max(len(a), len(b))


# --- matrices ---------------------------------------------------------------------------

def test_hamming_matrix_matches_oracle():
    rng = random.Random(5)
    L = 12
    a = [rand_seq(rng, AA, L, L) for _ in range(40)]
    b = [rand_seq(rng, AA, L, L) for _ in range(30)]
    d = np.asarray(hamming_matrix(a, b, threads=0))
    assert d.shape == (40, 30) and d.dtype == np.int32
    for i in range(40):
        for k in range(30):
            assert d[i, k] == ref_hamming(a[i], b[k])


def test_levenshtein_matrix_matches_oracle():
    rng = random.Random(6)
    a = [rand_seq(rng, "ACGT", 0, 15) for _ in range(35)]
    b = [rand_seq(rng, "ACGT", 0, 15) for _ in range(25)]
    d = np.asarray(levenshtein_matrix(a, b, threads=2))
    assert d.shape == (35, 25)
    for i in range(35):
        for k in range(25):
            assert d[i, k] == ref_levenshtein(a[i], b[k])


def test_matrix_self_has_zero_diagonal_and_is_symmetric():
    rng = random.Random(7)
    seqs = [rand_seq(rng, AA, 10, 10) for _ in range(20)]
    h = np.asarray(hamming_matrix(seqs, seqs))
    assert (h.diagonal() == 0).all() and (h == h.T).all()
    lev = [rand_seq(rng, AA, 0, 14) for _ in range(20)]
    lv = np.asarray(levenshtein_matrix(lev, lev))
    assert (lv.diagonal() == 0).all() and (lv == lv.T).all()


def test_matrix_thread_counts_agree():
    rng = random.Random(8)
    a = [rand_seq(rng, "ACGT", 0, 15) for _ in range(30)]
    b = [rand_seq(rng, "ACGT", 0, 15) for _ in range(20)]
    one = np.asarray(levenshtein_matrix(a, b, threads=1))
    many = np.asarray(levenshtein_matrix(a, b, threads=8))
    assert (one == many).all()


def test_hamming_matrix_length_mismatch_raises_cleanly():
    # A throw inside a worker thread must surface as a Python exception, not a SIGABRT.
    with pytest.raises(ValueError, match="equal-length"):
        hamming_matrix(["AAAA", "AAA"], ["AAAA"])


def test_matrix_empty_inputs():
    assert np.asarray(hamming_matrix([], ["A"])).shape == (0, 1)
    assert np.asarray(hamming_matrix(["A"], [])).shape == (1, 0)
    assert np.asarray(levenshtein_matrix([], [])).shape == (0, 0)


def test_module_is_exported():
    assert seqtree.distance is not None
    assert set(seqtree.distance.__all__) == {
        "hamming", "levenshtein", "hamming_matrix", "levenshtein_matrix",
        "neighbourhood", "neighbourhood_union", "union_size"}


# --- Hamming neighbourhoods -------------------------------------------------------------
#
# Sizes are exact, so these are equalities, not bounds. Over a k-letter alphabet a length-L
# sequence has (k-1)*L neighbours at distance exactly 1 -- one substitution per position,
# k-1 choices each -- hence 19*L over the 20 standard residues.

S = "CASSLGQYF"        # L = 9


def test_closed_ball_size_is_19L_plus_one():
    for s in [S, "CASSLAPGATNEKLFF", "A"]:
        assert len(neighbourhood(s, 1)) == 19 * len(s) + 1


def test_ball_without_the_centre_is_19L():
    for s in [S, "CASSLAPGATNEKLFF", "A"]:
        assert len(neighbourhood(s, 1, include_self=False)) == 19 * len(s)


def test_every_member_is_within_r_of_the_centre():
    # Cross-checked with the C++ hamming, not with the generator's own bookkeeping.
    for r in (1, 2):
        ball = neighbourhood(S, r)
        assert all(hamming(S, m) <= r for m in ball)
        assert max(hamming(S, m) for m in ball) == r


def test_ball_members_are_distinct_and_the_centre_comes_first():
    ball = neighbourhood(S, 2)
    assert len(set(ball)) == len(ball)
    assert ball[0] == S


def test_r2_ball_size_matches_the_closed_form():
    # |B_2| = 1 + 19*L + 19^2 * C(L, 2): choose the 2 substituted positions, 19 residues each.
    L = len(S)
    assert len(neighbourhood(S, 2)) == 1 + 19 * L + 19 ** 2 * (L * (L - 1) // 2)


def test_r0_is_just_the_centre():
    assert neighbourhood(S, 0) == [S]
    assert neighbourhood(S, 0, include_self=False) == []


def test_shell_is_the_sphere_at_exactly_r():
    shell = neighbourhood(S, 2, shell=True)
    assert all(hamming(S, m) == 2 for m in shell)
    assert len(shell) == 19 ** 2 * (len(S) * (len(S) - 1) // 2)
    # The ball is the disjoint union of its shells.
    assert (set(neighbourhood(S, 2))
            == set(neighbourhood(S, 0)) | set(neighbourhood(S, 1, shell=True)) | set(shell))


def test_custom_alphabet():
    assert len(neighbourhood("ACGT", 1, alphabet="ACGT")) == 3 * 4 + 1
    assert set(neighbourhood("AA", 1, alphabet="AB")) == {"AA", "BA", "AB"}


def test_negative_radius_raises():
    with pytest.raises(ValueError, match="non-negative"):
        neighbourhood(S, -1)


# --- the union --------------------------------------------------------------------------

def test_union_of_one_equals_the_single_neighbourhood():
    assert set(neighbourhood_union([S])) == set(neighbourhood(S))
    assert set(neighbourhood_union([S], 2)) == set(neighbourhood(S, 2))


def test_union_is_idempotent_over_a_repeated_sequence():
    assert set(neighbourhood_union([S, S])) == set(neighbourhood(S))
    assert len(neighbourhood_union([S, S, S])) == 19 * len(S) + 1


def test_union_of_a_hamming_1_pair_overlaps_by_exactly_20():
    # s and t differ at exactly one position p. A sequence u lies in B_1(s) & B_1(t) iff it
    # agrees with both outside p: if u differed from s at some q != p it would be 1 away from
    # s but 2 away from t (q and p), so only position p may vary. Varying p over the whole
    # alphabet gives 20 such u -- s and t themselves plus 18 others: the memo's "2 + 18".
    s, t = "CASSLGQYF", "CASSPGQYF"
    assert hamming(s, t) == 1
    a, b = neighbourhood(s), neighbourhood(t)
    u = neighbourhood_union([s, t])
    assert len(u) < len(a) + len(b)
    assert len(set(a) & set(b)) == 20
    assert len(u) == len(a) + len(b) - 20 == 38 * len(s) - 18
    assert set(u) == set(a) | set(b)


def test_union_over_mixed_lengths_cannot_overlap():
    # Substitution only: a length-9 ball and a length-10 ball are disjoint by construction.
    a, b = "CASSLGQYF", "CASSLGQYFF"
    assert len(neighbourhood_union([a, b])) == (19 * 9 + 1) + (19 * 10 + 1)


def test_union_without_the_centres_drops_every_centre():
    # include_self=False removes the r=0 shell, which for a union is exactly set(seqs) --
    # including t, even though t is a genuine distance-1 neighbour of s.
    s, t = "CASSLGQYF", "CASSPGQYF"
    u = neighbourhood_union([s, t], include_self=False)
    assert s not in u and t not in u
    assert len(u) == 38 * len(s) - 18 - 2


def test_union_shell_uses_the_distance_to_the_nearest_centre():
    # t is 1 from s, so t is a centre (shell 0) and never appears in the shell-1 set.
    s, t = "CASSLGQYF", "CASSPGQYF"
    shell = neighbourhood_union([s, t], 1, shell=True)
    assert t not in shell and s not in shell
    assert all(min(hamming(s, m), hamming(t, m)) == 1 for m in shell)
    assert len(shell) == 38 * len(s) - 18 - 2


def test_shells_partition_the_union_exactly():
    """The property downstream shell-profiling relies on: shells 0..r tile the closed ball.

    Disjoint, exhaustive, and each member sits at distance exactly d from its *nearest* centre
    -- checked against the C++ ``hamming``, not against the generator's own bookkeeping. Were a
    shell to silently drop or double-count members, a per-shell quantity summed back up would
    stop reproducing the ball, with nothing raised anywhere.
    """
    centres = ["CASSLGQ", "CATSLGQ", "CASSPGQ", "CASSLGY"]      # near-duplicates: balls overlap
    R = 2
    ball = neighbourhood_union(centres, R)
    shells = [set(neighbourhood_union(centres, d, shell=True)) for d in range(R + 1)]

    for i in range(R + 1):
        for j in range(i + 1, R + 1):
            assert not (shells[i] & shells[j]), f"shells {i} and {j} overlap"
    assert set().union(*shells) == set(ball)
    assert sum(len(s) for s in shells) == len(ball)
    assert union_size(centres, R) == sum(union_size(centres, d, shell=True) for d in range(R + 1))

    for d, members in enumerate(shells):
        for m in members:
            assert min(hamming(c, m) for c in centres) == d


def test_a_bare_string_is_rejected_rather_than_split_into_centres():
    """``set("CASSLGQYF")`` is 8 one-character centres, so the union would be 20, not 172.

    A plausible integer and no error -- exactly the shape of a silent wrong answer. Callers hold
    a single junction often enough that this has to raise rather than answer.
    """
    for fn in (neighbourhood_union, union_size):
        with pytest.raises(TypeError, match="not one string"):
            fn("CASSLGQYF")
    assert union_size(["CASSLGQYF"]) == 172      # the call that was meant


def test_union_size_agrees_with_the_materialised_union():
    rng = random.Random(11)
    seqs = [rand_seq(rng, AA, 8, 8) for _ in range(5)] + [S, S]
    for r in (0, 1, 2):
        for include_self in (True, False):
            for shell in (False, True):
                got = union_size(seqs, r, include_self=include_self, shell=shell)
                assert got == len(neighbourhood_union(seqs, r, include_self=include_self,
                                                      shell=shell))


def test_union_members_are_distinct():
    seqs = ["CASSLGQYF", "CASSPGQYF", "CASSLGQYY"]
    u = neighbourhood_union(seqs, 1)
    assert len(set(u)) == len(u)


def test_union_equals_the_python_set_oracle():
    rng = random.Random(12)
    seqs = [rand_seq(rng, AA, 6, 6) for _ in range(4)]
    oracle = set()
    for s in seqs:
        for i in range(len(s)):
            for a in AA:
                oracle.add(s[:i] + a + s[i + 1:])
    assert set(neighbourhood_union(seqs, 1)) == oracle
    assert union_size(seqs, 1) == len(oracle)


# --- degenerate inputs ------------------------------------------------------------------

def test_empty_input():
    assert neighbourhood_union([]) == []
    assert union_size([]) == 0
    assert neighbourhood_union([], 2, shell=True) == []


def test_empty_sequence():
    assert neighbourhood("", 1) == [""]
    assert neighbourhood("", 1, include_self=False) == []
    assert neighbourhood("", 1, shell=True) == []      # nowhere to substitute


def test_single_character_sequence():
    assert len(neighbourhood("A", 1)) == 20
    assert set(neighbourhood("A", 1)) == set(AA)
    assert len(neighbourhood("A", 2)) == 20            # radius past the diameter adds nothing
