"""Hamming and Levenshtein distances, checked against pure-Python oracles.

Unit costs, no matrix, no alphabet -- these are plain string distances. The oracles here are
trivial reference implementations, so the tests are differential over random data and need no
external dependency (unlike test_pairwise.py, which oracles against BioPython).
"""
import random

import numpy as np
import pytest

import seqtree
from seqtree.distance import hamming, hamming_matrix, levenshtein, levenshtein_matrix

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
        "hamming", "levenshtein", "hamming_matrix", "levenshtein_matrix"}
