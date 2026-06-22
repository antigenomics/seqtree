"""Substitution-matrix support: the built-in matrix list, custom matrices, and
seqtm matrix scoring. Built-ins: identity, BLOSUM62, PAM250, PAM100, structural."""
import pytest

import seqtree as st

BUILTINS = ["identity", "BLOSUM62", "PAM250", "PAM100", "structural"]


def test_amino_acid_order():
    aa = st.amino_acids()
    assert aa.startswith("ARNDCQEGHILKMFPSTWYV")  # BLOSUM62 / PAM column order
    assert len(aa) == 24
    assert st.alphabet_symbols("nt") == "ACGT"


@pytest.mark.parametrize("name", BUILTINS)
def test_named_matrix_exact_match_is_zero(name):
    idx = st.Index.build(["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF"], alphabet="aa")
    p = st.SearchParams(matrix=name, max_total_edits=3, max_penalty=80, engine="seqtrie")
    hits = {h.ref_id: h.score for h in idx.search("CASSLAPGATNEKLFF", p)}
    assert hits[0] == 0  # exact match -> zero penalty


# seqtm Hamming-ball score for the two-substitution ref (A->E, P->L), per matrix.
# Each is sum of Gram penalties sim(a,a)+sim(b,b)-2*sim(a,b); identity is plain edit cost.
MATRIX_SCORES = {"identity": 2, "structural": 18, "PAM250": 24, "BLOSUM62": 28, "PAM100": 30}


@pytest.mark.parametrize("name,expected", MATRIX_SCORES.items())
def test_named_matrix_seqtm_score(name, expected):
    idx = st.Index.build(["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF"], alphabet="aa")
    q = "CASSLAPGATNEKLFF"
    p = st.SearchParams(matrix=name, max_subs=3, engine="seqtm")
    score = {h.ref_id: h.score for h in idx.search(q, p)}
    assert score[0] == 0
    assert score[1] == expected
    # High gap cost forbids gaps, so the C++ global alignment matches the Hamming score.
    pa = st.SearchParams(matrix=name, gap_open=100, engine="seqtm")
    assert idx.align(1, q, pa).score == expected


def test_builtin_matrices_are_distinct():
    assert len(set(MATRIX_SCORES.values())) == len(MATRIX_SCORES)


def test_identity_is_edit_cost():
    # identity works on any alphabet and is plain unit (Hamming) cost.
    idx = st.Index.build(["ACGTACGT"], alphabet="nt")
    p = st.SearchParams(matrix="identity", max_subs=2, engine="seqtm")
    assert {h.ref_id: h.score for h in idx.search("ACGAACGT", p)}[0] == 1  # one substitution


@pytest.mark.parametrize("name", ["BLOSUM62", "PAM250", "PAM100", "structural"])
def test_aa_matrix_rejects_nt_alphabet(name):
    idx = st.Index.build(["ACGTACGT"], alphabet="nt")
    with pytest.raises(ValueError):
        idx.search("ACGTACGT", st.SearchParams(matrix=name, max_total_edits=1, engine="seqtrie"))


def test_unknown_matrix_name_rejected():
    with pytest.raises(ValueError):
        st.SearchParams(matrix="PAM50")  # removed from the built-in list


def test_custom_matrix_matches_builtin():
    # A custom SubstitutionMatrix built from BLOSUM62's object must score identically
    # to the named builtin.
    idx = st.Index.build(["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF", "CASSPRDGATNEKLFF"], alphabet="aa")
    q = "CASSLAPGATNEKLFF"
    obj = st.SubstitutionMatrix.blosum62()
    assert obj.size() == 24
    a = [h.score for h in idx.search(q, st.SearchParams(matrix=obj, max_total_edits=3, engine="seqtrie"))]
    b = [h.score for h in idx.search(q, st.SearchParams(matrix="BLOSUM62", max_total_edits=3, engine="seqtrie"))]
    assert a == b


def test_matrix_size_mismatch_raises():
    idx = st.Index.build(["ACGTACGT"], alphabet="nt")
    bad = st.SubstitutionMatrix.blosum62()  # 24 symbols, but nt alphabet has 4
    with pytest.raises(ValueError):
        idx.search("ACGTACGT", st.SearchParams(matrix=bad, max_total_edits=1, engine="seqtrie"))


def test_substitution_matrix_penalty_scalar():
    # Scalar amino-acid penalty lookup (Gram distance): 0 on identity, larger for dissimilar.
    bl = st.SubstitutionMatrix.blosum62()
    assert bl.penalty("L", "L") == 0
    assert bl.penalty("W", "W") == 0
    assert bl.penalty("L", "I") < bl.penalty("L", "P")  # conservative < dissimilar
    # The scalar penalties are exactly the engine's per-position Gram costs: the A->E, P->L
    # double substitution scores 28 under seqtm (see MATRIX_SCORES["BLOSUM62"]).
    assert bl.penalty("A", "E") + bl.penalty("P", "L") == 28
    with pytest.raises(ValueError):
        bl.penalty("L", "AB")   # not a single residue
    with pytest.raises(ValueError):
        bl.penalty("L", "?")    # unknown amino acid


def test_collisions_only_with_indels():
    import random

    rng = random.Random(0)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    refs = ["".join(rng.choice(aa) for _ in range(14)) for _ in range(5000)]
    idx = st.Index.build(refs, alphabet="aa")
    q = [rng.choice(refs) for _ in range(200)]
    # Substitution-only (Hamming) reaches each reference by exactly one alignment.
    assert sum(idx.collisions_batch(q, st.SearchParams(max_subs=3, engine="seqtm"))) == 0
    # seqtrie reaches each reference once (one leaf), so never collides.
    assert sum(idx.collisions_batch(q, st.SearchParams(max_total_edits=3, engine="seqtrie"))) == 0
    # Indels create multiple edit paths to the same reference -> collisions.
    c = idx.collisions_batch(
        q, st.SearchParams(max_subs=3, max_ins=2, max_dels=2, max_total_edits=3, engine="seqtm"))
    assert sum(c) > 0
