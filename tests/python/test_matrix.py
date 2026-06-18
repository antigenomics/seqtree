"""Substitution-matrix support: PAM50, custom matrices, and seqtm matrix scoring."""
import pytest

import seqtree as st


def test_amino_acid_order():
    aa = st.amino_acids()
    assert aa.startswith("ARNDCQEGHILKMFPSTWYV")  # BLOSUM62 / PAM50 column order
    assert len(aa) == 24
    assert st.alphabet_symbols("nt") == "ACGT"


@pytest.mark.parametrize("name", ["BLOSUM62", "PAM50"])
def test_named_matrix_exact_match_is_zero(name):
    idx = st.Index.build(["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF"], alphabet="aa")
    p = st.SearchParams(matrix=name, max_total_edits=3, engine="seqtrie")
    hits = {h.ref_id: h.score for h in idx.search("CASSLAPGATNEKLFF", p)}
    assert hits[0] == 0  # exact match -> zero penalty


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


def test_seqtm_matrix_score_is_sum_of_substitution_penalties():
    # PAM50 penalty(a,b) = max(sim(a,a), sim(b,b)) - sim(a,b).
    # For A->E: 7-(-1)=8; for P->L: 8-(-6)=14; total 22 for the two-substitution ref.
    idx = st.Index.build(["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF"], alphabet="aa")
    p = st.SearchParams(matrix="PAM50", max_subs=3, engine="seqtm")
    score = {h.ref_id: h.score for h in idx.search("CASSLAPGATNEKLFF", p)}
    assert score[0] == 0
    assert score[1] == 22

    # With a gap cost too high to ever use, the C++ global alignment of an equal-length
    # pair must agree with the seqtm Hamming-ball score (best score across alignments).
    pa = st.SearchParams(matrix="PAM50", gap_open=100, engine="seqtm")
    assert idx.align(1, "CASSLAPGATNEKLFF", pa).score == 22


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
