import pytest

import seqtree


def test_empty_index():
    idx = seqtree.Index.build([], alphabet="aa")
    assert len(idx) == 0
    assert idx.search("CAT", seqtree.SearchParams()) == []


def test_empty_query():
    idx = seqtree.Index.build(["A", "AA"], alphabet="aa")
    p = seqtree.SearchParams(max_ins=2, engine="seqtm")
    ids = {h.ref_id for h in idx.search("", p)}
    assert ids == {0, 1}  # both reachable by insertions from empty query


def test_query_longer_than_refs():
    idx = seqtree.Index.build(["CAT"], alphabet="aa")
    p = seqtree.SearchParams(max_dels=3, engine="seqtm")
    hits = idx.search("CATWWW", p)
    assert any(h.ref_id == 0 and h.n_dels == 3 for h in hits)


def test_duplicate_refs_distinct_ids():
    idx = seqtree.Index.build(["CAT", "CAT"], alphabet="aa")
    hits = idx.search("CAT", seqtree.SearchParams())
    assert sorted(h.ref_id for h in hits) == [0, 1]


def test_empty_reference_in_db():
    idx = seqtree.Index.build(["", "A", "AA"], alphabet="aa")
    assert len(idx) == 3
    assert idx.ref_seq(0) == ""
    assert {h.ref_id for h in idx.search("", seqtree.SearchParams())} == {0}


def test_long_sequence_one_substitution():
    aa = "ACDEFGHIKLMNPQRSTVWY"
    ref = "".join(aa[i % 20] for i in range(2000))
    idx = seqtree.Index.build([ref], alphabet="aa")
    mut = ref[:1000] + ("C" if ref[1000] != "C" else "A") + ref[1001:]
    for eng, p in (
        ("seqtm", seqtree.SearchParams(max_subs=1, engine="seqtm")),
        ("seqtrie", seqtree.SearchParams(max_total_edits=1, engine="seqtrie")),
    ):
        hits = idx.search(mut, p)
        assert any(h.ref_id == 0 for h in hits), eng


def test_homopolymer_indels():
    idx = seqtree.Index.build(["AAAAAAAA", "AAAAAAA", "AAAAAAAAA"], alphabet="aa")
    p = seqtree.SearchParams(max_ins=1, max_dels=1, engine="seqtm")
    assert {h.ref_id for h in idx.search("AAAAAAAA", p)} == {0, 1, 2}


def test_blosum_penalty_scale():
    # squared-distance scale: CASSLELGATNEKLFF is 2 subs (A->E, P->L) from ref 0.
    idx = seqtree.Index.build(["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF"], alphabet="aa")
    p = seqtree.SearchParams(matrix="BLOSUM62", max_subs=2, engine="seqtm")
    score = {h.ref_id: h.score for h in idx.search("CASSLAPGATNEKLFF", p)}
    assert score[0] == 0
    # A->E: 4+5-2*(-1)=11 ; P->L: 7+4-2*(-3)=17 ; total 28
    assert score[1] == 28


def test_invalid_engine():
    with pytest.raises(ValueError):
        seqtree.SearchParams(engine="bogus")


def test_invalid_mode_setter():
    p = seqtree.SearchParams()
    with pytest.raises(ValueError):
        p.mode = "sideways"


def test_invalid_matrix():
    # An unknown matrix name is rejected eagerly at construction (like engine/mode).
    with pytest.raises(ValueError):
        seqtree.SearchParams(matrix="PAM30")  # not a built-in


def test_blosum_requires_aa():
    idx = seqtree.Index.build(["ACGT"], alphabet="nt")
    p = seqtree.SearchParams(matrix="BLOSUM62", engine="seqtrie", max_penalty=10)
    with pytest.raises(ValueError):
        idx.search("ACGT", p)


def test_invalid_symbol_in_query():
    idx = seqtree.Index.build(["CAT"], alphabet="aa")
    with pytest.raises(Exception):  # invalid symbol -> error surfaces
        idx.search("CA1", seqtree.SearchParams())
