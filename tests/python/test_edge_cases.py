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


def test_invalid_engine():
    with pytest.raises(ValueError):
        seqtree.SearchParams(engine="bogus")


def test_invalid_mode_setter():
    p = seqtree.SearchParams()
    with pytest.raises(ValueError):
        p.mode = "sideways"


def test_invalid_matrix():
    idx = seqtree.Index.build(["CAT"], alphabet="aa")
    p = seqtree.SearchParams(matrix="PAM250")
    with pytest.raises(ValueError):
        idx.search("CAT", p)


def test_blosum_requires_aa():
    idx = seqtree.Index.build(["ACGT"], alphabet="nt")
    p = seqtree.SearchParams(matrix="BLOSUM62", engine="seqtrie", max_penalty=10)
    with pytest.raises(ValueError):
        idx.search("ACGT", p)


def test_invalid_symbol_in_query():
    idx = seqtree.Index.build(["CAT"], alphabet="aa")
    with pytest.raises(Exception):  # invalid symbol -> error surfaces
        idx.search("CA1", seqtree.SearchParams())
