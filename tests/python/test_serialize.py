import pytest

import seqtree


def test_save_load_roundtrip(tmp_path):
    refs = ["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF", "CASSPQGATNEKLFF", "CAT", "CAT"]
    idx = seqtree.Index.build(refs, alphabet="aa")
    path = str(tmp_path / "idx.sqtree")
    idx.save(path)

    idx2 = seqtree.Index.load(path)
    assert len(idx2) == len(idx)
    assert idx2.ref_seq(0) == idx.ref_seq(0)

    p = seqtree.SearchParams(max_subs=2, max_ins=1, max_dels=1, engine="seqtm")
    for q in refs:
        a = sorted((h.ref_id, h.score) for h in idx.search(q, p))
        b = sorted((h.ref_id, h.score) for h in idx2.search(q, p))
        assert a == b


def test_load_bad_file_raises(tmp_path):
    path = str(tmp_path / "bad.bin")
    with open(path, "wb") as f:
        f.write(b"not a seqtree index at all")
    with pytest.raises(Exception):
        seqtree.Index.load(path)


def test_load_missing_file_raises():
    with pytest.raises(Exception):
        seqtree.Index.load("/tmp/seqtree_definitely_missing_98765.sqtree")
