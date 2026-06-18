import seqtree


def test_bundled_control_loads_and_caches(tmp_path):
    idx = seqtree.load_control(size=5000, cache_dir=str(tmp_path))
    assert len(idx) == 5000
    # cache file written; second call reloads it
    assert any(p.suffix == ".sqtree" for p in tmp_path.iterdir())
    idx2 = seqtree.load_control(size=5000, cache_dir=str(tmp_path))
    assert len(idx2) == 5000


def test_control_members_are_findable(tmp_path):
    idx = seqtree.load_control(size=2000, cache_dir=str(tmp_path))
    p = seqtree.SearchParams(max_subs=1, engine="seqtm")
    member = idx.ref_seq(0)
    assert any(h.ref_id == 0 and h.score == 0 for h in idx.search(member, p))


def test_unknown_control_raises():
    import pytest

    with pytest.raises(ValueError):
        seqtree.load_control(name="klingon_trb_aa", size=10)
