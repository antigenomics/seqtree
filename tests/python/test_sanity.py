"""Property-style sanity checks that should hold for any search result."""
import random

import seqtree

AA = "ACDEFGHIKLMNPQRSTVWY"


def _db(n, length=14, seed=0):
    rng = random.Random(seed)
    return ["".join(rng.choice(AA) for _ in range(length)) for _ in range(n)]


def test_hits_respect_scope():
    idx = seqtree.Index.build(_db(2000, seed=1), alphabet="aa")
    p = seqtree.SearchParams(max_subs=2, max_ins=1, max_dels=1, engine="seqtm")
    for q in _db(50, seed=2):
        for h in idx.search(q, p):
            assert h.n_subs <= 2 and h.n_ins <= 1 and h.n_dels <= 1


def test_score_equals_edits_unit_cost():
    idx = seqtree.Index.build(_db(2000, seed=3), alphabet="aa")
    p = seqtree.SearchParams(max_subs=2, max_ins=2, max_dels=2, engine="seqtm")
    for q in _db(50, seed=4):
        for h in idx.search(q, p):
            assert h.score == h.n_subs + h.n_ins + h.n_dels


def test_larger_scope_is_superset():
    db = _db(3000, seed=5)
    idx = seqtree.Index.build(db, alphabet="aa")
    queries = _db(40, seed=6)
    p1 = seqtree.SearchParams(max_subs=1, engine="seqtm")
    p2 = seqtree.SearchParams(max_subs=2, engine="seqtm")
    for q in queries:
        s1 = {h.ref_id for h in idx.search(q, p1)}
        s2 = {h.ref_id for h in idx.search(q, p2)}
        assert s1 <= s2


def test_engines_agree_on_total_edit_budget():
    db = _db(3000, seed=7)
    idx = seqtree.Index.build(db, alphabet="aa")
    pm = seqtree.SearchParams(max_subs=2, max_ins=2, max_dels=2, max_total_edits=2, engine="seqtm")
    pt = seqtree.SearchParams(max_total_edits=2, engine="seqtrie")
    for q in _db(40, seed=8):
        m = {h.ref_id: h.score for h in idx.search(q, pm)}
        t = {h.ref_id: h.score for h in idx.search(q, pt)}
        assert m == t


def test_top_hit_is_global_minimum():
    db = _db(2000, seed=9)
    idx = seqtree.Index.build(db, alphabet="aa")
    p = seqtree.SearchParams(max_subs=3, engine="seqtm")
    for q in _db(40, seed=10):
        allhits = idx.search(q, p)
        top = idx.search_top(q, p, k=1)
        if allhits:
            assert top[0].score == min(h.score for h in allhits)


def test_blosum_score_is_monotone_in_budget():
    blosum_db = _db(2000, seed=11)
    idx = seqtree.Index.build(blosum_db, alphabet="aa")
    q = blosum_db[0]
    small = seqtree.SearchParams(matrix="BLOSUM62", max_penalty=8, engine="seqtrie", gap_open=4)
    big = seqtree.SearchParams(matrix="BLOSUM62", max_penalty=20, engine="seqtrie", gap_open=4)
    s = {h.ref_id for h in idx.search(q, small)}
    b = {h.ref_id for h in idx.search(q, big)}
    assert s <= b
