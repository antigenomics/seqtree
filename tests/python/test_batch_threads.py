import random
import os
import subprocess
import sys

import pytest

import seqtree


def _key(hits):
    return sorted((h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels) for h in hits)


def _random_db(n, length=15, seed=0):
    rng = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    return ["".join(rng.choice(aa) for _ in range(length)) for _ in range(n)]


def test_batch_equals_serial():
    db = _random_db(500, seed=1)
    queries = _random_db(200, seed=2)
    idx = seqtree.Index.build(db, alphabet="aa")
    p = seqtree.SearchParams(max_subs=2, engine="seqtm")

    serial = [idx.search(q, p) for q in queries]
    batch = idx.search_batch(queries, p, threads=2)
    assert len(serial) == len(batch)
    for s, b in zip(serial, batch):
        assert _key(s) == _key(b)


def test_threads_determinism():
    db = _random_db(800, seed=3)
    queries = _random_db(300, seed=4)
    idx = seqtree.Index.build(db, alphabet="aa")
    p = seqtree.SearchParams(max_subs=1, max_ins=1, max_dels=1, engine="seqtm")

    one = idx.search_batch(queries, p, threads=1)
    many = idx.search_batch(queries, p, threads=2)
    assert [_key(r) for r in one] == [_key(r) for r in many]


def test_pairwise_symmetry_of_membership():
    a = _random_db(50, seed=5)
    b = _random_db(120, seed=6)
    p = seqtree.SearchParams(max_subs=2, engine="seqtm")

    # pairwise must be a-major regardless of which side is indexed internally.
    res = seqtree.pairwise_batch(a, b, p, alphabet="aa")
    assert len(res) == len(a)

    # cross-check against explicit index-on-b
    idx_b = seqtree.Index.build(b, alphabet="aa")
    expected = idx_b.search_batch(a, p, threads=2)
    assert [_key(r) for r in res] == [_key(r) for r in expected]


@pytest.mark.skipif(not hasattr(os, "sched_setaffinity"), reason="Linux affinity contract")
def test_native_batch_respects_restricted_allocation():
    subprocess.run([sys.executable, "-c", '''
import os
os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
import seqtree
index = seqtree.Index.build(["AAAA", "AAAC"])
params = seqtree.SearchParams(max_subs=1, max_total_edits=1)
queries = ["AAAA"] * 8
assert index.edit_histogram_batch(queries, params, threads=0) == [[1, 1]] * 8
one = index.search_batch(queries, params, threads=1)
all_available = index.search_batch(queries, params, threads=0)
assert [[(h.ref_id, h.score) for h in row] for row in one] == [[(h.ref_id, h.score) for h in row] for row in all_available]
try:
    index.search_batch(queries, params, threads=2)
except ValueError:
    pass
else:
    raise AssertionError("explicit oversubscription was accepted")
'''], check=True, timeout=30)
