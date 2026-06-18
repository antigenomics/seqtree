import random

import seqtree
from seqtree.evalue import _poisson_sf, evalues

AA = "ACDEFGHIKLMNPQRSTVWY"


def _rand(n, rng, length=14):
    return ["".join(rng.choice(AA) for _ in range(length)) for _ in range(n)]


def test_poisson_sf_basics():
    assert _poisson_sf(0, 5.0) == 1.0
    assert _poisson_sf(1, 0.0) == 0.0
    # P(Poisson(1) >= 1) = 1 - e^-1
    assert abs(_poisson_sf(1, 1.0) - (1 - 2.718281828 ** -1)) < 1e-6


def test_planted_cluster_is_significant_background_is_not():
    rng = random.Random(0)
    control = seqtree.Index.build(_rand(20000, rng), alphabet="aa")

    q = "".join(rng.choice(AA) for _ in range(14))
    cluster = []
    for _ in range(50):
        s = list(q)
        s[rng.randrange(14)] = rng.choice(AA)
        cluster.append("".join(s))
    target = seqtree.Index.build(_rand(20000, rng) + cluster, alphabet="aa")

    p = seqtree.SearchParams(max_subs=1, engine="seqtm")
    r = evalues(target, control, [q], p)[0]
    assert r["n_target"] >= 40           # the planted neighbors are found
    assert r["p_enrichment"] < 1e-3      # far more neighbors than the background predicts

    qb = "".join(rng.choice(AA) for _ in range(14))  # a random background query
    rb = evalues(target, control, [qb], p)[0]
    assert rb["p_enrichment"] > 0.01     # not significant


def test_exclude_exact_drops_self_hit():
    rng = random.Random(7)
    q = "CASSPGTEAFF"
    mut = []
    for _ in range(2):
        s = list(q)
        s[rng.randrange(len(q))] = rng.choice(AA)
        mut.append("".join(s))
    target = seqtree.Index.build([q] + mut, alphabet="aa")  # q is a member of the target
    control = seqtree.Index.build(_rand(5000, rng), alphabet="aa")
    p = seqtree.SearchParams(max_subs=1, engine="seqtm")

    full = evalues(target, control, [q], p, exclude_exact=False)[0]
    punc = evalues(target, control, [q], p, exclude_exact=True)[0]
    assert punc["n_target"] == full["n_target"] - 1  # the exact self-match is removed


def test_evalue_monotone_in_scope():
    rng = random.Random(1)
    db = _rand(10000, rng)
    target = seqtree.Index.build(db, alphabet="aa")
    control = seqtree.Index.build(_rand(10000, rng), alphabet="aa")
    q = [db[0]]
    e1 = evalues(target, control, q, seqtree.SearchParams(max_subs=1, engine="seqtm"))[0]["E"]
    e2 = evalues(target, control, q, seqtree.SearchParams(max_subs=3, engine="seqtm"))[0]["E"]
    assert e2 >= e1
