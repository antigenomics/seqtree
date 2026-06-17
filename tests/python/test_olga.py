"""Sanity tests on OLGA-generated TCR CDR3 sequences (realistic data)."""
import random

import seqtree
from _olga import generate


def test_aa_self_match_and_recall():
    aa = sorted({aa for _, aa in generate(200)})
    idx = seqtree.Index.build(aa, alphabet="aa")
    assert len(idx) == len(aa)

    # every reference matches itself exactly
    p0 = seqtree.SearchParams()
    for i, s in enumerate(aa):
        assert any(h.ref_id == i and h.score == 0 for h in idx.search(s, p0))

    # a 1-substitution mutant still finds its parent
    rng = random.Random(0)
    p = seqtree.SearchParams(max_subs=1, engine="seqtm")
    found = 0
    for i, s in enumerate(aa[:100]):
        j = rng.randrange(len(s))
        mut = s[:j] + rng.choice("ACDEFGHIKLMNPQRSTVWY") + s[j + 1:]
        if i in {h.ref_id for h in idx.search(mut, p)}:
            found += 1
    assert found == 100  # parent always within scope


def test_nt_error_correction_neighbors():
    # nucleotide CDR3s: a single-base error should be found within 1 edit.
    nt = sorted({nt for nt, _ in generate(200)})
    idx = seqtree.Index.build(nt, alphabet="nt")
    rng = random.Random(1)
    p = seqtree.SearchParams(max_subs=1, max_ins=1, max_dels=1, engine="seqtm")
    for i, s in enumerate(nt[:50]):
        j = rng.randrange(len(s))
        mut = s[:j] + rng.choice("ACGT") + s[j + 1:]
        ids = {h.ref_id for h in idx.search(mut, p)}
        assert i in ids
