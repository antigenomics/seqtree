"""Regression tests for the engine-routing and scoring bugs fixed on feature/engine-bugfix.

Each test pins one previously-silent wrong answer. See ROADMAP.md and the commit message.
"""
import pytest

import seqtree as st

REFS = ["CASSLGQAYEQYF", "CASSPGTGVYGYTF", "CASRQGAWDTQYF", "CASSFRSSYNEQFF",
        "CAWSVSGGGTDTQYF", "CASSLAPGATNEKLFF", "CSARDRTGNGYTF", "CASSIRSSYEQYF"]
QUERY = "CASSLGQAYEQYF"


@pytest.fixture
def idx():
    return st.Index.build(list(REFS), "aa")


def test_matrix_with_indels_does_not_full_scan(idx):
    """B1: engine='auto' routed matrix+indel into seqtrie, whose budget defaulted to +inf,
    so every reference was returned regardless of the edit caps."""
    p = st.SearchParams(max_subs=1, max_ins=1, max_dels=0, matrix="BLOSUM62")
    hits = idx.search(QUERY, p)
    assert len(hits) < len(REFS)
    assert [idx.ref_seq(h.ref_id) for h in hits] == [QUERY]


def test_seqtrie_with_matrix_requires_explicit_budget(idx):
    """B1: seqtrie cuts only on accumulated penalty; without one it walks the whole trie."""
    with pytest.raises(ValueError, match="max_penalty"):
        idx.search(QUERY, st.SearchParams(max_ins=1, matrix="BLOSUM62", engine="seqtrie"))
    # ...and with a budget it works.
    assert idx.search(QUERY, st.SearchParams(matrix="BLOSUM62", max_penalty=1, engine="seqtrie"))


def test_auto_reports_exact_edit_counts_with_a_matrix():
    """B2: seqtrie ignores per-type caps and reports n_subs=n_ins=n_dels=0. Auto must not
    pick it, so a matrix + indel search still yields exact, per-type edit counts."""
    idx = st.Index.build(["CASSLGQAYEQYF", "CASSLGQAYEQYFF", "CASSLGQAYEQY"], "aa")
    p = st.SearchParams(max_subs=0, max_ins=1, max_dels=1, max_total_edits=1, matrix="BLOSUM62")
    got = {idx.ref_seq(h.ref_id): (h.n_subs, h.n_ins, h.n_dels) for h in idx.search(QUERY, p)}
    assert got == {"CASSLGQAYEQYF": (0, 0, 0),    # exact
                   "CASSLGQAYEQYFF": (0, 1, 0),   # ref longer -> insertion
                   "CASSLGQAYEQY": (0, 0, 1)}     # ref shorter -> deletion


def test_local_mode_is_gone():
    """B3: mode='local' was a silent no-op identical to 'all', with zero call sites."""
    with pytest.raises(ValueError, match="unknown mode"):
        st.SearchParams(mode="local")


def test_gap_extend_affects_the_score(idx):
    """B4: gap_extend was plumbed to C++ and read by no scoring code."""
    two_gaps = "CASSLGYEQYF"  # 11 aa vs a 13 aa ref -> one gap run of length 2
    cheap = idx.align(0, two_gaps, st.SearchParams(matrix="BLOSUM62", gap_open=11, gap_extend=1))
    dear = idx.align(0, two_gaps, st.SearchParams(matrix="BLOSUM62", gap_open=11, gap_extend=9))
    assert cheap.score < dear.score
    assert dear.score - cheap.score == 8  # one extension column, 9 - 1


def test_affine_reduces_to_linear_when_open_equals_extend(idx):
    """B4: gap_open == gap_extend must reproduce the linear-gap NW this replaced."""
    a = idx.align(0, "CASSLGYEQYF", st.SearchParams(gap_open=1, gap_extend=1))
    assert a.score == 2 and a.ops.count("I") + a.ops.count("D") == 2


def test_scale_is_the_matrix_unit(idx):
    """B6: the Gram transform puts a typical BLOSUM62 mismatch at ~14, so the default
    gap_open=1 made the aligner gap rather than substitute on an equal-length pair."""
    m = st.SubstitutionMatrix.blosum62()
    assert m.scale() == 14
    assert st.SubstitutionMatrix.unit(20).scale() == 1

    equal_len = "CASSIRSSYEQYF"  # same length as REFS[0]; needs no indel at all
    degenerate = idx.align(0, equal_len, st.SearchParams(matrix="BLOSUM62", gap_open=1))
    sane = idx.align(0, equal_len, st.SearchParams(matrix="BLOSUM62", gap_open=m.scale()))
    assert degenerate.ops.count("I") > 0          # the old default: gaps, zero substitutions
    assert sane.ops == "M" * 4 + "S" * 4 + "M" * 5


def test_pairwise_batch_transpose_swaps_ins_and_dels():
    """B8: the a-major transpose re-emitted n_ins/n_dels unswapped, so the reported edit
    direction was inverted -- but only when len(a) >= len(b) (the transpose path)."""
    a, b = ["CASSLGQAYEQYF"], ["CASSLGQAYEQYFF"]  # 13 vs 14
    p = st.SearchParams(max_subs=0, max_ins=1, max_dels=1, max_total_edits=1)

    truth = st.Index.build(list(b), "aa").search(a[0], p)[0]  # query a -> ref b
    assert (truth.n_ins, truth.n_dels) == (1, 0)

    ab = st.pairwise_batch(a, b, p, "aa")[0][0]  # hit describes a[0] -> b[0]
    assert (ab.n_ins, ab.n_dels) == (truth.n_ins, truth.n_dels)

    ba = st.pairwise_batch(b, a, p, "aa")[0][0]  # the mirror
    assert (ba.n_ins, ba.n_dels) == (ab.n_dels, ab.n_ins)


def test_align_validates_the_query_alphabet(idx):
    """B9: align() indexed the matrix with Codec::kInvalid for out-of-alphabet symbols."""
    with pytest.raises(ValueError, match="invalid symbol"):
        idx.align(0, "CASS!GQAYEQYF", st.SearchParams())


def test_negative_gap_costs_rejected(idx):
    """Scores must stay non-negative: the ball, scope monotonicity and trie pruning depend on it."""
    with pytest.raises(ValueError, match=">= 0"):
        idx.align(0, "CASSLGYEQYF", st.SearchParams(gap_open=-1))


def test_positional_matrix_from_tables_is_bound():
    """from_tables existed in C++ but was never exposed; the two-regime (flank | core)
    frame needs it."""
    size, width = 4, 3
    pm = st.PositionalMatrix.from_tables(size, width, [0] * (width * size * size), [1, 0, 0])
    assert (pm.size(), pm.width()) == (size, width)
    assert pm.masked(0) and not pm.masked(1)
