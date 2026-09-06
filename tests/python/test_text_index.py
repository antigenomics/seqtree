"""TextIndex: exact k-mismatch search over a concatenated text.

The C++ suite owns the completeness proof (a brute-force set-equality grid over L, max_subs
and k). These tests own the Python surface: the flat arrays and their zero-copy views, the
per-query object view, the named errors, and the guarantees a caller reads off the docs --
determinism across thread counts, ``truncated`` never being silent, and the whole nearest
shell coming back rather than one representative of it.
"""
import random

import pytest

from seqtree import TextIndex

REFS = ["CASSLGQYF", "CASSLGQYW", "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ", "CASSLGQWW"]


@pytest.fixture
def ix():
    return TextIndex.build(REFS, alphabet="aa", k=4)


def brute_force(refs, q, max_subs):
    """Every window of every record, scanned directly. The oracle, sharing no code."""
    out = set()
    for r, ref in enumerate(refs):
        for off in range(len(ref) - len(q) + 1):
            n = sum(a != b for a, b in zip(q, ref[off:off + len(q)]))
            if n <= max_subs:
                out.add((r, off, n))
    return out


def triples(res, i):
    return {(h.ref_id, h.offset, h.n_subs) for h in res[i]}


# --- the answer itself -------------------------------------------------------------------

def test_matches_a_brute_force_scan(ix):
    queries = ["CASSLGQYF", "CASSLGQYA", "KTAYIAKQR", "WWWWWWWWW"]
    for m in range(4):
        res = ix.search_batch(queries, max_subs=m)
        for i, q in enumerate(queries):
            assert triples(res, i) == brute_force(REFS, q, m), (q, m)


def test_matches_a_brute_force_scan_on_random_text():
    rng = random.Random(20260906)
    alpha = "ACDEFG"
    refs = ["".join(rng.choice(alpha) for _ in range(rng.randint(25, 60))) for _ in range(8)]
    ixr = TextIndex.build(refs, alphabet="aa", k=4)
    queries = [refs[rng.randrange(len(refs))][3:3 + L] for L in (6, 8, 11, 14, 20)]
    for m in range(4):
        res = ixr.search_batch(queries, max_subs=m)
        for i, q in enumerate(queries):
            assert triples(res, i) == brute_force(refs, q, m), (q, m)


def test_one_index_serves_every_query_length():
    """The whole point: Index needs one build per length, TextIndex needs one build."""
    ixr = TextIndex.build(REFS, alphabet="aa", k=4)
    lengths = [len(q) for q in ("KTAY", "KTAYI", "KTAYIAKQR", "KTAYIAKQRQISFVKSHFSRQ")]
    res = ixr.search_batch(["KTAY", "KTAYI", "KTAYIAKQR", "KTAYIAKQRQISFVKSHFSRQ"], max_subs=1)
    assert len(res) == 4
    assert all(res[i][0].ref_id == 2 for i in range(4)), lengths


# --- the guarantees ----------------------------------------------------------------------

def test_best_only_returns_the_whole_nearest_shell():
    # Two records sit at distance 2 from this query; best_only must return BOTH, not one.
    ixr = TextIndex.build(["CASSLGQYF", "CASSLGQAA", "CASSLGQWW"], alphabet="aa", k=4)
    res = ixr.search_batch(["CASSLGQYF"], max_subs=3, best_only=True, exclude_exact=True)
    assert [h.n_subs for h in res[0]] == [2, 2]
    assert {h.ref_id for h in res[0]} == {1, 2}


def test_a_cap_is_reported_and_keeps_the_best_hits(ix):
    res = ix.search_batch(["CASSLGQYF"], max_subs=2)
    assert res.num_hits == 3
    assert list(res.truncated) == [0]

    capped = ix.search_batch(["CASSLGQYF"], max_subs=2, max_hits=1)
    assert res.num_hits == 3 and capped.num_hits == 1
    assert list(capped.truncated) == [1]
    assert capped[0][0].n_subs == 0  # sorted before the cap, so the exact hit survives


def test_hit_order_is_deterministic_across_thread_counts():
    rng = random.Random(4)
    alpha = "ACDEFG"
    refs = ["".join(rng.choice(alpha) for _ in range(50)) for _ in range(60)]
    ixr = TextIndex.build(refs, alphabet="aa", k=4)
    queries = [r[:12] for r in refs]
    base = ixr.search_batch(queries, max_subs=2, threads=1).to_numpy()
    for t in (2, 4, 8):
        got = ixr.search_batch(queries, max_subs=2, threads=t).to_numpy()
        for key in ("query_begin", "ref_id", "offset", "n_subs", "mm_pos"):
            assert (got[key] == base[key]).all(), key


def test_a_hit_never_spans_a_record_boundary():
    ixr = TextIndex.build(["AAAA", "AAAA"], alphabet="aa", k=3)
    assert ixr.search_batch(["AAAAAAAA"], max_subs=2).num_hits == 0


def test_a_matrix_scores_hits_without_changing_which_are_returned():
    from seqtree import SubstitutionMatrix
    ixr = TextIndex.build(["CASSLGQYF"], alphabet="aa", k=4)
    plain = ixr.search_batch(["CASSIGQYF"], max_subs=1)
    scored = ixr.search_batch(["CASSIGQYF"], max_subs=1, matrix=SubstitutionMatrix.blosum62())
    assert triples(plain, 0) == triples(scored, 0)
    assert plain[0][0].score == 0
    assert scored[0][0].score == SubstitutionMatrix.blosum62().similarity("I", "L")


# --- the alphabet, which is where the silent wrong answer lives --------------------------

def test_a_query_symbol_outside_the_alphabet_names_the_query(ix):
    with pytest.raises(ValueError, match=r"queries\[1\] \('CASSUGQYF'\).*'U' is not in the"):
        ix.search_batch(["CASSLGQYF", "CASSUGQYF"], max_subs=1)


def test_a_reference_symbol_outside_the_alphabet_is_a_counted_hole():
    """The human proteome has 36 U residues and the mouse 33; refusing the build over them
    would make TextIndex unusable on the very data it exists for. They become holes instead --
    no hit crosses one -- and are reported rather than dropped silently."""
    ixr = TextIndex.build(["CASSLGQYF", "CASSUGQYF"], alphabet="aa", k=4)
    assert ixr.num_unknown == 1
    res = ixr.search_batch(["CASSLGQYF"], max_subs=1)
    assert [h.ref_id for h in res[0]] == [0]  # record 1 is unreachable through its hole


def test_ambiguity_codes_match_only_themselves():
    # X is a real symbol in this codec, not a wildcard. The proteomes carry 8,453 human and
    # 4,336 mouse non-standard residues, so this is not hypothetical.
    ixr = TextIndex.build(["XAA"], alphabet="aa", k=3)
    assert ixr.search_batch(["XAA"]).num_hits == 1
    assert ixr.search_batch(["AAA"]).num_hits == 0
    assert ixr.search_batch(["AAA"], max_subs=1).num_hits == 1


def test_a_query_shorter_than_k_is_refused_not_answered_incompletely(ix):
    with pytest.raises(ValueError, match="cannot answer a shorter query completely"):
        ix.search_batch(["CAS"])


def test_k_is_capped_so_a_direct_addressed_table_cannot_ask_for_gigabytes():
    with pytest.raises(ValueError, match="seed buckets"):
        TextIndex.build(REFS, alphabet="aa", k=8)
    assert TextIndex.build(REFS, alphabet="aa", k=5).k == 5


def test_group_ids_must_match_the_records():
    with pytest.raises(ValueError, match="group_ids has 2 entries for 4 records"):
        TextIndex.build(REFS, alphabet="aa", k=4, group_ids=[1, 2])


# --- output ------------------------------------------------------------------------------

def test_mismatch_detail_names_the_pair(ix):
    res = ix.search_batch(["CASSLGQYA"], max_subs=1)
    hit = next(h for h in res[0] if h.ref_id == 0)
    assert hit.mismatches == [(8, "A", "F")]


def test_group_fold_makes_a_tie_a_first_class_output():
    ixr = TextIndex.build(REFS, alphabet="aa", k=4, group_ids=[7, 7, 9, 11])
    res = ixr.search_batch(["CASSLGQYF"], max_subs=2, group_by=True)
    assert res.groups(0) == [(7, 0, 2), (11, 2, 1)]  # (group_id, min_subs, n_hits)


def test_groups_without_group_by_says_so(ix):
    res = ix.search_batch(["CASSLGQYF"], max_subs=1)
    with pytest.raises(ValueError, match="group_by=True"):
        res.groups(0)


def test_arrays_are_zero_copy_views_that_keep_the_result_alive(ix):
    res = ix.search_batch(["CASSLGQYF"], max_subs=2)
    view = res.arrays()["ref_id"]
    mv = memoryview(view)
    assert mv.format == "I" and mv.readonly and len(mv) == res.num_hits
    del res  # the view holds a strong reference, so this must not free the buffer
    assert list(memoryview(view)) == [0, 1, 3]


def test_to_numpy_shares_memory_and_carries_the_csr(ix):
    np = pytest.importorskip("numpy")
    res = ix.search_batch(["CASSLGQYF", "KTAYIAKQR"], max_subs=1)
    a = res.to_numpy()
    assert a["query_begin"].tolist() == [0, 2, 3]
    assert a["ref_id"].dtype == np.uint32 and a["n_subs"].dtype == np.uint16
    assert a["mm_query_aa"].tolist() == [b"F"]
    assert not a["ref_id"].flags["OWNDATA"]  # a view, not a copy


def test_indexing_is_bounds_checked_and_accepts_negatives(ix):
    res = ix.search_batch(["CASSLGQYF"], max_subs=1)
    assert res[-1] == res[0] or [tuple(h) for h in res[-1]] == [tuple(h) for h in res[0]]
    with pytest.raises(IndexError):
        res[5]


def test_an_empty_batch_is_an_empty_result(ix):
    res = ix.search_batch([], max_subs=1)
    assert len(res) == 0 and res.num_hits == 0


# --- persistence -------------------------------------------------------------------------

@pytest.mark.parametrize("mmap", [True, False])
def test_save_load_round_trips(ix, tmp_path, mmap):
    path = str(tmp_path / "text.sti")
    ix.save(path)
    back = TextIndex.load(path, mmap=mmap)
    assert (back.num_refs, back.num_residues, back.k, back.alphabet) == (
        ix.num_refs, ix.num_residues, ix.k, ix.alphabet)
    assert back.ref_seq(2) == REFS[2]
    a, b = ix.search_batch(["CASSLGQYF"], max_subs=2), back.search_batch(["CASSLGQYF"], max_subs=2)
    assert triples(a, 0) == triples(b, 0)


def test_load_rejects_a_file_that_is_not_a_text_index(tmp_path):
    path = tmp_path / "junk.sti"
    path.write_bytes(b"NOPE" + bytes(200))
    with pytest.raises(Exception):
        TextIndex.load(str(path))


# --- the search scheme, from Python -------------------------------------------------------

@pytest.mark.parametrize("k", [3, 4, 5])
def test_the_oracle_still_holds_across_every_dispatch_boundary(k):
    """`b = min(max_subs + 1, L // k)` changes at every multiple of k, and the per-block error
    budget changes with it. Sweep across those boundaries rather than sampling near one."""
    rng = random.Random(20260906 + k)
    alpha = "ACDEFG"
    refs = ["".join(rng.choice(alpha) for _ in range(rng.randint(40, 80))) for _ in range(10)]
    ixr = TextIndex.build(refs, alphabet="aa", k=k)
    queries = []
    for L in range(k, 4 * k + 3):
        r = refs[rng.randrange(len(refs))]
        off = rng.randrange(len(r) - L + 1)
        queries.append(r[off:off + L])
        mut = list(queries[-1])
        mut[rng.randrange(L)] = rng.choice(alpha)
        queries.append("".join(mut))
    for m in range(4):
        res = ixr.search_batch(queries, max_subs=m)
        for i, q in enumerate(queries):
            assert triples(res, i) == brute_force(refs, q, m), (k, q, m)


def test_k_changes_the_work_and_never_the_answer():
    rng = random.Random(1234)
    alpha = "ACDEFG"
    refs = ["".join(rng.choice(alpha) for _ in range(60)) for _ in range(8)]
    queries = [refs[rng.randrange(8)][5:5 + L] for L in (5, 6, 8, 10, 13, 19)]
    for m in range(4):
        got = [TextIndex.build(refs, alphabet="aa", k=k).search_batch(queries, max_subs=m)
               for k in (3, 4, 5)]
        for i in range(len(queries)):
            assert triples(got[0], i) == triples(got[1], i) == triples(got[2], i)


def test_a_hit_reachable_from_two_blocks_comes_back_once(ix):
    """An exact query matches every block's seed, so each block proposes the same start. Dedup
    is by block ownership rather than a sort, and a duplicated hit is the way that fails."""
    res = ix.search_batch(["MKTAYIAKQRQISFVKSHFSRQ"], max_subs=3)
    seen = [(h.ref_id, h.offset, h.n_subs) for h in res[0]]
    assert len(seen) == len(set(seen))
    assert set(seen) == brute_force(REFS, "MKTAYIAKQRQISFVKSHFSRQ", 3)


# --- limits that would otherwise be silent ------------------------------------------------

def test_a_query_too_long_for_a_16_bit_mismatch_position_is_refused(ix):
    with pytest.raises(ValueError, match="16-bit and caps a query at 65535"):
        ix.search_batch(["A" * 65536])


def test_ref_seq_reports_what_the_index_holds_not_what_was_handed_in():
    """A residue the codec never named comes back as 'X' -- and 'X' is a real symbol here, so
    the round-trip is lossy in a way that would search differently. Pinned, not discovered."""
    ixr = TextIndex.build(["CASSUGQYF"], alphabet="aa", k=4)
    assert ixr.ref_seq(0) == "CASSXGQYF"
    assert ixr.num_unknown == 1
    assert ixr.search_batch(["CASSXGQYF"]).num_hits == 0  # the hole is not an X, it is a hole


def test_an_index_file_from_an_older_format_is_refused(ix, tmp_path):
    path = tmp_path / "v1.sti"
    ix.save(str(path))
    raw = bytearray(path.read_bytes())
    raw[4] = 1  # version 2 -> 1
    path.write_bytes(bytes(raw))
    with pytest.raises(Exception, match="text index"):
        TextIndex.load(str(path))
