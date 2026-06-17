import seqtree


def test_self_match_exact(idx, cdr3_db):
    p = seqtree.SearchParams()  # all caps 0 -> exact
    for i, seq in enumerate(cdr3_db):
        hits = idx.search(seq, p)
        ids = {h.ref_id for h in hits}
        assert i in ids
        me = next(h for h in hits if h.ref_id == i)
        assert me.score == 0
        assert (me.n_subs, me.n_ins, me.n_dels) == (0, 0, 0)


def test_one_substitution_finds_parent(idx):
    # CASSLELGATNEKLFF (id 1) is 2 subs from CASSLAPGATNEKLFF (id 0)
    p = seqtree.SearchParams(max_subs=2, engine="seqtm")
    hits = {h.ref_id: h for h in idx.search("CASSLAPGATNEKLFF", p)}
    assert 1 in hits
    assert hits[1].n_subs == 2
    assert hits[1].score == 2


def test_one_deletion_finds_shorter_parent(idx):
    # CASSLAPGATNEKLF (id 3) is the query with one extra char vs nothing... query is full,
    # ref id 3 is one char shorter -> reachable by a single deletion from the query.
    p = seqtree.SearchParams(max_dels=1, engine="seqtm")
    hits = {h.ref_id: h for h in idx.search("CASSLAPGATNEKLFF", p)}
    assert 3 in hits
    assert hits[3].n_dels == 1


def test_tuple_unpacking(idx):
    p = seqtree.SearchParams()
    hit = idx.search("CASSPQGATNEKLFF", p)[0]
    ref_id, score, n_subs, n_ins, n_dels = hit
    assert (ref_id, score, n_subs, n_ins, n_dels) == (hit.ref_id, hit.score, 0, 0, 0)


def test_top_hit(idx):
    p = seqtree.SearchParams(max_subs=3, engine="seqtm")
    top = idx.search_top("CASSLAPGATNEKLFF", p, k=3)
    assert top[0].ref_id == 0  # exact match ranks first
    assert top[0].score == 0
    scores = [h.score for h in top]
    assert scores == sorted(scores)


def test_seqtrie_total_edits_only(idx):
    # seqtrie ignores per-type caps; max_total_edits alone drives the budget.
    p = seqtree.SearchParams(max_total_edits=1, engine="seqtrie")
    hits = {h.ref_id: h for h in idx.search("CASSLWPGATNEKLFF", p)}
    assert 0 in hits
    assert hits[0].score == 1


def test_matrix_budget(idx):
    p = seqtree.SearchParams(matrix="BLOSUM62", max_penalty=5, engine="seqtrie", gap_open=100)
    hits = {h.ref_id: h for h in idx.search("CASSLAPGATNEKLFF", p)}
    assert 0 in hits and hits[0].score == 0


def test_alignment_consistency(idx):
    p = seqtree.SearchParams(max_subs=2, engine="seqtm")
    hit = {h.ref_id: h for h in idx.search("CASSLAPGATNEKLFF", p)}[1]
    aln = idx.align(1, "CASSLAPGATNEKLFF", p)
    assert aln.ops.count("S") == hit.n_subs
    assert aln.ops.count("I") == hit.n_ins
    assert aln.ops.count("D") == hit.n_dels
    assert aln.score == hit.score


def test_nucleotide_alphabet():
    idx = seqtree.Index.build(["ACGTACGT", "ACGAACGT", "ACGTACGA"], alphabet="nt")
    p = seqtree.SearchParams(max_subs=1, engine="seqtm")
    hits = {h.ref_id for h in idx.search("ACGTACGT", p)}
    assert hits == {0, 1, 2}
