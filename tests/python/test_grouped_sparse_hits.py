"""Sparse native affine hits against an independent short-sequence DP oracle."""
from types import SimpleNamespace

import pytest
from seqtree import SubstitutionMatrix, gapblock

Q = ['WAAAAG', 'WRRAAG', 'WAGGAG', 'WAAAAG']
R = ['WAAAAG', 'WRRAAG', 'WAGGAG', 'WARRAAG', 'WAAAAG']
LABELS = [[0], [0, 1], [1], [1], [0]]
TARGETS = [0, 0, 1, 1]
MATRIX = SubstitutionMatrix.blosum62()


def options():
    return dict(reference_target_offsets=[0, 1, 3, 4, 5, 6],
        reference_target_ids=[0, 0, 1, 1, 1, 0], target_peptide_ids=[0, 1],
        reference_exclusive_peptide_ids=[0, -1, 1, 1, 0],
        query_gene_ids=[0]*len(Q), reference_gene_ids=[0]*len(R),
        matrix=MATRIX, gap_open=5, gap_extend=2, alignment_scale=3,
        trim_left=1, trim_right=1, query_group_ids=[0, 1, 0, 1],
        reference_group_ids=[1, 0, 1, 0, 0], group_distances=[[0, 7], [9, 0]],
        feature_offsets=gapblock.PairedFeatureOffsets([[0, 0]]*len(Q),
            [[1, 1]]*len(R), [[0, 2], [2, 0]], scale=4),
        additional_feature_offsets=[gapblock.PairedFeatureOffsets([[0, 0]]*len(Q),
            [[1, 1]]*len(R), [[0, 3], [3, 0]], scale=2)])


def similarity(q, r, mode):
    # Explicit three-state affine DP, independent of seqtree's native scorer.
    h = [[0]*(len(r)+1) for _ in range(len(q)+1)]
    e = [[-10**9]*(len(r)+1) for _ in range(len(q)+1)]
    f = [[-10**9]*(len(r)+1) for _ in range(len(q)+1)]
    if mode == 'global':
        for i in range(1, len(q)+1): h[i][0] = -5-2*(i-1)
        for j in range(1, len(r)+1): h[0][j] = -5-2*(j-1)
    for i, a in enumerate(q, 1):
        for j, b in enumerate(r, 1):
            e[i][j] = max(h[i-1][j]-5, e[i-1][j]-2)
            f[i][j] = max(h[i][j-1]-5, f[i][j-1]-2)
            h[i][j] = max(h[i-1][j-1]+MATRIX.similarity(a, b), e[i][j], f[i][j])
            if mode == 'local': h[i][j] = max(0, h[i][j])
    return max(map(max, h)) if mode == 'local' else h[-1][-1]


def expected(mode, exclude, self_ids, quantum):
    expected = []
    for i, fullq in enumerate(Q):
        for j, fullr in enumerate(R):
            if TARGETS[i] not in LABELS[j]: continue
            if self_ids is not None and j <= self_ids[i]: continue
            if exclude and fullq == fullr: continue
            q, r = fullq[1:-1], fullr[1:-1]
            diagonal = sum(MATRIX.similarity(x, x) for x in q+r)
            raw = 3*max(0, diagonal-2*similarity(q, r, mode))
            offset = [[0, 7], [9, 0]][[0, 1, 0, 1][i]][[1, 0, 1, 0, 0][j]]+14
            score = raw//quantum*quantum+offset
            if score <= 120: expected.append([i, j, score])
    return expected


@pytest.mark.parametrize('threads', [1, 2])
@pytest.mark.parametrize('mode', ['local', 'global'])
@pytest.mark.parametrize('exclude', [False, True])
@pytest.mark.parametrize('self_ids', [None, [0, 1, 2, -1]])
@pytest.mark.parametrize('quantum', [1, 7])
def test_sparse_hits_oracle_filtering_threads_and_chunks(mode, exclude, self_ids, quantum, threads, require_thread_budget):
    want = expected(mode, exclude, self_ids, quantum)
    for chunk_size in [1, 3, 4]:
        chunks = list(gapblock.grouped_count_batch(Q, R, [0, 30, 120], **options(),
            alignment_mode=mode, exclude_exact=exclude, return_hits=True,
            query_target_ids=TARGETS, query_reference_ids=self_ids,
            threads=threads, chunk_size=chunk_size, score_quantum=quantum))
        assert [hit for chunk in chunks for hit in memoryview(chunk.hits).tolist()] == want
        for chunk in chunks:
            assert memoryview(chunk.hits).readonly
            assert memoryview(chunk.hits).shape[1] == 3
            assert memoryview(chunk.counts).nbytes == 0
            assert memoryview(chunk.sizes).nbytes == 0
            assert memoryview(chunk.exact_counts).nbytes == 0
            assert all(chunk.start <= row[0] < chunk.stop for row in memoryview(chunk.hits).tolist())


def test_sparse_cap_is_actual_output_and_never_silent_truncation():
    args = dict(reference_target_offsets=[0, 1, 2, 3], reference_target_ids=[0]*3,
        target_peptide_ids=[0], reference_exclusive_peptide_ids=[0]*3,
        query_gene_ids=[0], reference_gene_ids=[0]*3, query_target_ids=[0],
        return_hits=True, max_hit_bytes=12, threads=1)
    # Three candidates but only one hit: a worst-case candidate preallocation would fail.
    chunk = next(gapblock.grouped_count_batch(['AA'], ['AA', 'RR', 'WW'], [0], **args))
    assert memoryview(chunk.hits).tolist() == [[0, 0, 0]]
    with pytest.raises(OverflowError, match='max_hit_bytes'):
        next(gapblock.grouped_count_batch(['AA'], ['AA']*3, [0], **args))


@pytest.mark.parametrize('threads', [1, 2])
def test_sparse_cap_is_shared_across_workers(threads, require_thread_budget):
    with pytest.raises(OverflowError, match='max_hit_bytes'):
        next(gapblock.grouped_count_batch(['AA']*4, ['AA'], [0],
            reference_target_offsets=[0, 1], reference_target_ids=[0],
            target_peptide_ids=[0], reference_exclusive_peptide_ids=[0],
            query_gene_ids=[0]*4, reference_gene_ids=[0], query_target_ids=[0]*4,
            return_hits=True, max_hit_bytes=36, threads=threads))


@pytest.mark.parametrize('changes', [
    {'query_target_ids': None}, {'query_target_ids': [0]}, {'query_target_ids': [0, 0, 0, 2]},
    {'query_reference_ids': [0]}, {'query_reference_ids': [-2]*4},
    {'query_reference_ids': [5]*4}, {'max_hit_bytes': 11}, {'max_hit_bytes': True},
    {'return_hits': False},
])
def test_sparse_input_validation(changes):
    args = dict(return_hits=True, query_target_ids=TARGETS, threads=1, alignment_mode='local')
    args.update(changes)
    with pytest.raises(ValueError):
        gapblock.grouped_count_batch(Q, R, [120], **options(), **args)


def test_native_sparse_capability_guard(monkeypatch):
    monkeypatch.setattr(gapblock, '_GroupedGeometry', SimpleNamespace)
    with pytest.raises(RuntimeError, match='native return_hits'):
        gapblock.grouped_count_batch(['AA'], ['AA'], [0],
            reference_target_offsets=[0, 1], reference_target_ids=[0], target_peptide_ids=[0],
            reference_exclusive_peptide_ids=[0], query_gene_ids=[0], reference_gene_ids=[0],
            return_hits=True, query_target_ids=[0], threads=1)


@pytest.mark.parametrize('cuts', [[], [-1]])
def test_sparse_empty_hits_keep_h_by_three_shape(cuts):
    chunk = next(gapblock.grouped_count_batch(Q, R, cuts, **options(),
        alignment_mode='local', return_hits=True, query_target_ids=TARGETS, threads=1))
    assert memoryview(chunk.hits).shape == (0, 3)


def test_default_counts_stay_available_and_hits_empty():
    chunk = next(gapblock.grouped_count_batch(Q, R, [120], **options(),
        alignment_mode='local', threads=1))
    assert memoryview(chunk.counts).shape == (4, 7, 1)
    assert memoryview(chunk.hits).shape == (0, 3)


def test_sparse_empty_axes():
    args = dict(reference_target_offsets=[0], reference_target_ids=[],
        target_peptide_ids=[], reference_exclusive_peptide_ids=[], reference_gene_ids=[],
        return_hits=True, threads=1)
    assert list(gapblock.grouped_count_batch([], [], [0], query_gene_ids=[],
        query_target_ids=[], **args)) == []
    args.update(target_peptide_ids=[0])
    chunk = next(gapblock.grouped_count_batch(['AA'], [], [0], query_gene_ids=[0],
        query_target_ids=[0], **args))
    assert memoryview(chunk.hits).shape == (0, 3)


def test_same_sequence_distinct_reference_records_retain_zero_distance_edges():
    # Reference index carries caller identity; sequence equality does not merge rows.
    args = dict(reference_target_offsets=[0, 1, 2, 3], reference_target_ids=[0]*3,
        target_peptide_ids=[0], reference_exclusive_peptide_ids=[0]*3,
        query_gene_ids=[0]*3, reference_gene_ids=[0, 1, 2], return_hits=True,
        query_target_ids=[0]*3, query_reference_ids=[0, 1, 2], threads=1,
        matrix=MATRIX, gap_open=5, gap_extend=2, alignment_mode='local')
    chunks = list(gapblock.grouped_count_batch(['AAAA']*3, ['AAAA']*3, [0], **args))
    assert [hit for chunk in chunks for hit in memoryview(chunk.hits).tolist()] == [
        [0, 1, 0], [0, 2, 0], [1, 2, 0]]
