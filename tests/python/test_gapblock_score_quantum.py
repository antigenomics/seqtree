"""Generic alignment quantization before additive offsets, no domain inputs."""
import pytest
from seqtree import gapblock, SubstitutionMatrix


def metadata(nq, nr):
    return dict(reference_target_offsets=list(range(nr + 1)),
        reference_target_ids=[0] * nr, target_peptide_ids=[0],
        reference_exclusive_peptide_ids=[0] * nr,
        query_gene_ids=[0] * nq, reference_gene_ids=[0] * nr)


def rows(buffer):
    # Native buffer interface without a runtime NumPy dependency.
    return memoryview(buffer).tolist()


@pytest.mark.parametrize('threads', [1, 4])
@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('factored', [False, True])
@pytest.mark.parametrize('withheld', [False, True])
def test_floor_each_alignment_before_offsets_and_pair_sum(threads, paired, factored, withheld, require_thread_budget):
    q, r = ['AAAA', 'RRRA'], ['ARAA', 'RRRA', 'AAAA', 'AA']
    cuts = [0, 3, 5, 6, 9, 20, 100]
    matrix = SubstitutionMatrix.blosum62()
    weights = {2: [3, 3], 4: [3, 3, 3, 3]}
    raw = gapblock.score_matrix(q, r, matrix=matrix, gap_open=7, gap_extend=2,
        position_weights_by_length=weights, threads=1)
    opts = dict(matrix=matrix, gap_open=7, gap_extend=2, threads=threads,
        exclude_exact=True, **metadata(len(q), len(r)))
    geometry = dict(position_weights_by_length=weights,
        query_group_ids=[0, 1], reference_group_ids=[1, 0, 0, 1],
        group_distances=[[0, 3], [5, 0]])
    offset = [[3, 0, 0, 3], [0, 5, 5, 0]]
    if factored:
        geometry['feature_offsets'] = gapblock.PairedFeatureOffsets(
            [[0, 0]] * len(q), [[1, 1]] * len(r), [[0, 2], [2, 0]])
        offset = [[value + 2 for value in row] for row in offset]
    if withheld:
        opts.update(query_withheld_peptide_ids=[-1] * len(q),
            reference_peptide_offsets=list(range(len(r) + 1)), reference_peptide_ids=[0] * len(r))
    if paired:
        kwargs = {key + suffix: value for suffix in ('_alpha', '_beta')
                  for key, value in geometry.items()}
        iterator = gapblock.grouped_paired_sum_count_batch(q, q, r, r, cuts,
            score_quantum_alpha=5, score_quantum_beta=7, **opts, **kwargs)
    else:
        iterator = gapblock.grouped_count_batch(q, r, cuts, score_quantum=5,
            **opts, **geometry)
    chunk = next(iterator)
    expected = []
    for i, query in enumerate(q):
        scores = [(raw[i, j] // 5) * 5 + offset[i][j] for j in range(len(r))]
        if paired:
            scores = [value + (raw[i, j] // 7) * 7 + offset[i][j]
                      for j, value in enumerate(scores)]
        curve = [sum(score <= cutoff and query != r[j]
                     for j, score in enumerate(scores)) for cutoff in cuts]
        expected.append([curve] * 4)
    assert rows(chunk.counts) == expected
    assert rows(chunk.sizes) == [[3] * 4, [3] * 4]
    assert rows(chunk.exact_counts) == [[1], [1]]


@pytest.mark.parametrize('bad', [0, -1, True, 1.5, 2**31])
def test_invalid_quantum_rejected_even_empty_axes(bad):
    with pytest.raises((ValueError, OverflowError), match='quantum'):
        gapblock.grouped_count_batch([], [], [], score_quantum=bad, **metadata(0, 0))


def test_default_quantum_one_is_exactly_original_and_empty_axes():
    q, r = ['AAAA'], ['ARAA', 'AAAA']
    kwargs = metadata(1, 2)
    a = next(gapblock.grouped_count_batch(q, r, [0, 1, 2, 3], **kwargs))
    b = next(gapblock.grouped_count_batch(q, r, [0, 1, 2, 3], score_quantum=1, **kwargs))
    assert rows(a.counts) == rows(b.counts)
    assert rows(a.sizes) == rows(b.sizes)
    assert list(gapblock.grouped_count_batch([], [], [], score_quantum=600, **metadata(0, 0))) == []


def test_native_quantum_validation_is_not_only_python(monkeypatch):
    original = gapblock._grouped_geometry
    def invalid(*args, **kwargs):
        result = original(*args, **kwargs)
        result.score_quantum = 0
        return result
    monkeypatch.setattr(gapblock, '_grouped_geometry', invalid)
    with pytest.raises(ValueError, match='quantum'):
        gapblock.grouped_count_batch([], [], [], **metadata(0, 0))


def test_max_int32_quantum_and_empty_reference_cutoffs():
    chunk = next(gapblock.grouped_count_batch(['AAAA'], ['RRRR'], [0, 1],
        score_quantum=2**31 - 1, **metadata(1, 1)))
    assert rows(chunk.counts) == [[[1, 1]] * 4]
    chunk = next(gapblock.grouped_count_batch(['AAAA'], [], [],
        score_quantum=600, **metadata(1, 0)))
    assert rows(chunk.sizes) == [[0] * 4]
