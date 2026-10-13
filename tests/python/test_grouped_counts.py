"""Generic grouped counts equal explicit existing count populations, without domain data."""
import gc
import numpy as np
import pytest
from seqtree import gapblock as gb

Q = ['AAA', 'AAC', 'TTT', 'AAA']
R = ['AAA', 'ACA', 'AAA', 'CCCC', 'TTT']
QB = ['CCC', 'CCA', 'AAA', 'CCC']
RB = ['CCC', 'CCC', 'CCA', 'CCCA', 'AAA']
CUTS = [-1, 0, 2, 5, 20]
MEMBERS = [[0, 1], [2], [0, 2], [], []]
PARENTS = [0, 0, 1]
EXCLUSIVE = [0, 1, -1, 0, -1]
QGENE, RGENE = [0, 1, 1, 0], [1, 0, 1, 0, 1]


def metadata():
    return dict(reference_target_offsets=[0, 2, 3, 5, 5, 5],
        reference_target_ids=[0, 1, 2, 0, 2], target_peptide_ids=PARENTS,
        reference_exclusive_peptide_ids=EXCLUSIVE, query_gene_ids=QGENE,
        reference_gene_ids=RGENE, chunk_size=2)


def geometry(features, indices=range(5), lane=0):
    indices = list(indices)
    qg, rg = [0, 1, 0, 1], [1, 0, 0, 1, 1]
    result = dict(gap_open=2, gap_extend=1, query_group_ids=qg,
                  reference_group_ids=[rg[i] for i in indices], group_distances=[[0, 3], [3, 2]])
    if features:
        qc = [[0, 1], [1, 0], [1, 1], [0, 1]]
        rc = [[1, 0], [0, 0], [0, 1], [1, 1], [1, 0]]
        result['feature_offsets'] = gb.PairedFeatureOffsets(qc, [rc[i] for i in indices],
            [[0, 3], [3, 0]], weight=.5 if lane == 0 else 1.5, scale=1)
    return result


def paired_geometry(features, indices=range(5)):
    result = {}
    for lane, suffix in enumerate(('alpha', 'beta')):
        for key, value in geometry(features, indices, lane).items():
            result[key if key in ('gap_open', 'gap_extend') else key+'_'+suffix] = value
    return result


@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('exclude', [False, True])
@pytest.mark.parametrize('features', [False, True])
@pytest.mark.parametrize('threads', [1, 4])
def test_every_grouped_curve_and_population_matches_original(paired, exclude, features, threads):
    options = paired_geometry(features) if paired else geometry(features)
    batch = (gb.grouped_paired_sum_count_batch(Q, QB, R, RB, CUTS, **metadata(), **options,
               threads=threads, exclude_exact=exclude) if paired else
             gb.grouped_count_batch(Q, R, CUTS, **metadata(), **options, threads=threads, exclude_exact=exclude))
    assert (batch.query_count, batch.reference_count, batch.target_count, batch.peptide_count) == (4, 5, 3, 2)
    chunks = list(batch)
    assert [(c.start, c.stop) for c in chunks] == [(0, 2), (2, 4)]
    counts = np.concatenate([np.asarray(c.counts) for c in chunks])
    sizes = np.concatenate([np.asarray(c.sizes) for c in chunks])
    exact = np.concatenate([np.asarray(c.exact_counts) for c in chunks])
    assert counts.dtype == sizes.dtype == exact.dtype == np.uint32
    assert not np.asarray(chunks[0].counts).flags.writeable
    for i, q in enumerate(Q):
        identical = [q == r and (not paired or QB[i] == RB[j]) for j, r in enumerate(R)]
        subsets = [list(range(5))]
        subsets += [[j for j in range(5) if t in MEMBERS[j]] for t in range(3)]
        subsets += [[j for j in range(5) if t in MEMBERS[j] and QGENE[i] == RGENE[j]] for t in range(3)]
        subsets += [[j for j in range(5) if EXCLUSIVE[j] == p] for p in range(2)]
        for channel, refs in enumerate(subsets):
            subopt = paired_geometry(features, refs) if paired else geometry(features, refs)
            native = (gb.paired_sum_count_batch(Q, QB, [R[j] for j in refs], [RB[j] for j in refs],
                CUTS, **subopt, shared_thresholds=True, exclude_exact=exclude, threads=threads)[0] if paired else
                gb.count_batch(Q, [R[j] for j in refs], CUTS, **subopt,
                    shared_thresholds=True, exclude_exact=exclude, threads=threads))
            np.testing.assert_array_equal(counts[i, channel], native[i])
            assert sizes[i, channel] == sum(not (exclude and identical[j]) for j in refs)
        for t in range(3):
            assert exact[i, t] == sum(identical[j] and t in MEMBERS[j] for j in range(5))
        # All-source complement retains shared-label keys and unselected references.
        for p in range(2):
            refs = [j for j in range(5) if EXCLUSIVE[j] != p]
            subopt = paired_geometry(features, refs) if paired else geometry(features, refs)
            native = (gb.paired_sum_count_batch(Q, QB, [R[j] for j in refs], [RB[j] for j in refs],
                CUTS, **subopt, shared_thresholds=True, exclude_exact=exclude)[0] if paired else
                gb.count_batch(Q, [R[j] for j in refs], CUTS, **subopt,
                    shared_thresholds=True, exclude_exact=exclude))
            np.testing.assert_array_equal(counts[i, 0]-counts[i, 7+p], native[i])
    held = np.asarray(chunks[0].counts)
    del chunks, batch
    gc.collect()
    assert held.shape == (2, 9, 5) and not held.flags.writeable


@pytest.mark.parametrize('change', [dict(reference_target_offsets=[]),
    dict(reference_target_ids=[0, 0, 2, 0, 2]), dict(reference_target_ids=[0, 9, 2, 0, 2]),
    dict(reference_exclusive_peptide_ids=[1, 1, -1, 0, -1]),
    dict(reference_exclusive_peptide_ids=[0, 1, -2, 0, -1]),
    dict(query_gene_ids=[0]), dict(reference_gene_ids=[-1]*5), dict(chunk_size=0)])
def test_invalid_group_metadata_rejected_before_iteration(change):
    args = metadata(); args.update(change)
    with pytest.raises((ValueError, TypeError, OverflowError)):
        gb.grouped_count_batch(Q, R, CUTS, **args)


@pytest.mark.parametrize('cuts', [[1, 0], [0, 0], [True]])
def test_invalid_cutoffs_rejected(cuts):
    with pytest.raises((ValueError, TypeError)):
        gb.grouped_count_batch(Q, R, cuts, **metadata())


def test_empty_query_reference_and_threshold_axes():
    empty = dict(reference_target_offsets=[0], reference_target_ids=[], target_peptide_ids=[],
        reference_exclusive_peptide_ids=[], reference_gene_ids=[], query_gene_ids=[])
    assert list(gb.grouped_count_batch([], [], [], **empty)) == []
    empty['query_gene_ids'] = [0]
    chunk = next(gb.grouped_count_batch(['AAA'], [], [], **empty))
    assert np.asarray(chunk.counts).shape == (1, 1, 0)
    assert np.asarray(chunk.sizes).tolist() == [[0]]
    assert np.asarray(chunk.exact_counts).shape == (1, 0)


def test_prepared_ownership_and_weighted_prior_geometry():
    queries, refs = ['AAA', 'AAAA'], ['AAA', 'AAAA', 'ACA']
    args = dict(reference_target_offsets=[0, 1, 2, 3], reference_target_ids=[0]*3,
        target_peptide_ids=[0], reference_exclusive_peptide_ids=[0]*3,
        query_gene_ids=[0, 1], reference_gene_ids=[1, 0, 1], chunk_size=1)
    options = dict(query_position_weights=[[1, 2, 3], [3, 2, 1, 3]],
        reference_position_weights=[[1, 1, 2], [2, 1, 3, 1], [2, 2, 2]],
        gap_prior=gb.central_prior(2), gap_open=2, gap_extend=1)
    expected = gb.count_batch(queries, refs, CUTS, shared_thresholds=True, **options)
    batch = gb.grouped_count_batch(queries, refs, CUTS, **args, **options)
    queries.append('CCC'); refs.clear(); args['reference_target_ids'][0] = 99
    options['query_position_weights'][0][0] = 99
    actual = np.concatenate([np.asarray(chunk.counts)[:, 0] for chunk in batch])
    np.testing.assert_array_equal(actual, expected)
