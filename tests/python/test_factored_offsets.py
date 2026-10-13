"""Factored offsets preserve the explicit code-pair penalty oracle without Q*N storage."""
import math
import pytest
from seqtree.gapblock import (
    PairedFeatureOffsets, score_matrix, count_batch, ball_batch,
    group_audit_batch, paired_sum_count_batch,
)


@pytest.mark.parametrize("threads", [1, 4])
def test_factored_offsets_explicit_oracle(threads, require_thread_budget):
    q = ['AAC', 'AC', 'AAA', 'AAC']
    r = ['AAC', 'ACC', 'AC', 'AAC', 'CCC']
    qc = [[0, 1, 2, 0], [1, 2, 0, 1], [2, 0, 1, 2], [0, 1, 2, 0]]
    rc = [[0, 1, 2, 0], [2, 1, 0, 2], [1, 2, 0, 1], [0, 1, 2, 0], [2, 2, 2, 2]]
    costs = [[0, 1, 3], [2, 0, 5], [4, 3, 0]]
    v = [[0, 3], [7, 1]]
    groups = dict(query_group_ids=[0, 1, 0, 0], reference_group_ids=[0, 1, 1, 0, 1], group_distances=v)
    thresholds = [[-1, 0, 1, 5, 12, 30, 12], [], [30, 2, 7], [4, 13]]
    component = score_matrix(q, r, gap_open=2, threads=1)
    for weight in [0.0, 0.1, 0.5, 1.0, 1.25, 2.5]:
        offsets = PairedFeatureOffsets(qc, rc, costs, weight=weight, scale=2)
        full = [[component[i, j] + v[groups['query_group_ids'][i]][groups['reference_group_ids'][j]] +
                 2 * round(weight * sum(max(costs[a][b], costs[c][d])
                     for a, c, b, d in zip(qc[i][::2], qc[i][1::2], rc[j][::2], rc[j][1::2])) / 2)
                 for j in range(len(r))] for i in range(len(q))]
        args = dict(gap_open=2, threads=threads, feature_offsets=offsets, **groups)
        matrix = score_matrix(q, r, **args)
        assert [[matrix[i, j] for j in range(len(r))] for i in range(len(q))] == full
        for excluded in [False, True]:
            counts, mass = count_batch(q, r, thresholds, linear_mass=True, exclude_exact=excluded, **args)
            expected = [[sum(s <= t for j, s in enumerate(row) if not excluded or q[i] != r[j]) for t in thresholds[i]] for i, row in enumerate(full)]
            expected_mass = [[sum(max(0, t-s) for j, s in enumerate(row) if not excluded or q[i] != r[j]) for t in thresholds[i]] for i, row in enumerate(full)]
            assert counts == expected
            assert mass == expected_mass
            assert count_batch(q, r, thresholds, exclude_exact=excluded, **args) == expected
            cuts = [10, 6, 15, -1]
            balls = ball_batch(q, r, cuts, exclude_exact=excluded, **args)
            assert [[(h.ref_id, h.score) for h in row] for row in balls] == [
                [(j, s) for j, s in enumerate(row) if s <= cuts[i] and (not excluded or q[i] != r[j])]
                for i, row in enumerate(full)]
            before, after, nearest = group_audit_batch(q, r, cuts, exclude_exact=excluded, **args)
            for i in range(len(q)):
                eligible = [j for j in range(len(r)) if component[i, j] <= cuts[i] and (not excluded or q[i] != r[j])]
                assert before[i] == len(eligible)
                assert after[i] == sum(full[i][j] <= cuts[i] for j in eligible)
                assert nearest[i] == min((full[i][j] for j in eligible), default=-1)
            pair_counts, pair_mass = paired_sum_count_batch(q, q, r, r, thresholds,
                gap_open=2, threads=threads, exclude_exact=excluded,
                feature_offsets_alpha=offsets, feature_offsets_beta=offsets)
            pair_scores = [[2 * (s - v[groups['query_group_ids'][i]][groups['reference_group_ids'][j]]) for j, s in enumerate(row)] for i, row in enumerate(full)]
            assert pair_counts == [[sum(s <= t for j, s in enumerate(row) if not excluded or q[i] != r[j]) for t in thresholds[i]] for i, row in enumerate(pair_scores)]
            assert pair_mass == [[sum(max(0,t-s) for j,s in enumerate(row) if not excluded or q[i] != r[j]) for t in thresholds[i]] for i,row in enumerate(pair_scores)]


def test_factored_rounding_and_owned_inputs():
    # 0.5 -> 0, 1.5 -> 2, 2.5 -> 2; decimal arithmetic follows multiply/divide.
    costs = [[0, 1, 3, 5, 25]]
    codes = [[1, 1], [2, 2], [3, 3], [4, 4]]
    f = PairedFeatureOffsets([[0, 0]], codes, costs, weight=0.5, scale=200)
    costs[0][1] = 1000
    codes[0][0] = 4
    sm = score_matrix(['A'], ['A']*4, feature_offsets=f, threads=1)
    assert [sm[0,j] for j in range(4)] == [0, 400, 400, 2400]
    f = PairedFeatureOffsets([[0,0]], [[1,1]], [[0,25]], weight=0.1)
    assert score_matrix(['A'], ['A'], feature_offsets=f, threads=1)[0,0] == 2
    assert count_batch([], [], [], feature_offsets=PairedFeatureOffsets([], [], [[0]]), threads=1) == []
    assert count_batch(['A'], [], [[1]], feature_offsets=PairedFeatureOffsets([[0,0]], [], [[0]]), threads=1) == [[0]]


@pytest.mark.parametrize('qc,rc,costs,weight,scale', [
    ([], [[0,0]], [[0]], 1, 1), ([[0]], [[0]], [[0]], 1, 1),
    ([[]], [[]], [[0]], 1, 1), ([[0,0]], [[0,0,0,0]], [[0]], 1, 1),
    ([[-1,0]], [[0,0]], [[0]], 1, 1), ([[1,0]], [[0,0]], [[0]], 1, 1),
    ([[0,0]], [[1,0]], [[0]], 1, 1), ([[0,0]], [[0,0]], [], 1, 1),
    ([[0,0]], [[0,0]], [[0], []], 1, 1), ([[0,0]], [[0,0]], [[-1]], 1, 1),
    ([[0,0]], [[0,0]], [[0]], -1, 1), ([[0,0]], [[0,0]], [[0]], math.inf, 1),
    ([[0,0]], [[0,0]], [[0]], math.nan, 1), ([[0,0]], [[0,0]], [[0]], 1, -1),
    ([[0,0]], [[0,0]], [[1]], 1e308, 1), ([[0,0]], [[0,0]], [[2**30]], 1, 2),
])
def test_factored_validation(qc,rc,costs,weight,scale):
    f = PairedFeatureOffsets(qc,rc,costs,weight=weight,scale=scale)
    with pytest.raises(ValueError):
        count_batch(['A'], ['A'], [[10]], feature_offsets=f, threads=1)


def test_factored_combined_overflow():
    f = PairedFeatureOffsets([[0,0]], [[0,0]], [[2**30]])
    with pytest.raises(ValueError, match='overflow'):
        score_matrix(['A'], ['A'], feature_offsets=f, query_group_ids=[0],
                     reference_group_ids=[0], group_distances=[[2**30]], threads=1)


def test_factored_numpy_rounding_oracle():
    import numpy as np
    for depth in [1, 2, 3, 5, 9]:
        costs = [list(range(100))]
        qcodes = [[0,0] * depth]
        rcodes = [[n,n] * depth for n in range(100)]
        for weight in [0.1, 0.2, 0.5, 0.7, 1.1, 1.25, 3.141592653589793]:
            descriptor = PairedFeatureOffsets(qcodes, rcodes, costs, weight=weight, scale=200)
            scores = score_matrix(['A'], ['A']*100, threads=1, feature_offsets=descriptor)
            expected = 200 * np.rint(weight * (np.arange(100, dtype=np.int64) * depth) / depth)
            assert [scores[0,j] for j in range(100)] == expected.tolist()


def test_factored_paired_identity_and_independent_lanes():
    a = PairedFeatureOffsets([[0,0]], [[0,0], [1,1]], [[0,3], [2,0]], scale=2)
    b = PairedFeatureOffsets([[0,0]], [[1,1], [0,0]], [[0,5], [4,0]], weight=0.5)
    args = dict(threads=1, feature_offsets_alpha=a, feature_offsets_beta=b,
                query_group_ids_alpha=[0], reference_group_ids_alpha=[0,0], group_distances_alpha=[[1]],
                query_group_ids_beta=[0], reference_group_ids_beta=[0,0], group_distances_beta=[[2]])
    assert paired_sum_count_batch(['A'], ['A'], ['A','A'], ['C','A'], [[3,4,6,9]], **args) == ([[0,0,1,2]], [[0,0,0,3]])
    assert paired_sum_count_batch(['A'], ['A'], ['A','A'], ['C','A'], [[9]], exclude_exact=True, **args) == ([[1]], [[3]])


def test_factored_zero_controls():
    for weight, scale in [(0, 200), (1e308, 0)]:
        f = PairedFeatureOffsets([[0,0]], [[0,0]], [[2**30]], weight=weight, scale=scale)
        assert score_matrix(['A'], ['C'], feature_offsets=f, threads=1)[0,0] == 1
