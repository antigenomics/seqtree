"""Compact cutoff bins match exact dense scores, including integer boundaries."""
import numpy as np
import pytest
from seqtree import gapblock as gb

@pytest.mark.parametrize('cuts', [[], [-2147483648], [-9, -3], [-2147483648, 0],
    [0], [0, 4, 8, 16, 32], [-8, -4, 0, 4, 8], [0, 2147483647],
    [-2147483648, -1, 0, 2, 2147483647], [-2147483648, 2147483647]])
@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('quantum', [1, 3])
@pytest.mark.parametrize('exclude', [False, True])
def test_cutoff_lookup_matches_full_score_oracle(cuts, paired, quantum, exclude):
    q = ['', 'A', 'AC', 'AAC', 'AAAA', 'CCCCCC']
    r = ['', 'A', 'CA', 'AAC', 'CCCC', 'CCCCCC']
    raw = np.asarray(gb.score_matrix(q, r, gap_open=2, gap_extend=1))
    expected = raw // quantum * quantum + 5
    lane = dict(query_group_ids=[0]*len(q), reference_group_ids=[0]*len(r),
        group_distances=[[3]], score_quantum=quantum,
        additional_feature_offsets=[gb.PairedFeatureOffsets([[0, 0]]*len(q),
            [[0, 0]]*len(r), [[1]], scale=2)])
    options = dict(gap_open=2, gap_extend=1)
    if paired:
        expected = expected + expected[:, ::-1]
        options.update({key+'_'+suffix: value for suffix in ('alpha', 'beta') for key, value in lane.items()})
    else:
        options.update(lane)
    metadata = dict(reference_target_offsets=list(range(len(r)+1)), reference_target_ids=[0]*len(r),
        target_peptide_ids=[0], reference_exclusive_peptide_ids=[0]*len(r),
        query_gene_ids=[0]*len(q), reference_gene_ids=[0]*len(r),
        threads=1, exclude_exact=exclude, chunk_size=2)
    batch = (gb.grouped_paired_sum_count_batch(q, q, r, r[::-1], cuts, **options, **metadata) if paired
        else gb.grouped_count_batch(q, r, cuts, **options, **metadata))
    actual = np.concatenate([np.asarray(c.counts)[:, 0] for c in batch])
    eligible = np.asarray([[not (exclude and a == b and (not paired or a == r[-1-j]))
        for j, b in enumerate(r)] for a in q])
    oracle = np.stack([((expected <= cut) & eligible).sum(axis=1) for cut in cuts], axis=1) if cuts else np.empty((len(q), 0), dtype=np.uint32)
    np.testing.assert_array_equal(actual, oracle)
