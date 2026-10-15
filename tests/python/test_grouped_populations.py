"""Preaggregated populations preserve duplicate full identities and missing genes."""
import numpy as np
import pytest
from seqtree import gapblock as gb

@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('exclude', [False, True])
@pytest.mark.parametrize('cuts', [[], [0, 1, 10]])
def test_population_sizes_and_exact_duplicates(paired, exclude, cuts):
    qa, qb = ['', 'A', 'AA', 'C'], ['AA', 'A', '', 'A']
    ra, rb = ['', 'A', 'A', 'AA', 'C'], ['AA', 'A', 'C', '', 'A']
    genes, rgenes = [9, 0, 1, 0], [0, 0, 1, 1, 1]
    members, exclusive = [[0, 1], [0], [1], [], [0]], [-1, 0, 1, -1, 0]
    kwargs = dict(reference_target_offsets=[0, 2, 3, 4, 4, 5],
        reference_target_ids=[0, 1, 0, 1, 0], target_peptide_ids=[0, 1],
        reference_exclusive_peptide_ids=exclusive, query_gene_ids=genes,
        reference_gene_ids=rgenes, chunk_size=2, exclude_exact=exclude, threads=1)
    batch = (gb.grouped_paired_sum_count_batch(qa, qb, ra, rb, cuts, **kwargs) if paired
        else gb.grouped_count_batch(qa, ra, cuts, **kwargs))
    chunks = list(batch)
    sizes = np.concatenate([np.asarray(c.sizes) for c in chunks])
    exact = np.concatenate([np.asarray(c.exact_counts) for c in chunks])
    for i, q in enumerate(qa):
        same = [q == r and (not paired or qb[i] == rb[j]) for j, r in enumerate(ra)]
        subsets = [list(range(len(ra)))]
        subsets += [[j for j in range(len(ra)) if t in members[j]] for t in range(2)]
        subsets += [[j for j in range(len(ra)) if t in members[j] and genes[i] == rgenes[j]] for t in range(2)]
        subsets += [[j for j in range(len(ra)) if exclusive[j] == t] for t in range(2)]
        assert sizes[i].tolist() == [sum(not (exclude and same[j]) for j in refs) for refs in subsets]
        assert exact[i].tolist() == [sum(same[j] and t in members[j] for j in range(len(ra))) for t in range(2)]
