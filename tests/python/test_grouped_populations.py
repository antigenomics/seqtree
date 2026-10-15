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
        assert sizes[i].tolist() == [len(refs) for refs in subsets]
        assert exact[i].tolist() == [sum(same[j] and t in members[j] for j in range(len(ra))) for t in range(2)]


@pytest.mark.parametrize('paired',[False,True])
@pytest.mark.parametrize('withheld',[False,True])
def test_exact_hit_mask_keeps_fixed_model_populations(paired,withheld):
    from seqtree import _core
    assert _core.exact_exclusion_hits_only is True
    options=dict(reference_target_offsets=[0,1,2,3],reference_target_ids=[0,0,1],
        target_peptide_ids=[0,1],reference_exclusive_peptide_ids=[0,0,1],
        query_gene_ids=[0],reference_gene_ids=[0,0,0],threads=1)
    if withheld:options.update(query_withheld_peptide_ids=[1],reference_peptide_offsets=[0,1,2,3],reference_peptide_ids=[0,0,1])
    def run(exclude):
        args=(['AAA'],['GGG'],['AAA','AAT','CCC'],['GGG','GGG','TTT'],[0,1,5]) if paired else (['AAA'],['AAA','AAT','CCC'],[0,1,5])
        call=gb.grouped_paired_sum_count_batch if paired else gb.grouped_count_batch
        return next(call(*args,**options,exclude_exact=exclude))
    full,masked=run(False),run(True)
    np.testing.assert_array_equal(np.asarray(masked.sizes),np.asarray(full.sizes))
    np.testing.assert_array_equal(np.asarray(masked.exact_counts),np.asarray(full.exact_counts))
    assert np.asarray(masked.counts)[0,1,0]==0 and np.asarray(full.counts)[0,1,0]==1
    assert np.asarray(masked.sizes)[0,1]==2
