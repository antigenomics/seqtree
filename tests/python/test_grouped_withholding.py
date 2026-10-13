"""Withheld label counts equal independent full fold reconstruction."""
import numpy as np
import pytest
from seqtree import gapblock as gb

Q = ['AAA','AAA','CCC','AAA','GGG','TTT']
R = ['AAA','AAA','AAC','CCC','GGG','TTT']
QB = ['CCC','CCC','AAA','CCC','TTT','CCC']
RB = ['CCC','AAC','CCC','AAA','TTT','CCC']
MEMBERS = [[0,1,2],[0],[1,2],[3],[4],[]]
PARENTS = [0,0,1,2,3]
FULL = [[0,1],[0],[0,1],[2],[3,4],[4]]
WITHHELD = [0,1,-1,2,4,0]
CUTS = [-1,0,2,5,20]
QGENES = [0,1,1,0,0,1]
RGENES = [1,0,1,0,0,1]


def csr(rows):
    offsets=[0];ids=[]
    for row in rows:ids.extend(row);offsets.append(len(ids))
    return offsets,ids


def metadata(indices=range(6), withheld=None):
    indices=list(indices)
    full=[FULL[j] if withheld is None else [p for p in FULL[j] if p!=withheld] for j in indices]
    members=[[t for t in MEMBERS[j] if withheld is None or PARENTS[t]!=withheld] for j in indices]
    offsets,ids=csr(members)
    return dict(reference_target_offsets=offsets,reference_target_ids=ids,
        target_peptide_ids=PARENTS,
        reference_exclusive_peptide_ids=[row[0] if len(row)==1 and row[0]<4 else -1 for row in full],
        query_gene_ids=QGENES,reference_gene_ids=[RGENES[j] for j in indices],chunk_size=2)


def options(indices, paired, features):
    indices=list(indices)
    result=dict(gap_open=2,gap_extend=1)
    for lane in range(2 if paired else 1):
        suffix=('_alpha' if lane==0 else '_beta') if paired else ''
        result['query_group_ids'+suffix]=[0,1,0,1,1,0]
        result['reference_group_ids'+suffix]=[[1,0,1,0,1,0][j] for j in indices]
        result['group_distances'+suffix]=[[0,3],[3,2]]
        if features:
            result['feature_offsets'+suffix]=gb.PairedFeatureOffsets(
                [[0,1],[1,0],[1,1],[0,1],[1,0],[0,0]],
                [[[1,0],[0,0],[0,1],[1,1],[1,0],[0,1]][j] for j in indices],
                [[0,3],[3,0]],weight=.5 if lane==0 else 1.5)
    return result


def run(paired, indices, args, features, threads, exclude):
    opt=options(indices,paired,features)
    if paired:
        chunks=list(gb.grouped_paired_sum_count_batch(Q,QB,[R[j] for j in indices],
            [RB[j] for j in indices],CUTS,**args,**opt,threads=threads,exclude_exact=exclude))
    else:
        chunks=list(gb.grouped_count_batch(Q,[R[j] for j in indices],CUTS,
            **args,**opt,threads=threads,exclude_exact=exclude))
    return tuple(np.concatenate([np.asarray(getattr(c,key)) for c in chunks])
                 for key in ('counts','sizes','exact_counts'))


@pytest.mark.parametrize('paired',[False,True])
@pytest.mark.parametrize('exclude',[False,True])
@pytest.mark.parametrize('features',[False,True])
@pytest.mark.parametrize('threads',[1,4])
def test_every_fold_matches_independently_rebuilt_population(paired,exclude,features,threads):
    args=metadata();offsets,ids=csr(FULL)
    args.update(query_withheld_peptide_ids=WITHHELD,
                reference_peptide_offsets=offsets,reference_peptide_ids=ids)
    actual=run(paired,list(range(6)),args,features,threads,exclude)
    for i, withheld in enumerate(WITHHELD):
        indices=[j for j,row in enumerate(FULL) if any(p!=withheld for p in row)]
        expected=run(paired,indices,metadata(indices,withheld),features,threads,exclude)
        for output,oracle in zip(actual,expected):
            np.testing.assert_array_equal(output[i],oracle[i])
    # Same sequence in different folds changes surviving target/exclusive channels.
    assert actual[1][0,1]==0
    assert actual[1][1,2]>0


def test_no_withholding_vector_matches_existing_default():
    args=metadata();offsets,ids=csr(FULL)
    expected=run(False,list(range(6)),args,False,1,False)
    args.update(query_withheld_peptide_ids=[-1]*6,
        reference_peptide_offsets=offsets,reference_peptide_ids=ids)
    for output,oracle in zip(run(False,list(range(6)),args,False,1,False),expected):
        np.testing.assert_array_equal(output,oracle)


@pytest.mark.parametrize('change',[
    {'query_withheld_peptide_ids':[0]},
    {'query_withheld_peptide_ids':[-2]*6},
    {'query_withheld_peptide_ids':[True]*6},
    {'reference_peptide_offsets':None},
    {'reference_peptide_ids':None},
    {'reference_peptide_offsets':[0]},
    {'reference_peptide_ids':[0,0,0,0,1,2,3,4,4]},
    {'reference_peptide_ids':[0,2,0,0,1,2,3,4,4]},
    {'reference_exclusive_peptide_ids':[0,0,-1,2,-1,-1]},
    {'reference_peptide_ids':[-1,1,0,0,1,2,3,4,4]},
    {'reference_peptide_ids':[0,1,0,0,1,2,3,4,1<<31]},
])
def test_invalid_withholding_metadata_rejected_before_iteration(change):
    args=metadata();offsets,ids=csr(FULL)
    args.update(query_withheld_peptide_ids=WITHHELD,
        reference_peptide_offsets=offsets,reference_peptide_ids=ids)
    args.update(change)
    with pytest.raises((ValueError,TypeError,OverflowError)):
        gb.grouped_count_batch(Q,R,CUTS,**args)


def test_full_memberships_require_explicit_query_fold_axis():
    args=metadata();offsets,ids=csr(FULL)
    with pytest.raises(ValueError):
        gb.grouped_count_batch(Q,R,CUTS,**args,reference_peptide_offsets=offsets,reference_peptide_ids=ids)


def test_empty_fold_axes_and_all_memberships_removed():
    empty=dict(reference_target_offsets=[0],reference_target_ids=[],target_peptide_ids=[],
        reference_exclusive_peptide_ids=[],query_gene_ids=[],reference_gene_ids=[],
        query_withheld_peptide_ids=[],reference_peptide_offsets=[0],reference_peptide_ids=[])
    assert list(gb.grouped_count_batch([],[],[],**empty))==[]
    empty.update(query_gene_ids=[0],query_withheld_peptide_ids=[0])
    chunk=next(gb.grouped_count_batch(['AAA'],[],[],**empty))
    assert np.asarray(chunk.counts).shape==(1,1,0)
    assert np.asarray(chunk.sizes).tolist()==[[0]]
    args=dict(reference_target_offsets=[0,1],reference_target_ids=[0],target_peptide_ids=[0],
        reference_exclusive_peptide_ids=[0],query_gene_ids=[0],reference_gene_ids=[0],
        query_withheld_peptide_ids=[0],reference_peptide_offsets=[0,1],reference_peptide_ids=[0])
    chunk=next(gb.grouped_count_batch(['AAA'],['AAA'],CUTS,**args))
    assert not np.asarray(chunk.counts).any()
    assert not np.asarray(chunk.sizes).any()
    assert not np.asarray(chunk.exact_counts).any()
    args.update(reference_target_ids=[],reference_target_offsets=[0,0],
        reference_exclusive_peptide_ids=[-1],reference_peptide_offsets=[0,0],reference_peptide_ids=[])
    chunk=next(gb.grouped_count_batch(['AAA'],['AAA'],CUTS,**args))
    assert not np.asarray(chunk.sizes).any()


def test_fold_metadata_owned_and_unselected_only_survivors():
    args=metadata();offsets,ids=csr(FULL)
    folds=WITHHELD.copy()
    args.update(query_withheld_peptide_ids=folds,reference_peptide_offsets=offsets,reference_peptide_ids=ids)
    expected=run(False,list(range(6)),args,False,1,True)
    batch=gb.grouped_count_batch(Q,R,CUTS,**args,**options(range(6),False,False),exclude_exact=True)
    folds.clear();offsets.clear();ids.clear()
    chunks=list(batch)
    for key,oracle in zip(('counts','sizes','exact_counts'),expected):
        np.testing.assert_array_equal(np.concatenate([np.asarray(getattr(c,key)) for c in chunks]),oracle)
    assert not np.asarray(chunks[0].counts).flags.writeable
    args=dict(reference_target_offsets=[0,0],reference_target_ids=[],target_peptide_ids=[],
        reference_exclusive_peptide_ids=[-1],query_gene_ids=[0],reference_gene_ids=[0],
        query_withheld_peptide_ids=[0],reference_peptide_offsets=[0,1],reference_peptide_ids=[4])
    chunk=next(gb.grouped_count_batch(['AAA'],['CCC'],CUTS,**args))
    assert np.asarray(chunk.sizes).tolist()==[[1]]
    assert np.asarray(chunk.counts)[0,0,-1]==1
