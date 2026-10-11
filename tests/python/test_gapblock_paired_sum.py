"""Linked sum reduction agrees with dense lane sums without int32 truncation."""
import numpy as np
import pytest
from seqtree import gapblock


@pytest.mark.parametrize('exclude_exact',[False,True])
def test_paired_sum_dense_oracle_and_threads(exclude_exact):
    qa=['AC','ACC','','AC'];qb=['CA','AC','AC','CA']
    ra=['AC','AC','ACC','','CA'];rb=['CA','AC','AC','AC','CA']
    cut=[[-1,0,3,20,5,20],[],[0,3,20],[20,0]]
    common=dict(gap_open=3,gap_extend=2,gap_prior=gapblock.positions_prior((1,)))
    alpha=dict(position_weights_by_length={0:[],2:[0,2],3:[3,1,2]},
        query_group_ids=[0,1,1,0],reference_group_ids=[1,0,1,0,0],group_distances=[[7,2],[0,11]])
    beta=dict(position_weights_by_length={2:[1,3]},
        query_group_ids=[1,0,1,1],reference_group_ids=[1,0,0,1,0],group_distances=[[0,6],[4,0]])
    A=np.asarray(gapblock.score_matrix(qa,ra,**common,**alpha)).astype(np.int64)
    B=np.asarray(gapblock.score_matrix(qb,rb,**common,**beta)).astype(np.int64)
    expected_counts=[];expected_mass=[]
    for i,row in enumerate(cut):
        scores=[int(A[i,j]+B[i,j]) for j in range(len(ra))
                if not(exclude_exact and qa[i]==ra[j] and qb[i]==rb[j])]
        expected_counts.append([sum(d<=t for d in scores) for t in row])
        expected_mass.append([sum(max(0,t-d) for d in scores) for t in row])
    opts={**common,**{key+'_alpha':value for key,value in alpha.items()},
                   **{key+'_beta':value for key,value in beta.items()}}
    for threads in (1,4):
        assert gapblock.paired_sum_count_batch(qa,qb,ra,rb,cut,threads=threads,
            exclude_exact=exclude_exact,**opts)==(expected_counts,expected_mass)


def test_paired_sum_full_identity_puncture_and_ties():
    # Alpha alone matches every reference; only the first linked full pair is exact.
    assert gapblock.paired_sum_count_batch(['AC'],['CA'],['AC','AC'],['CA','AC'],[[0,2,3]],
        exclude_exact=True)==([[0,1,1]],[[0,0,1]])
    # Existing max contract remains different, unchanged.
    assert gapblock.paired_count_batch(['AC'],['AC'],['CA'],['CA'],[[2]])==[[1]]
    assert gapblock.paired_sum_count_batch(['AC'],['AC'],['CA'],['CA'],[[2]])==([[0]],[[0]])


def test_paired_sum_int32_overflow_and_empty_axes():
    maximum=2**31-1
    groups=dict(query_group_ids_alpha=[0],reference_group_ids_alpha=[0],group_distances_alpha=[[maximum]],
                query_group_ids_beta=[0],reference_group_ids_beta=[0],group_distances_beta=[[maximum]])
    assert gapblock.paired_sum_count_batch([''],[''],[''],[''],[[-1,0,maximum]],
        gap_open=0,gap_extend=0,**groups)==([[0,0,0]],[[0,0,0]])
    groups['group_distances_beta']=[[0]]
    assert gapblock.paired_sum_count_batch([''],[''],[''],[''],[[maximum]],
        gap_open=0,gap_extend=0,**groups)==([[1]],[[0]])
    assert gapblock.paired_sum_count_batch([],[],['AC'],['CA'],[])==([],[])
    assert gapblock.paired_sum_count_batch(['AC'],['CA'],[],[],[[0,3]])==([[0,0]],[[0,0]])


@pytest.mark.parametrize('args,opts',[
    ((['AC'],[],['AC'],['CA'],[[2]]),{}),
    ((['AC'],['CA'],['AC'],[],[[2]]),{}),
    ((['AC'],['CA'],['AC'],['CA'],[]),{}),
    ((['AC'],['CA'],['AC'],['CA'],[[2]]),{'position_weights_by_length_beta':{2:[-1,1]}}),
    ((['AC'],['CA'],['AC'],['CA'],[[2]]),{'query_group_ids_beta':[0]}),
])
def test_paired_sum_validation(args,opts):
    with pytest.raises((ValueError,TypeError)):
        gapblock.paired_sum_count_batch(*args,**opts)


@pytest.mark.parametrize('args,opts',[
    ((['AC'],[],['AC'],['CA'],[[2]]),{}),
    ((['AC'],['CA'],['AC'],['CA'],[[2]]),{'gap_open':-1}),
    ((['AC'],['CA'],['AC'],['CA'],[[2]]),{'query_group_ids_beta':[0]}),
    ((['AC'],['CA'],['AC'],['CA'],[[2]]),{'position_weights_by_length_beta':[[],[],[-1,1]]}),
])
def test_paired_sum_direct_native_validation(args,opts):
    from seqtree import _core
    with pytest.raises((ValueError,TypeError)):
        _core.gapblock_paired_sum_count_batch(*args,**opts)


def test_paired_sum_mass_uint64_and_puncture_with_zero_flanks():
    maximum=2**31-1
    assert gapblock.paired_sum_count_batch(['AC'],['CA'],['AC']*4,['CA']*4,[[maximum]])==(
        [[4]],[[4*maximum]])
    # Identical weighted cores are not identical full strings in both lanes.
    assert gapblock.paired_sum_count_batch(['AC'],['CA'],['CC'],['CA'],[[1]],
        position_weights_by_length_alpha={2:[0,1]},exclude_exact=True)==([[1]],[[1]])


def test_paired_sum_cell_addition_is_int64():
    # No group bound or first-lane bound rejects: only the actual lane sum exceeds int32.
    weight=2**30
    options=dict(gap_open=0,gap_extend=0,position_weights_by_length_alpha={1:[weight]},
                 position_weights_by_length_beta={1:[weight]})
    assert gapblock.paired_sum_count_batch(['A'],['A'],['C'],['C'],[[2**31-1]],**options)==([[0]],[[0]])
