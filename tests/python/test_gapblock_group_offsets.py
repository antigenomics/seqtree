"""Generic additive group penalties share the existing gapblock cell and reducers."""
import numpy as np
import pytest
from seqtree import gapblock


def test_group_offsets_matrix_reducers_and_puncture():
    q=['AC','ACC','','AC']; r=['AC','CA','ACC','CA','']
    qg=[0,1,1,0]; rg=[1,0,1,0,0]; distances=[[7,2],[0,11]]
    opts=dict(gap_open=3,gap_extend=2,position_weights_by_length={0:[],2:[0,2],3:[3,1,2]})
    base=np.asarray(gapblock.score_matrix(q,r,threads=1,**opts))
    expected=base+np.array([[distances[a][b] for b in rg] for a in qg])
    group=dict(query_group_ids=qg,reference_group_ids=rg,group_distances=distances)
    thresholds=[[-1,0,2,5,20,5],[],[0,3,20],[20]]
    for t in (1,4):
        np.testing.assert_array_equal(np.asarray(gapblock.score_matrix(q,r,threads=t,**opts,**group)),expected)
        for excluded in (False,True):
            counts=gapblock.count_batch(q,r,thresholds,threads=t,exclude_exact=excluded,**opts,**group)
            balls=gapblock.ball_batch(q,r,[5]*len(q),threads=t,exclude_exact=excluded,**opts,**group)
            for i in range(len(q)):
                eligible=lambda j:not(excluded and q[i]==r[j])
                assert counts[i]==[sum(eligible(j) and expected[i,j]<=c for j in range(len(r))) for c in thresholds[i]]
                assert [(h.ref_id,h.score) for h in balls[i]]==[(j,int(expected[i,j])) for j in range(len(r)) if eligible(j) and expected[i,j]<=5]
    # A full-junction exact match remains punctured even with a positive group penalty.
    assert gapblock.count_batch(['AC'],['AC'],[[100]],exclude_exact=True,query_group_ids=[0],reference_group_ids=[1],group_distances=distances)==[[0]]
    # Offset alone rejects the entire group even when its junction distance is zero.
    assert gapblock.ball_batch(['AC'],['AC'],[1],query_group_ids=[0],reference_group_ids=[1],group_distances=distances)==[[]]


@pytest.mark.parametrize('api',['score_matrix','count_batch','ball_batch'])
@pytest.mark.parametrize('groups',[
 dict(query_group_ids=[0]),
 dict(query_group_ids=[0],reference_group_ids=[0],group_distances=[]),
 dict(query_group_ids=[],reference_group_ids=[0],group_distances=[[0]]),
 dict(query_group_ids=[-1],reference_group_ids=[0],group_distances=[[0]]),
 dict(query_group_ids=[1],reference_group_ids=[0],group_distances=[[0]]),
 dict(query_group_ids=[0],reference_group_ids=[1],group_distances=[[0]]),
 dict(query_group_ids=[0.5],reference_group_ids=[0],group_distances=[[0]]),
 dict(query_group_ids=[False],reference_group_ids=[0],group_distances=[[0]]),
 dict(query_group_ids=[0],reference_group_ids=[0],group_distances=[[0],[]]),
 dict(query_group_ids=[0],reference_group_ids=[0],group_distances=[[-1]]),
 dict(query_group_ids=[0],reference_group_ids=[0],group_distances=[[1.5]]),
 dict(query_group_ids=[0],reference_group_ids=[0],group_distances=[[2**31-1]]),
])
def test_group_validation(api,groups):
    args=(['AC'],['CA'])
    if api=='count_batch':args+=([[10]],)
    if api=='ball_batch':args+=([10],)
    with pytest.raises((ValueError,TypeError,OverflowError)):
        getattr(gapblock,api)(*args,threads=1,**groups)


def test_defaults_and_empty_axes():
    q=['AC',''];r=['AC','CA']
    base=np.asarray(gapblock.score_matrix(q,r,threads=1))
    zero=dict(query_group_ids=[0,0],reference_group_ids=[0,0],group_distances=[[0]])
    np.testing.assert_array_equal(base,np.asarray(gapblock.score_matrix(q,r,threads=-1,**zero)))
    assert gapblock.count_batch([],r,[],query_group_ids=[],reference_group_ids=[0,0],group_distances=[[5]])==[]
    assert gapblock.count_batch(q,[],[[0],[5]],query_group_ids=[0,0],reference_group_ids=[],group_distances=[[5]])==[[0],[0]]


@pytest.mark.parametrize('api',['gapblock_matrix','gapblock_count_batch','gapblock_ball_batch'])
@pytest.mark.parametrize('groups',[
 dict(query_group_ids=[],reference_group_ids=[0],group_distances=[[0]]),
 dict(query_group_ids=[-1],reference_group_ids=[0],group_distances=[[0]]),
 dict(query_group_ids=[0],reference_group_ids=[0],group_distances=[[0],[]]),
 dict(query_group_ids=[0],reference_group_ids=[0],group_distances=[[-1]]),
 dict(query_group_ids=[0],reference_group_ids=[0],group_distances=[[2**31-1]]),
])
def test_native_validation_before_workers(api,groups):
    from seqtree import _core
    args=(['AC'],['CA'])
    if api=='gapblock_count_batch':args+=([[10]],)
    if api=='gapblock_ball_batch':args+=([10],)
    with pytest.raises(ValueError):
        getattr(_core,api)(*args,threads=4,**groups)


def test_rectangular_group_table():
    opts=dict(query_group_ids=[0],reference_group_ids=[2,0],group_distances=[[7,4,1]])
    observed=np.asarray(gapblock.score_matrix(['AC'],['AC','AC'],threads=1,**opts))
    np.testing.assert_array_equal(observed,[[1,7]])


def test_int32_boundary_group_only_score():
    maximum=2**31-1
    opts=dict(gap_open=0,gap_extend=0,query_group_ids=[0],reference_group_ids=[0],group_distances=[[maximum]],threads=4)
    assert gapblock.score_matrix([''],[''],**opts)[0,0]==maximum
    assert gapblock.count_batch([''],[''],[[maximum-1,maximum]],**opts)==[[0,1]]
    assert [(h.ref_id,h.score) for h in gapblock.ball_batch([''],[''],[maximum],**opts)[0]]==[(0,maximum)]
    assert gapblock.count_batch([''],[''],[[maximum]],exclude_exact=True,**opts)==[[0]]
