"""Integer linear taper shares exact counts, scoring and full-identity puncture."""
import numpy as np
import pytest
from seqtree import gapblock, SubstitutionMatrix


def scalar(scores, thresholds, queries, refs, exclude_exact):
    counts=[];masses=[]
    for i,row in enumerate(thresholds):
        eligible=[int(scores[i,j]) for j in range(len(refs))
                  if not(exclude_exact and queries[i]==refs[j])]
        counts.append([sum(d<=t for d in eligible) for t in row])
        masses.append([sum(max(0,t-d) for d in eligible) for t in row])
    return counts,masses


@pytest.mark.parametrize('exclude_exact',[False,True])
def test_linear_mass_scalar_oracle_and_threads(exclude_exact):
    q=['AC','ACC','','AC'];r=['AC','CA','ACC','CA','']
    thresholds=[[-1,0,2,5,20,5],[],[0,3,20],[20,0]]
    options=dict(gap_open=3,gap_extend=2,position_weights_by_length={0:[],2:[0,2],3:[3,1,2]},
                 query_group_ids=[0,1,1,0],reference_group_ids=[1,0,1,0,0],group_distances=[[7,2],[0,11]])
    scores=np.asarray(gapblock.score_matrix(q,r,threads=1,**options))
    expected=scalar(scores,thresholds,q,r,exclude_exact)
    for threads in (1,4):
        actual=gapblock.count_batch(q,r,thresholds,threads=threads,exclude_exact=exclude_exact,
                                   linear_mass=True,**options)
        assert actual==expected
        assert actual[0]==gapblock.count_batch(q,r,thresholds,threads=threads,
                                              exclude_exact=exclude_exact,**options)


def test_linear_mass_boundaries_and_uint64():
    assert gapblock.count_batch(['AC'],['AC','CA'],[[0,2,1,-1,2]],linear_mass=True)==(
        [[1,2,1,0,2]],[[0,2,1,0,2]])
    maximum=2**31-1
    assert gapblock.count_batch(['AC'],['AC']*4,[[maximum]],linear_mass=True)==([[4]],[[4*maximum]])
    assert gapblock.count_batch(['AC'],['AC'],[[maximum]],exclude_exact=True,linear_mass=True)==([[0]],[[0]])
    assert gapblock.count_batch(['AC'],['AC'],[[5]],linear_mass=True,
        query_group_ids=[0],reference_group_ids=[0],group_distances=[[6]])==([[0]],[[0]])


def test_linear_mass_empty_axes_and_prior():
    assert gapblock.count_batch([],['AC'],[],linear_mass=True)==([],[])
    assert gapblock.count_batch(['AC'],[],[[-1,0,2]],linear_mass=True)==([[0,0,0]],[[0,0,0]])
    q=['AC','ACC'];r=['ACC','AC']
    opts=dict(gap_prior=gapblock.positions_prior((1,)),matrix=SubstitutionMatrix.blosum62(),gap_open=28,gap_extend=14)
    scores=np.asarray(gapblock.score_matrix(q,r,**opts))
    cut=[[0,28,50],[50,28,0]]
    assert gapblock.count_batch(q,r,cut,linear_mass=True,**opts)==scalar(scores,cut,q,r,False)


@pytest.mark.parametrize('value',[1,None,'yes'])
def test_linear_mass_requires_explicit_boolean(value):
    with pytest.raises((ValueError,TypeError)):
        gapblock.count_batch(['AC'],['CA'],[[2]],linear_mass=value)


@pytest.mark.parametrize('options',[
    {'gap_open':-1}, {'gap_extend':-1},
    {'query_group_ids':[0],'reference_group_ids':[0],'group_distances':[[-1]]},
    {'position_weights_by_length':[[],[],[-1,1]]},
])
def test_linear_mass_native_validation(options):
    from seqtree import _core
    with pytest.raises((ValueError,TypeError,OverflowError)):
        _core.gapblock_count_mass_batch(['AC'],['CA'],[[2]],**options)
