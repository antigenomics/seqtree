"""Shared cutoffs reproduce expanded-row batch curves with one validation pass."""
import pytest
from seqtree import _core
from seqtree.gapblock import PairedFeatureOffsets, count_batch, paired_sum_count_batch


@pytest.mark.parametrize('threads',[1,4])
@pytest.mark.parametrize('exclude',[False,True])
def test_shared_thresholds_match_expanded(threads,exclude, require_thread_budget):
    q=['AAC','AC','AAC']; r=['AAC','ACC','AC','AAC']
    qb=['ACC','AC','AAC']; rb=['ACC','AAC','AC','AAC']
    cuts=[-1,0,3,1,3,12,20]
    f=PairedFeatureOffsets([[0,0],[1,1],[0,0]],[[0,0],[1,1],[1,1],[0,0]],[[0,3],[2,0]],weight=.5,scale=2)
    opts=dict(query_group_ids=[0,1,0],reference_group_ids=[0,1,1,0],group_distances=[[0,3],[2,0]],feature_offsets=f,gap_open=2,threads=threads,exclude_exact=exclude)
    for mass in [False,True]:
        assert count_batch(q,r,cuts,shared_thresholds=True,linear_mass=mass,**opts)==count_batch(q,r,[cuts]*len(q),linear_mass=mass,**opts)
    pair_opts=dict(threads=threads,exclude_exact=exclude,gap_open=2,feature_offsets_alpha=f,feature_offsets_beta=f,
        query_group_ids_alpha=opts['query_group_ids'],reference_group_ids_alpha=opts['reference_group_ids'],group_distances_alpha=opts['group_distances'])
    assert paired_sum_count_batch(q,qb,r,rb,cuts,shared_thresholds=True,**pair_opts)==paired_sum_count_batch(q,qb,r,rb,[cuts]*len(q),**pair_opts)


def test_shared_empty_axes_cutoffs_and_generators():
    assert count_batch([],['A'],[0,3],shared_thresholds=True,threads=1)==[]
    assert count_batch([],[],[],shared_thresholds=True,linear_mass=True,threads=1)==([],[])
    assert count_batch(['A','C'],[],[2,0,2],shared_thresholds=True,linear_mass=True,threads=1)==([[0,0,0],[0,0,0]],[[0,0,0],[0,0,0]])
    assert count_batch(['A','C'],['A'],[],shared_thresholds=True,linear_mass=True,threads=1)==([[],[]],[[],[]])
    assert count_batch(['A'],['C'],iter([1,2]),shared_thresholds=True,threads=1)==[[1,1]]
    assert paired_sum_count_batch([],[],['A'],['C'],[0,3],shared_thresholds=True,threads=1)==([],[])
    assert paired_sum_count_batch(['A'],['C'],[],[],[],shared_thresholds=True,threads=1)==([[]],[[]])


def test_shared_validates_values_once_per_call():
    class Integer:
        def __init__(self,value): self.value=value; self.calls=0
        def __index__(self): self.calls+=1; return self.value
    cuts=[Integer(i) for i in range(153)]
    q=['A']*100
    count_batch(q,[],cuts,shared_thresholds=True,threads=1)
    assert [v.calls for v in cuts]==[1]*153
    count_batch(q,[],cuts,shared_thresholds=True,linear_mass=True,threads=1)
    assert [v.calls for v in cuts]==[2]*153
    paired_sum_count_batch(q,q,[],[],cuts,shared_thresholds=True,threads=1)
    assert [v.calls for v in cuts]==[3]*153


@pytest.mark.parametrize('cuts',[[True],[1.5],[[1]],[-(1<<31)-1],[1<<31]])
def test_shared_bad_values(cuts):
    with pytest.raises(ValueError): count_batch(['A'],['C'],cuts,shared_thresholds=True,threads=1)
    with pytest.raises(ValueError): paired_sum_count_batch(['A'],['C'],['A'],['C'],cuts,shared_thresholds=True,threads=1)


@pytest.mark.parametrize('flag',[1,None,'yes'])
def test_shared_bad_flag(flag):
    with pytest.raises(TypeError,match='boolean'): count_batch(['A'],['C'],[[1]],shared_thresholds=flag,threads=1)
    with pytest.raises(TypeError,match='boolean'): paired_sum_count_batch(['A'],['C'],['A'],['C'],[[1]],shared_thresholds=flag,threads=1)


def test_shared_native_row_shape():
    for rows in [[],[[0],[1]]]:
        with pytest.raises(ValueError,match='exactly one row'):
            _core.gapblock_count_batch(['A'],['C'],rows,threads=1,shared_thresholds=True)
    assert _core.gapblock_count_batch([],[],[[0,1]],threads=1,shared_thresholds=True)==[]
    assert _core.gapblock_count_batch(['A','C'],['A'],[[]],threads=1,shared_thresholds=True)==[[],[]]
