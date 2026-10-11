"""Streaming group audit preserves exact component/full balls without hit materialization."""
import numpy as np
import pytest
from seqtree import gapblock

def test_group_audit_matches_dense_components_and_threads():
    q=['AA','CA','aa',''];r=['AA','CC','CAA','AA','']
    opts=dict(gap_open=2,gap_extend=1,position_weights_by_length={0:[],2:[0,3],3:[2,1,0]})
    group=dict(query_group_ids=[0,1,0,1],reference_group_ids=[0,1,0,1,1],group_distances=[[0,5],[3,0]])
    base=np.asarray(gapblock.score_matrix(q,r,**opts));full=np.asarray(gapblock.score_matrix(q,r,**opts,**group))
    limits=[0,3,10,-1]
    for threads in [1,2]:
        for excluded in [False,True]:
            before,after,nearest=gapblock.group_audit_batch(q,r,limits,threads=threads,exclude_exact=excluded,**opts,**group)
            for i,t in enumerate(limits):
                valid=np.array([not excluded or q[i].upper()!=s.upper() for s in r]);selected=valid&(base[i]<=t)
                assert before[i]==int(selected.sum())
                assert after[i]==int((valid&(full[i]<=t)).sum())
                assert nearest[i]==(int(full[i,selected].min()) if selected.any() else -1)
    assert gapblock.group_audit_batch(['AA'],[],[0])==([0],[0],[-1])
    with pytest.raises(ValueError):gapblock.group_audit_batch(q,r,[])
