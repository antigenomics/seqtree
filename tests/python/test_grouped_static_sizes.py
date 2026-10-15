"""Static population sizes equal the unchanged full-membership traversal."""
import numpy as np
import pytest
from seqtree import gapblock as gb


@pytest.mark.parametrize('paired',[False,True])
@pytest.mark.parametrize('exclude',[False,True])
@pytest.mark.parametrize('threads',[1,4])
@pytest.mark.parametrize('cuts',[[],[-1],[-1,0,3,10,40]])
def test_all_buffers_match_no_withholding_full_membership_oracle(paired,exclude,threads,cuts,require_thread_budget):
    q=['AAA','CCC','ATA'];r=['AAA','AAA','ACA','TTT','CCC']
    qb=['CCC','AAA','CCC'];rb=['CCC','CCA','CCC','AAA','AAA']
    metadata=dict(reference_target_offsets=[0,1,3,5,5,5],
        reference_target_ids=[0,1,2,0,2],target_peptide_ids=[0,1,1],
        reference_exclusive_peptide_ids=[0,1,-1,0,-1],
        query_gene_ids=[0,42,2147483647],reference_gene_ids=[1,0,1,2147483647,1],chunk_size=2)
    full=dict(query_withheld_peptide_ids=[-1]*3,
        reference_peptide_offsets=[0,1,2,4,5,7],reference_peptide_ids=[0,1,0,1,0,0,1])
    features=gb.PairedFeatureOffsets([[0,0],[1,1],[0,1]],
        [[0,0],[1,1],[1,0],[0,1],[1,1]],[[0,2],[2,0]],weight=.5,scale=2)
    extra=gb.PairedFeatureOffsets([[0,0]]*3,[[0,0]]*5,[[2]],scale=2)
    lane=dict(query_group_ids=[0,1,2],reference_group_ids=[0,1,1,0,1],
        group_distances=[[0,3],[2,0],[1,4]],feature_offsets=features,
        additional_feature_offsets=[extra],query_position_weights=[[1,2,3],[3,1,2],[2,3,1]],
        reference_position_weights=[[1,1,1],[2,1,1],[1,3,1],[1,1,3],[3,1,1]])
    options=dict(gap_open=2,gap_extend=1,threads=threads,exclude_exact=exclude)
    if paired:
        options.update({key+'_'+suffix:value for suffix in ('alpha','beta') for key,value in lane.items()})
        options.update(score_quantum_alpha=2,score_quantum_beta=3)
        call=lambda folds:gb.grouped_paired_sum_count_batch(q,qb,r,rb,cuts,**metadata,**options,**folds)
    else:
        options.update(lane,score_quantum=2)
        call=lambda folds:gb.grouped_count_batch(q,r,cuts,**metadata,**options,**folds)
    def buffers(batch):
        chunks=list(batch)
        return [np.concatenate([np.asarray(getattr(c,name)) for c in chunks])
            for name in ('counts','sizes','exact_counts')]
    actual=buffers(call({}));original=buffers(call(full))
    for new,old in zip(actual,original):
        assert new.dtype==old.dtype==np.uint32
        assert new.shape==old.shape and new.tobytes()==old.tobytes()
    members=[[0],[1,2],[0,2],[],[]]
    for i in range(3):
        exact=[j for j in range(5) if q[i]==r[j] and (not paired or qb[i]==rb[j])]
        retained=[j for j in range(5) if not (exclude and j in exact)]
        target=[[j for j in retained if t in members[j]] for t in range(3)]
        same_gene=[[j for j in target[t] if metadata['query_gene_ids'][i]==metadata['reference_gene_ids'][j]] for t in range(3)]
        exclusive=[[j for j in retained if metadata['reference_exclusive_peptide_ids'][j]==p] for p in range(2)]
        np.testing.assert_array_equal(actual[1][i],[len(retained)]+list(map(len,target+same_gene+exclusive)))
        np.testing.assert_array_equal(actual[2][i],[sum(t in members[j] for j in exact) for t in range(3)])
