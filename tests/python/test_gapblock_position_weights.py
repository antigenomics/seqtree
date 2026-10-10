"""Synthetic correctness checks for the queued positional single-block API only."""
import numpy as np
import pytest
from seqtree import Index,SearchParams,SubstitutionMatrix,PositionalMatrix,gapblock
B=SubstitutionMatrix.blosum62()

def exhaustive(q,r,weights,go=37,ge=11,prior=None):
    length=min(len(q),len(r));longer=max(len(q),len(r));d=longer-length
    scores=[]
    for start in range(length+1):
        value=(go+(d-1)*ge if d else 0)+(prior(start,d,longer) if prior and d else 0)
        for j in range(length):
            qi=j+(d if len(q)>len(r) and j>=start else 0)
            ri=j+(d if len(r)>len(q) and j>=start else 0)
            value+=weights[longer][max(qi,ri)]*B.penalty(q[qi],r[ri])
        scores.append(value)
    return min(scores)

@pytest.mark.parametrize('prior',[None,lambda i,d,m:3*abs(i-2)])
def test_positions_symmetry_reducers(prior):
    queries=['','A','AC','ACA','CACA','ACAAA']
    refs=queries+['ACA']
    weights={n:[1+j*17 for j in range(n)] for n in range(6)}
    options=dict(matrix=B,gap_open=37,gap_extend=11,gap_prior=prior,position_weights_by_length=weights)
    expected=np.array([[exhaustive(q,r,weights,prior=prior) for r in refs] for q in queries])
    for threads in (1,4):
        observed=np.asarray(gapblock.score_matrix(queries,refs,threads=threads,**options))
        np.testing.assert_array_equal(observed,expected)
        np.testing.assert_array_equal(np.asarray(gapblock.score_matrix(refs,queries,threads=threads,**options)).T,expected)
        thresholds=[[-1,0,37,100,500,100,10000] for q in queries]
        counts=gapblock.count_batch(queries,refs,thresholds,threads=threads,exclude_exact=True,**options)
        balls=gapblock.ball_batch(queries,refs,[100]*len(queries),threads=threads,exclude_exact=True,**options)
        for i,q in enumerate(queries):
            assert counts[i]==[sum(q!=r and expected[i,j]<=cut for j,r in enumerate(refs)) for cut in thresholds[i]]
            assert [(h.ref_id,h.score) for h in balls[i]]==[(j,int(expected[i,j])) for j,r in enumerate(refs) if q!=r and expected[i,j]<=100]
            assert all(h.n_subs==h.n_ins==h.n_dels==0 for h in balls[i])


def test_equal_length_positional_matrix():
    queries=['CASSLGQAYEQYF','CASALGQAYEQYF']
    refs=['CASSLGQAYEQYF','AASSLGQAYEQYF','CASSAGQAYEQYF','CASSLGQAYEQAF']
    length=len(queries[0]);weights={length:[117,62,56,85,117,156,153,152,98,62,71,74,98]}
    p=SearchParams(max_subs=5,max_ins=0,max_dels=0,max_total_edits=5,engine='seqtm');p.pos_matrix=PositionalMatrix.from_weights(B,weights[length])
    hits=Index.build(refs,'aa').search_batch(queries,p,threads=1)
    observed=np.asarray(gapblock.score_matrix(queries,refs,matrix=B,position_weights_by_length=weights,threads=1))
    for i,row in enumerate(hits):
        for hit in row:
            assert observed[i,hit.ref_id]==hit.score

@pytest.mark.parametrize('weights',[{2:[1]},{2:[1,-1]},{2:[1,1.5]},{2:[True,1]},{2:[1,2**31]},{3:[1,1,1]}, {2:[2**31-1]*2}])
@pytest.mark.parametrize('api',['score_matrix','count_batch','ball_batch'])
def test_validation(weights,api):
    args=(['AA'],['AC'])
    if api=='count_batch':args+=([[100]],)
    if api=='ball_batch':args+=([100],)
    with pytest.raises((ValueError,TypeError,OverflowError)):
        getattr(gapblock,api)(*args,matrix=B,position_weights_by_length=weights,threads=1)


def test_defaults_zero_weight_identity_duplicates_empty():
    queries=['AA'];refs=['AA','AC','AC']
    plain=np.asarray(gapblock.score_matrix(queries,refs,matrix=B,threads=1))
    np.testing.assert_array_equal(plain,np.asarray(gapblock.score_matrix(queries,refs,matrix=B,position_weights_by_length={},threads=4)))
    assert gapblock.count_batch(queries,refs,[[0]],matrix=B,position_weights_by_length={2:[0,0]},exclude_exact=True,threads=1)==[[2]]
    assert gapblock.count_batch([],refs,[],matrix=B,threads=1)==[]
    assert gapblock.count_batch(queries,[],[[0]],matrix=B,threads=1)==[[0]]
    assert gapblock.ball_batch(queries,refs,[-1],matrix=B,threads=1)==[[]]
    with pytest.raises(ValueError):gapblock.ball_batch(queries,refs,[],threads=1)


def test_matrix_negative_threads_preserves_automatic_concurrency():
    q=['AC','CA'];r=['AA','AC','CA']
    for weights in (None,{2:[3,7]}):
        args=dict(matrix=B,position_weights_by_length=weights)
        expected=np.asarray(gapblock.score_matrix(q,r,threads=1,**args))
        np.testing.assert_array_equal(expected,np.asarray(gapblock.score_matrix(q,r,threads=-1,**args)))
