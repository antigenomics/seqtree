"""Per-input positional geometry: exhaustive symmetric oracle and shared reductions."""
import numpy as np
import pytest
from seqtree import gapblock


def oracle(q, r, qw, rw, go=7, ge=3):
    d = abs(len(q)-len(r))
    scores = []
    for start in range(min(len(q), len(r))+1):
        score = go+(d-1)*ge if d else 0
        for i in range(min(len(q), len(r))):
            qi = i+(d if len(q)>len(r) and i>=start else 0)
            ri = i+(d if len(r)>len(q) and i>=start else 0)
            score += (q[qi]!=r[ri])*max(qw[qi], rw[ri])
        scores.append(score)
    return min(scores)


def test_sequence_weights_oracle_and_linked_reductions():
    q = ['', 'AC', 'ACA', 'CAAA']
    r = ['', 'CA', 'ACA', 'ACAAA', 'CA']
    qw = [[], [10,100], [100,10,10], [10,100,100,10]]
    rw = [[], [100,10], [10,100,10], [10,10,100,100,10], [10,100]]
    expected = np.array([[oracle(a,b,aw,bw) for b,bw in zip(r,rw)] for a,aw in zip(q,qw)])
    cuts = [[-1,0,10,100,250,100] for _ in q]
    opts = dict(gap_open=7, gap_extend=3, query_position_weights=qw, reference_position_weights=rw)
    for threads in (1,4):
        np.testing.assert_array_equal(np.asarray(gapblock.score_matrix(q,r,threads=threads,**opts)), expected)
        reverse = gapblock.score_matrix(r,q,gap_open=7,gap_extend=3,threads=threads,
            query_position_weights=rw,reference_position_weights=qw)
        np.testing.assert_array_equal(np.asarray(reverse).T,expected)
        counts,masses = gapblock.count_batch(q,r,cuts,threads=threads,linear_mass=True,exclude_exact=True,**opts)
        paircounts,pairmasses = gapblock.paired_sum_count_batch(q,q,r,r,cuts,gap_open=7,gap_extend=3,
            threads=threads,exclude_exact=True,query_position_weights_alpha=qw,
            query_position_weights_beta=qw,reference_position_weights_alpha=rw,reference_position_weights_beta=rw)
        for i,a in enumerate(q):
            eligible = np.array([a!=b for b in r])
            assert counts[i] == [int(np.sum(eligible & (expected[i]<=c))) for c in cuts[i]]
            assert masses[i] == [int(np.sum(np.maximum(0,c-expected[i][eligible]))) for c in cuts[i]]
            assert paircounts[i] == [int(np.sum(eligible & (2*expected[i]<=c))) for c in cuts[i]]
            assert pairmasses[i] == [int(np.sum(np.maximum(0,c-2*expected[i][eligible]))) for c in cuts[i]]
    uniform = dict(query_position_weights=[[100]*len(a) for a in q],reference_position_weights=[[100]*len(b) for b in r])
    old = gapblock.score_matrix(q,r,gap_open=7,gap_extend=3,position_weights_by_length={n:[100]*n for n in range(6)})
    np.testing.assert_array_equal(np.asarray(old),np.asarray(gapblock.score_matrix(q,r,gap_open=7,gap_extend=3,**uniform)))


@pytest.mark.parametrize('qw,rw,by_length', [([[1]],[[1,1]],None), ([[1,1]],None,None),
    ([[1,-1]],[[1,1]],None), ([[True,1]],[[1,1]],None), ([[1,1.5]],[[1,1]],None),
    ([[2**31-1]*2],[[1,1]],None), ([[1,1]],[[1,1]],{2:[1,1]})])
def test_sequence_weight_validation(qw,rw,by_length):
    with pytest.raises((ValueError,TypeError,OverflowError)):
        gapblock.count_batch(['AA'],['AC'],[[100]],query_position_weights=qw,
            reference_position_weights=rw,position_weights_by_length=by_length)
