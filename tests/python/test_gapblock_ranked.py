"""Exhaustive, small generic-sequence oracles for bounded gapblock reductions."""
import random
import pytest
import seqtree as st
from seqtree.gapblock import topk_batch, paired_topk_batch, count_batch, paired_count_batch, score_matrix, central_prior, positions_prior


def tuples(rows): return [[(h.ref_id,h.score) for h in row] for row in rows]


@pytest.mark.parametrize('matrix',[None,st.SubstitutionMatrix.blosum62()])
@pytest.mark.parametrize('prior',[None,central_prior(7)])
@pytest.mark.parametrize('exclude',[False,True])
def test_topk_and_counts_equal_exhaustive_matrix(matrix,prior,exclude):
    rng=random.Random(42)
    q=['','A','AC','ACD','ac']+[''.join(rng.choice('ACDE') for _ in range(rng.randrange(1,9))) for _ in range(12)]
    r=['','A','AC','AC','ACDE','ac']+[''.join(rng.choice('ACDE') for _ in range(rng.randrange(1,11))) for _ in range(31)]
    opts=dict(matrix=matrix,gap_open=9,gap_extend=2,gap_prior=prior)
    sm=score_matrix(q,r,**opts,threads=1)
    for k in (1,5,10,100):
        expected=[sorted([(j,sm[i,j]) for j,s in enumerate(r) if not exclude or s.upper()!=x.upper()],
                         key=lambda h:(h[1],h[0]))[:k] for i,x in enumerate(q)]
        for threads in (1,4,0):
            assert tuples(topk_batch(q,r,k,**opts,threads=threads,exclude_exact=exclude))==expected
    thresholds=[[30,0,-1,9,9,200] if i%2 else [] for i in range(len(q))]
    expected=[[sum(sm[i,j]<=t for j,s in enumerate(r) if not exclude or s.upper()!=x.upper())
               for t in thresholds[i]] for i,x in enumerate(q)]
    assert count_batch(q,r,thresholds,**opts,threads=1,exclude_exact=exclude)==expected
    assert count_batch(q,r,thresholds,**opts,threads=4,exclude_exact=exclude)==expected


def test_paired_topk_max_lane_oracle_and_joint_exact_puncture():
    qa,qb=['AC','DE',''],['ACD','A','D']
    ra,rb=['AC','AC','AD','DE',''],['ACD','AC','ACD','A','D']
    opts=dict(gap_open=3,gap_extend=2,gap_prior=central_prior(4))
    a,b=score_matrix(qa,ra,**opts),score_matrix(qb,rb,**opts)
    for exclude in (False,True):
        expected=[sorted([(j,a[i,j],b[i,j]) for j in range(len(ra))
                          if not exclude or (qa[i].upper(),qb[i].upper())!=(ra[j].upper(),rb[j].upper())],
                         key=lambda h:(max(h[1],h[2]),h[0]))[:3] for i in range(len(qa))]
        for threads in (1,4):
            assert paired_topk_batch(qa,qb,ra,rb,3,**opts,threads=threads,exclude_exact=exclude)==expected


def test_exclusion_before_selection_and_zero_penalty_nonidentity():
    matrix=st.SubstitutionMatrix.from_similarity([[0]*24 for _ in range(24)])
    assert tuples(topk_batch(['A'],['A','C','A','D'],1,matrix=matrix,exclude_exact=True))==[[(1,0)]]
    assert count_batch(['A'],['A','C','A','D'],[[0]],matrix=matrix,exclude_exact=True)==[[2]]
    assert paired_topk_batch(['A'],['C'],['A','A','D'],['C','D','C'],1,
                             matrix=matrix,exclude_exact=True)==[[(1,0,0)]]


def test_empty_inputs_and_alphabets():
    assert topk_batch([],['AC'],5)==[]
    assert topk_batch(['AC'],[],5)==[[]]
    assert count_batch(['AC'],[],[[0,2]])==[[0,0]]
    assert count_batch([],['AC'],[])==[]
    assert paired_topk_batch(['AC'],['A'],[],[],5)==[[]]
    assert tuples(topk_batch(['acgt'],['ACGT','ACGA'],2,alphabet='nt',exclude_exact=True))==[[(1,1)]]
    assert tuples(topk_batch(['N'],['N','R'],2,alphabet='iupac',exclude_exact=True))==[[(1,1)]]


@pytest.mark.parametrize('k',[0,-1,1.5,True,'5',2**40])
def test_invalid_k(k):
    with pytest.raises((ValueError,TypeError)): topk_batch(['AC'],['AC'],k)


@pytest.mark.parametrize('args',[
    dict(gap_open=-1),dict(gap_extend=-1),dict(threads=-1),dict(threads=1.5),
    dict(alphabet='bad'),dict(alphabet='nt',matrix=st.SubstitutionMatrix.blosum62()),
    dict(gap_prior=lambda i,d,m:-1),dict(gap_open=2**31-1),
])
def test_invalid_score_parameters(args):
    with pytest.raises((ValueError,TypeError,OverflowError)): topk_batch(['AC'],['A'],2,**args)


@pytest.mark.parametrize('thresholds',[[],[[1],[2]],[[1.5]],[[True]],[[2**40]]])
def test_invalid_threshold_shape_or_values(thresholds):
    with pytest.raises((ValueError,TypeError)): count_batch(['AC'],['A'],thresholds)


def test_invalid_symbols_even_when_other_axis_empty():
    for fn in (lambda:topk_batch(['#'],[],1),lambda:topk_batch([],['#'],1),
               lambda:count_batch(['#'],[],[[]]),lambda:count_batch([],['#'],[])):
        with pytest.raises(ValueError): fn()
    with pytest.raises(ValueError): paired_topk_batch(['A'],[],['A'],['A'])
    with pytest.raises(ValueError): paired_topk_batch(['A'],['A'],['A'],[])


@pytest.mark.parametrize('prior',[central_prior(7),positions_prior([1,-1])])
def test_paired_count_matches_all_matrix_ties_and_threshold_order(prior):
    qa,qb=['AC','DE',''],['ACD','A','D']
    ra,rb=['AC','AC','AD','DE',''],['ACD','AC','ACD','A','D']
    opts=dict(gap_open=3,gap_extend=2,gap_prior=prior)
    a,b=score_matrix(qa,ra,**opts),score_matrix(qb,rb,**opts)
    thresholds=[[20,-1,0,5,5],[],[30,2,0]]
    for exclude in (False,True):
        expected=[[sum(max(a[i,j],b[i,j])<=t for j in range(len(ra))
                       if not exclude or (qa[i],qb[i])!=(ra[j],rb[j])) for t in ts]
                  for i,ts in enumerate(thresholds)]
        for threads in (1,4):
            assert paired_count_batch(qa,qb,ra,rb,thresholds,**opts,threads=threads,exclude_exact=exclude)==expected
        selected=paired_topk_batch(qa,qb,ra,rb,3,**opts,exclude_exact=exclude)
        expected_top=[[j for j in sorted(range(len(ra)),key=lambda j:(max(a[i,j],b[i,j]),j))
                        if not exclude or (qa[i],qb[i])!=(ra[j],rb[j])][:3] for i in range(len(qa))]
        assert [[h[0] for h in row] for row in selected]==expected_top


def test_native_reference_accessor_preserves_ref_id_order_and_duplicates():
    refs=['AC','AC','AG','aC']
    ix=st.Index.build(refs,'aa')
    assert ix.ref_seqs()==[ix.ref_seq(i) for i in range(len(ix))]
    assert ix.ref_seqs()==refs
