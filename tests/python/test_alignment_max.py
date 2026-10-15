"""Methods-style floating NW is independent of candidate scope and final coefficients."""
import math
import numpy as np
import pytest
from seqtree import Index, SearchParams


def nw(q,r,symbols,matrix):
    """Independent full-matrix three-state oracle; strict ties preserve first state."""
    lookup={a:i for i,a in enumerate(symbols)};gap=max(abs(x) for row in matrix for x in row)
    bad=(-math.inf,0.,0); m,n=len(q),len(r)
    states=[[[bad for _ in range(n+1)] for _ in range(m+1)] for _ in range(3)]
    states[0][0][0]=(0.,0.,0)
    for j in range(1,n+1):states[2][0][j]=(-gap*j,0.,j)
    for i in range(1,m+1):states[1][i][0]=(-gap*i,0.,i)
    pick=lambda values:max(values,key=lambda x:x[0])
    for i in range(1,m+1):
        for j in range(1,n+1):
            s=matrix[lookup[q[i-1]]][lookup[r[j-1]]]
            a=pick([states[k][i-1][j-1] for k in (0,1,2)])
            states[0][i][j]=(a[0]+s,a[1]+s,a[2])
            for k,ii,jj,order in [(1,i-1,j,(0,1,2)),(2,i,j-1,(0,2,1))]:
                a=pick([states[t][ii][jj] for t in order]);states[k][i][j]=(a[0]-gap,a[1],a[2]+1)
    _,raw,gaps=pick([states[k][m][n] for k in (0,1,2)])
    self_score=lambda seq:sum(matrix[lookup[a]][lookup[a]] for a in seq)
    return raw,raw-max(self_score(q),self_score(r)),gaps


def run(queries,refs,*,matrix=None,symbols='ACGT',target_ids=None,requests=None,exclude=False,groups=None,other=None,roles=None,coef=None,params=None,cap=2,threads=1,updates=None):
    matrix=matrix or [[2.25 if a==b else -.75 for b in symbols] for a in symbols]
    target_ids=target_ids or list(range(len(refs)));targets=max(target_ids)+1
    requests=requests or [list(range(targets)) for _ in queries]
    qoff=np.cumsum([0]+[len(r) for r in requests]).tolist();qt=[t for row in requests for t in row]
    roles=roles or [[t] for t in target_ids]
    roff=np.cumsum([0]+[len(r) for r in roles]).tolist();rc=[t for row in roles for t in row]
    kwargs=dict(query_target_offsets=qoff,query_target_ids=qt,reference_target_offsets=list(range(len(refs)+1)),reference_target_ids=target_ids,
        target_class_ids=list(range(targets)),reference_class_offsets=roff,reference_class_ids=rc,
        query_group_ids=groups or [[0]*len(queries),[0]*len(queries)],reference_group_ids=[[0]*len(refs),[0]*len(refs)],group_offsets=[[[.1]],[[.2]]],
        reference_other_counts=other or [0]*len(refs),reference_multiplicities=[1]*len(refs),coefficients=coef or [.5,.15,-.7],threads=threads,exclude_exact=exclude)
    if updates:kwargs.update(updates)
    result=Index.build(refs,alphabet='nt').alignment_max_batch(queries,params or SearchParams(max_subs=4,max_ins=2,max_dels=2,max_total_edits=6),matrix,symbols,cap,**kwargs)
    return {key:np.asarray(value) for key,value in result.arrays().items()}


def test_fractional_raw_self_normalization_and_independent_nw():
    symbols='ACGT';matrix=[[1.25,-.75,-1.,-.5],[-.6,3.5,-.8,-.2],[-.7,-.9,2.,-.3],[-.4,-.5,-.6,1.]]
    q='AAAAACCCCC';r='CCCCCAAAAA'
    for coef in ([.5,.15,-100.],[.5,2.,0.]):
        out=run([q],[r,'TTTTTTTTTT'],matrix=matrix,requests=[[0]],coef=coef,params=SearchParams(max_subs=20,max_ins=10,max_dels=10,max_total_edits=30),cap=20)
        raw,norm,gap=nw(q,r,symbols,matrix)
        assert out['raw_similarity'][0]==pytest.approx(raw,abs=1e-12)
        assert out['normalized_similarity'][0]==pytest.approx(norm,abs=1e-12)
        assert out['gap_counts'][0]==gap
        eta=coef[0]+.3+coef[1]*norm+coef[2]*gap
        pfit=-math.expm1(-math.exp(eta))
        assert out['pfit'][0]==pytest.approx(pfit,rel=1e-13,abs=1e-15)
        assert out['scores'][0]==pytest.approx(pfit*math.log10(2),rel=1e-13,abs=1e-15)
    raw,norm,gap=nw('AA','AC',symbols,matrix)
    assert norm==-4.25 and norm!=-(2.5+4.75-2*raw)/2


def test_requested_target_filter_counts_puncture_and_unsupported():
    out=run(['ACGT'],['ACGT','ACGA','TTTT'],requests=[[1]],exclude=True,roles=[[0,1],[1],[2]],other=[1,1,0],params=SearchParams(max_subs=1,max_total_edits=1),cap=0)
    assert out['candidate_counts'].tolist()==[2]
    assert out['nw_pairs'].tolist()==[1]
    assert out['counts'].tolist()==[1] and out['aa_counts'].tolist()==[1]
    assert out['informativeness'][0]==pytest.approx(math.log10(2)) # N=4roles−2; k=1−1
    assert out['best_ref_ids'].tolist()==[1]
    unsupported=run(['ACGT'],['ACGT','ACGA'],groups=[[-1],[0]],requests=[[1]])
    assert np.isnan(unsupported['scores'][0]) and unsupported['counts'][0]==1
    assert unsupported['nw_pairs'][0]==0 and unsupported['query_supported'][0]==0


def candidate(q,r):
    todo=[(0,0,0,0)];seen=set(todo)
    while todo:
        i,j,s,g=todo.pop()
        if i==len(q) and j==len(r):return True
        next_states=[]
        if i<len(q) and j<len(r):next_states.append((i+1,j+1,s+(q[i]!=r[j]),g))
        if i<len(q):next_states.append((i+1,j,s,g+1))
        if j<len(r):next_states.append((i,j+1,s,g+1))
        for state in next_states:
            if state[2]<=4 and state[3]<=2 and state[2]+state[3]<=6 and state not in seen:
                seen.add(state);todo.append(state)
    return False


def test_aggregate_indel_existence_not_chosen_broad_alignment():
    q='ACGTTGCAAGTCCTAGGCAT';refs=[q[1:]+q[:1],q[2:]+q[:2]]
    assert candidate(q,refs[0]) and not candidate(q,refs[1])
    out=run([q],refs,requests=[[0,1]])
    assert out['counts'].tolist()==[1,0] and out['candidate_counts'].tolist()==[1]
    assert out['nw_pairs'].tolist()==[1] and out['gap_counts'][0]==2
    index=Index.build(refs,alphabet='nt');p=SearchParams(max_subs=4,max_ins=2,max_dels=2,max_total_edits=6)
    assert sum(index.edit_histogram_batch([q],p,1,False,max_indels=2)[0])==1
    assert sum(index.edit_histogram_batch([q],p,1,False)[0])==2
    p=SearchParams(engine='seqtrie',max_subs=4,max_ins=2,max_dels=2,max_total_edits=6)
    with pytest.raises(ValueError,match='seqtm'):
        index.edit_histogram_batch([q],p,1,False,max_indels=2)


@pytest.mark.parametrize('threads',[1,2])
def test_ties_and_serial_parallel_equality(threads,require_thread_budget):
    refs=['ACGA','ACGC','TTTT'];queries=['ACGT','ACGG']
    expected=run(queries,refs,threads=1)
    actual=run(queries,refs,threads=threads)
    for key in expected:np.testing.assert_equal(actual[key],expected[key])


@pytest.mark.parametrize('updates',[
    {'query_target_offsets':[1,1]}, {'query_target_offsets':[0,3]}, {'query_target_offsets':[0,2],'query_target_ids':[0,0]},
    {'reference_target_offsets':[0,2]}, {'reference_target_ids':[99,1]}, {'reference_class_ids':[0,0]},
    {'reference_other_counts':[99,0]}, {'reference_multiplicities':[0,1]}, {'coefficients':[0.,float('nan'),0.]},
    {'group_offsets':[[[float('inf')]],[[0.]]]}, {'query_group_ids':[[1],[0]]}, {'reference_group_ids':[[-2,0],[0,0]]},
])
def test_boundary_validation_rejects_inconsistent_numeric_inputs(updates):
    with pytest.raises((ValueError,TypeError)):
        run(['ACGT'],['ACGA','TTTT'],requests=[[0]],updates=updates)


@pytest.mark.parametrize('matrix',[[[0.]*4 for _ in range(4)],[[1.]*3 for _ in range(4)],[[float('nan')]*4 for _ in range(4)]])
def test_invalid_matrix_or_selection_gap(matrix):
    with pytest.raises(ValueError):run(['ACGT'],['ACGA','TTTT'],matrix=matrix)


def test_counter_duplicate_symbols_and_overflow_errors():
    with pytest.raises(ValueError):run(['ACGT'],['ACGT','acgt'])
    with pytest.raises(ValueError):run(['ACGT'],['ACGA','TTTT'],symbols='AACG')
    with pytest.raises(ValueError):run(['ACGT'],['ACGA','TTTT'],symbols='ACG')
    with pytest.raises(ValueError):run(['ACGT'],['ACGA','TTTT'],cap=-1)
    with pytest.raises(OverflowError):run(['ACGT'],['ACGA','TTTT'],updates={'group_offsets':[[[1e308]],[[1e308]]]})
    index=Index.build(['ACGT'],alphabet='nt');params=SearchParams(max_subs=1,max_total_edits=1)
    with pytest.raises(ValueError):index.edit_histogram_batch(['ACGT'],params,max_indels=-2)


def test_random_fractional_nw_and_all_candidate_memberships():
    import random
    rng=random.Random(17)
    refs=['','A','AC','GT','ACG','TGC','TTTT','ACGT','ACGTT','GGACGT']
    queries=refs+['AGC','CCTA','GACGT']
    matrix=[[rng.randint(-8,9)/7 for _ in 'ACGT'] for _ in 'ACGT']
    out=run(queries,refs,matrix=matrix)
    index=Index.build(refs,alphabet='nt')
    params=SearchParams(max_subs=4,max_ins=2,max_dels=2,max_total_edits=6)
    capped=[sum(row) for row in index.edit_histogram_batch(queries,params,threads=1,max_indels=2)]
    rectangles=[]
    for ins,dels in [(2,0),(0,2),(1,1),(1,0),(0,1)]:
        p=SearchParams(max_subs=4,max_ins=ins,max_dels=dels,max_total_edits=6)
        rectangles.append(np.asarray([sum(row) for row in index.edit_histogram_batch(queries,p,threads=1)]))
    np.testing.assert_equal(capped,rectangles[0]+rectangles[1]+rectangles[2]-rectangles[3]-rectangles[4])
    np.testing.assert_equal(capped,out['candidate_counts'])
    for i,q in enumerate(queries):
        expected_candidates=sum(candidate(q,r) for r in refs)
        assert out['candidate_counts'][i]==expected_candidates
        for j,r in enumerate(refs):
            at=i*len(refs)+j
            if not candidate(q,r):assert out['counts'][at]==0 and out['scores'][at]==0
            else:
                raw,norm,gaps=nw(q,r,'ACGT',matrix)
                assert out['raw_similarity'][at]==pytest.approx(raw,abs=1e-12)
                assert out['normalized_similarity'][at]==pytest.approx(norm,abs=1e-12)
                assert out['gap_counts'][at]==gaps
    import gc
    gc.collect();assert not out['scores'].flags.writeable


def test_redundant_aggregate_cap_preserves_unit_count_fastpath_results():
    index=Index.build(['ACGT','ACGA','ACGG','ACGT','TTTT'],alphabet='nt')
    queries=['ACGT','TTTA']
    for ins,dels,total in [(0,0,1),(1,0,2),(1,1,2)]:
        p=SearchParams(max_subs=1,max_ins=ins,max_dels=dels,max_total_edits=total)
        for exclude in (False,True):
            expected=index.edit_histogram_batch(queries,p,threads=1,exclude_exact=exclude)
            for cap in (ins+dels,ins+dels+1):
                assert index.edit_histogram_batch(queries,p,threads=1,exclude_exact=exclude,max_indels=cap)==expected


@pytest.mark.parametrize('threads',[1,2])
@pytest.mark.parametrize('exclude',[False,True])
def test_candidate_backend_preserves_trie_fallback_buffers(threads,exclude,require_thread_budget):
    import random
    rng=random.Random(41)
    queries=['','A']+[''.join(rng.choice('ACGT') for _ in range(n)) for n in (8,20,63,64)]
    refs=list(dict.fromkeys(queries+[q[1:]+q[:1] for q in queries]+[q+'A' for q in queries]))
    target_ids=[i%3 for i in range(len(refs))]
    updates={'target_class_ids':[0,0,0],'reference_multiplicities':[i%4+1 for i in range(len(refs))],
        'reference_group_ids':[[(-1 if i%5==0 else 0) for i in range(len(refs))],[0]*len(refs)]}
    kwargs=dict(target_ids=target_ids,roles=[[0] for _ in refs],requests=[[i%3] for i in range(len(queries))],
        groups=[[-1,0,0,0,0,0],[0]*len(queries)],exclude=exclude,threads=threads,updates=updates)
    bounded=run(queries,refs,params=SearchParams(max_subs=4,max_ins=2,max_dels=2,max_total_edits=6),**kwargs)
    # Total7 is redundant: substitutions<=4 plus aggregate indels<=2 implies total<=6.
    # It deliberately exercises the existing trie fallback with exactly the same set.
    fallback=run(queries,refs,params=SearchParams(max_subs=4,max_ins=2,max_dels=2,max_total_edits=7),**kwargs)
    for key in bounded:
        assert bounded[key].dtype==fallback[key].dtype
        assert bounded[key].tobytes()==fallback[key].tobytes(),key
