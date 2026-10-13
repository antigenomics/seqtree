"""Exact installed-API oracle for the private unit histogram fast path."""
import random
import pytest
from seqtree import Index, SearchParams, SubstitutionMatrix, PositionalMatrix


def _expected(index, queries, params, exclude):
    out=[]
    for hits in index.search_batch(queries,params,threads=1):
        bins=[0]*(params.max_total_edits+1)
        for hit in hits:
            total=hit.n_subs+hit.n_ins+hit.n_dels
            if not exclude or total:bins[total]+=1
        out.append(bins)
    return out


def test_exhaustive_binary_and_duplicate_case_counts():
    sequences=['']+[''.join('A' if mask&(1<<i) else 'C' for i in range(n))
                       for n in range(1,7) for mask in range(1<<n)]
    index=Index.build(sequences+['AC','AC','ac'])
    for cap in range(1,6):
        params=SearchParams(max_subs=cap,max_ins=2,max_dels=2,max_total_edits=cap)
        for exclude in (False,True):
            expected=_expected(index,sequences,params,exclude)
            for threads in (1,2):
                assert index.edit_histogram_batch(sequences,params,threads,exclude)==expected


def test_boundaries_iupac_and_existing_engine_fallbacks():
    rng=random.Random(17);aa='ACDEFGHIKLMNPQRSTVWY'
    sequences=[''.join(rng.choice(aa) for _ in range(n)) for n in (1,2,10,20,63,64,65)]
    index=Index.build(sequences+[s+'AC' for s in sequences]+['acde','ACDE'])
    queries=sequences+['','ACDE','EFGHIKLM','ACDEFGHIK','INVALID!']
    for cap in range(1,6):
        params=SearchParams(max_subs=cap,max_ins=2,max_dels=2,max_total_edits=cap)
        assert index.edit_histogram_batch(queries[:-1],params,threads=2,exclude_exact=True)==_expected(index,queries[:-1],params,True)
        with pytest.raises(ValueError):index.edit_histogram_batch(queries,params,threads=2,exclude_exact=True)
    q=['ACDE','ACD'];i=Index.build(['ACDE','ACD','ACGDE','ACDF'])
    for kwargs in ({'max_ins':1},{'gap_open':0},{'max_penalty':2},{'matrix':SubstitutionMatrix.unit(24)}):
        args=dict(max_subs=5,max_ins=2,max_dels=2,max_total_edits=5);args.update(kwargs)
        p=SearchParams(**args)
        assert i.edit_histogram_batch(q,p,threads=2,exclude_exact=True)==_expected(i,q,p,True)
    p=SearchParams(max_subs=5,max_ins=0,max_dels=0,max_total_edits=5)
    p.pos_matrix=PositionalMatrix.from_weights(SubstitutionMatrix.unit(24),[0,1,1,0])
    assert i.edit_histogram_batch(['ACDE'],p,threads=2,exclude_exact=True)==_expected(i,['ACDE'],p,True)
    n=Index.build(['ACGTN','RYACGT','ryacgt','RYYACGT'],'iupac')
    p=SearchParams(max_subs=5,max_ins=2,max_dels=2,max_total_edits=5)
    assert n.edit_histogram_batch(['acgtN','RYACGT'],p,threads=2,exclude_exact=True)==_expected(n,['acgtN','RYACGT'],p,True)


@pytest.mark.parametrize("cap", [3,4,5])
@pytest.mark.parametrize("exclude", [False,True])
def test_bounded_rejection_preserves_per_type_caps_and_dense_hits(cap,exclude):
    refs=['ACDEFGHIKLMNPQRST','acdefghiklmnpqrst',
          'ACDEFACDGHIKLMNPQRST','ACDEFGHIKLMNP',
          'ACDEFGHIKLMNPACDRST','ACDGHIKLMNPQRST',
          'AAAAAAAAAAAAAAAA','CCCCCCCCCCCCCCCC']
    index=Index.build(refs+refs[:2])
    queries=refs+['ACDEFGHIKLMNPQRS','ACDEFGHIKLMNPQACDRST']
    for limits in ({},{'max_subs':1},{'max_ins':1},{'max_dels':1},
                   {'max_ins':3},{'max_dels':3},{'max_ins':0,'max_dels':0}):
        args=dict(max_subs=cap,max_ins=2,max_dels=2,max_total_edits=cap)
        args.update(limits)
        params=SearchParams(**args)
        expected=_expected(index,queries,params,exclude)
        for threads in (1,4):
            assert index.edit_histogram_batch(queries,params,threads,exclude)==expected
