"""Additive generic factored terms and independent paired priors."""
import pytest
from seqtree import gapblock, SubstitutionMatrix


def rows(buffer):
    return memoryview(buffer).tolist()


def metadata():
    return dict(reference_target_offsets=[0,1,2,4],reference_target_ids=[0,1,0,1],
        target_peptide_ids=[0,1],reference_exclusive_peptide_ids=[0,1,-1],
        query_gene_ids=[0,1],reference_gene_ids=[0,1,0])


def feature(cost, weight=1, scale=200):
    return gapblock.PairedFeatureOffsets([[0,0],[1,1]],[[1,1],[0,0],[1,1]],
        [[0,cost],[cost,0]],weight=weight,scale=scale)


@pytest.mark.parametrize('withheld',[False,True])
@pytest.mark.parametrize('paired',[False,True])
@pytest.mark.parametrize('threads',[1,4])
def test_independent_terms_and_lane_priors_against_explicit_expression(paired,threads,withheld):
    q,r=['AAAA','AAAR'],['AAAA','AARR','AAAAAA']
    qb,rb=['RRRR','RRAA'],['RRRR','AARRR','AAAA']
    matrix=SubstitutionMatrix.blosum62(); cuts=list(range(0,900,20))
    pa=gapblock.frame_prior(0,26);pb=gapblock.frame_prior(3,26)
    opts=dict(matrix=matrix,gap_open=13,gap_extend=13,threads=threads,
        exclude_exact=True,chunk_size=1,**metadata())
    folds=[-1,0] if withheld else [-1,-1]
    if withheld:
        opts.update(query_withheld_peptide_ids=folds,reference_peptide_offsets=[0,1,2,4],
                    reference_peptide_ids=[0,1,0,1])
    primary=feature(1,weight=.5)
    extras=[feature(1,weight=.5),feature(1,weight=.5)]
    if paired:
        batch=gapblock.grouped_paired_sum_count_batch(q,qb,r,rb,cuts,**opts,
            gap_prior_alpha=pa,gap_prior_beta=pb,feature_offsets_alpha=primary,
            additional_feature_offsets_alpha=extras,additional_feature_offsets_beta=[feature(3)])
    else:
        batch=gapblock.grouped_count_batch(q,r,cuts,**opts,gap_prior=pa,
            feature_offsets=primary,additional_feature_offsets=extras)
    base=gapblock.score_matrix(q,r,matrix,gap_open=13,gap_extend=13,gap_prior=pa)
    other=gapblock.score_matrix(qb,rb,matrix,gap_open=13,gap_extend=13,gap_prior=pb) if paired else None
    # Three half-unit terms round separately to zero. Beta's integer term is600.
    scores=[[base[i,j]+(other[i,j]+600*int(i%2 != (1,0,1)[j]) if paired else 0)
             for j in range(3)] for i in range(2)]
    expected=[];sizes=[];exact=[]
    members=[[0],[1],[0,1]]
    for i in range(2):
        channels=[[] for _ in range(7)];ex=[0,0]
        for j,original_labels in enumerate(members):
            labels=[t for t in original_labels if t!=folds[i]]
            if not labels:continue
            identical=q[i]==r[j] and (not paired or qb[i]==rb[j])
            if identical:
                for t in labels:ex[t]+=1
                continue
            channels[0].append(scores[i][j])
            for t in labels:
                channels[1+t].append(scores[i][j])
                if i%2==j%2:channels[3+t].append(scores[i][j])
            if len(labels)==1:channels[5+labels[0]].append(scores[i][j])
        expected.append([[sum(score<=cut for score in row) for cut in cuts] for row in channels])
        sizes.append([len(row) for row in channels]);exact.append(ex)
    chunks=list(batch)
    assert [row for chunk in chunks for row in rows(chunk.counts)]==expected
    assert [row for chunk in chunks for row in rows(chunk.sizes)]==sizes
    assert [row for chunk in chunks for row in rows(chunk.exact_counts)]==exact


@pytest.mark.parametrize('lane',['alpha','beta','both'])
def test_paired_prior_overrides_inherit_shared_independently(lane):
    q,r=['AAAA','AAAR'],['AAAA','AARR','AAAAAA']
    qb,rb=['RRRR','RRAA'],['RRRR','AARRR','AAAA']
    matrix=SubstitutionMatrix.blosum62();shared=gapblock.frame_prior(1,26)
    pa=gapblock.frame_prior(0,26);pb=gapblock.frame_prior(3,26)
    kwargs={}
    if lane in ('alpha','both'):kwargs['gap_prior_alpha']=pa
    if lane in ('beta','both'):kwargs['gap_prior_beta']=pb
    a=gapblock.score_matrix(q,r,matrix,gap_open=13,gap_extend=13,
        gap_prior=pa if 'gap_prior_alpha' in kwargs else shared)
    b=gapblock.score_matrix(qb,rb,matrix,gap_open=13,gap_extend=13,
        gap_prior=pb if 'gap_prior_beta' in kwargs else shared)
    cuts=list(range(0,300,5))
    batch=gapblock.grouped_paired_sum_count_batch(q,qb,r,rb,cuts,matrix=matrix,
        gap_open=13,gap_extend=13,gap_prior=shared,**kwargs,**metadata())
    got=rows(next(batch).counts)
    assert [row[0] for row in got]==[[sum(a[i,j]+b[i,j]<=cut for j in range(3))
        for cut in cuts] for i in range(2)]


@pytest.mark.parametrize('bad',[
    [object()], [gapblock.PairedFeatureOffsets([],[],[[0]])],
    [feature(-1)], [feature(1,weight=float('nan'))], [feature(1,scale=-1)],
    [feature(1,scale=2**31-1)],
])
def test_additional_term_validation(bad):
    with pytest.raises((ValueError,OverflowError)):
        gapblock.grouped_count_batch(['AAAA','AAAR'],['AAAA','AARR','AAAAAA'],[0],
            additional_feature_offsets=bad,**metadata())


def test_native_term_and_prior_validation(monkeypatch):
    original=gapblock._additional_geometry
    def corrupt(geometry,*args):
        original(geometry,*args)
        geometry.additional_features=[gapblock.PairedFeatureOffsets([],[],[[0]])]
    monkeypatch.setattr(gapblock,'_additional_geometry',corrupt)
    with pytest.raises(ValueError,match='axes'):
        gapblock.grouped_count_batch(['AAAA','AAAR'],['AAAA','AARR','AAAAAA'],[0],**metadata())
    monkeypatch.setattr(gapblock,'_additional_geometry',original)
    original_prior=gapblock._lane_prior
    def corrupt_prior(geometry,*args):
        original_prior(geometry,*args);geometry.gap_prior=[0]
    monkeypatch.setattr(gapblock,'_lane_prior',corrupt_prior)
    with pytest.raises(ValueError,match='prior table'):
        gapblock.grouped_paired_sum_count_batch(['AAAA','AAAR'],['AAAA','AAAR'],
            ['AAAA','AARR','AAAAAA'],['AAAA','AARR','AAAAAA'],[0],**metadata())


@pytest.mark.parametrize('mode',['gapblock','global','local'])
def test_extra_terms_remain_outside_lane_quantization(mode):
    q,r=['AAAA'],['RRRR'];matrix=SubstitutionMatrix.blosum62()
    opts=dict(matrix=matrix,gap_open=13,gap_extend=13)
    if mode=='gapblock':
        raw=gapblock.score_matrix(q,r,**opts)[0,0]
    else:
        from seqtree import pairwise
        raw=max(0,sum(matrix.similarity(a,a) for a in q[0]+r[0])-
            2*pairwise.score(q[0],r[0],matrix,mode=mode,gap_open=13,gap_extend=13))*3
        opts.update(alignment_mode=mode,alignment_scale=3)
    expected=raw//5*5+raw//7*7+1000
    def term(cost):
        return gapblock.PairedFeatureOffsets([[0,0]],[[1,1]],[[0,cost],[cost,0]],scale=200)
    batch=gapblock.grouped_paired_sum_count_batch(q,q,r,r,[expected-1,expected],**opts,
        score_quantum_alpha=5,score_quantum_beta=7,additional_feature_offsets_alpha=[term(2)],
        additional_feature_offsets_beta=[term(3)],reference_target_offsets=[0,1],
        reference_target_ids=[0],target_peptide_ids=[0],reference_exclusive_peptide_ids=[0],
        query_gene_ids=[0],reference_gene_ids=[0])
    assert rows(next(batch).counts)==[[[0,1]]*4]
