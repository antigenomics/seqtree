"""Bounded affine grouped counts against the existing pairwise-score oracle."""
import pytest
from seqtree import SubstitutionMatrix, gapblock, pairwise


def buf(value):
    return memoryview(value).tolist()


@pytest.mark.parametrize('mode', ['local', 'global'])
@pytest.mark.parametrize('paired', [False, True])
@pytest.mark.parametrize('withheld', [False, True])
@pytest.mark.parametrize('threads', [1, 4])
def test_all_curves_sizes_exact_against_independent_raw_diagonal_oracle(mode, paired, withheld, threads, require_thread_budget):
    q = ['CAAAAF', 'CRRRRF', 'WAAAAW']
    r = ['CAAAAF', 'CRRAAF', 'CRRRRF', 'CAAARRAAF', 'WAAAAW']
    qb = ['CWWWWF', 'CAAAAF', 'CSSSSF']
    rb = ['CWWWWF', 'CWWAAF', 'CAAAAF', 'CNNAAF', 'CSSSSF']
    labels = [[0], [0, 1], [1], [], [1]]
    full = [[0], [0, 1], [1], [2], [1]]
    exclusive = [0, -1, 1, -1, 1]
    qgene, rgene = [0, 1, 0], [0, 0, 1, 0, 2]
    groupsq, groupsr = [0, 1, 0], [1, 0, 0, 1, 0]
    offsets = [[0, 3], [5, 0]]
    cuts = [0, 3, 100, 6000, 18000, 48000, 120000]
    matrix = SubstitutionMatrix.blosum62()
    csr = [0]; ids = []
    for row in labels:
        ids.extend(row); csr.append(len(ids))
    folds = [-1, 0, 1] if withheld else [-1] * len(q)
    common = dict(reference_target_offsets=csr, reference_target_ids=ids,
        target_peptide_ids=[0, 1], reference_exclusive_peptide_ids=exclusive,
        query_gene_ids=qgene, reference_gene_ids=rgene, threads=threads,
        exclude_exact=True, chunk_size=2, matrix=matrix, gap_open=13, gap_extend=13,
        alignment_mode=mode, alignment_scale=600, trim_left=1, trim_right=1)
    if withheld:
        csrfull = [0]; fullids = []
        for row in full:
            fullids.extend(row); csrfull.append(len(fullids))
        common.update(query_withheld_peptide_ids=folds, reference_peptide_offsets=csrfull,
                      reference_peptide_ids=fullids)
    geometry = dict(query_group_ids=groupsq, reference_group_ids=groupsr, group_distances=offsets,
        feature_offsets=gapblock.PairedFeatureOffsets([[0, 0]] * len(q),
            [[1, 1]] * len(r), [[0, 2], [2, 0]], scale=200))
    if paired:
        opts = {key + suffix: value for suffix in ('_alpha', '_beta') for key, value in geometry.items()}
        batch = gapblock.grouped_paired_sum_count_batch(q, qb, r, rb, cuts, **common, **opts)
    else:
        batch = gapblock.grouped_count_batch(q, r, cuts, **common, **geometry)
    # Dense test oracle only: production calls the grouped population once.
    def distances(left, right):
        l, rr = [s[1:-1] for s in left], [s[1:-1] for s in right]
        sim = pairwise.score_matrix(l, rr, matrix, mode=mode, gap_open=13, gap_extend=13, threads=1)
        diagl = [sum(matrix.similarity(a, a) for a in s) for s in l]
        diagr = [sum(matrix.similarity(a, a) for a in s) for s in rr]
        return [[600 * max(0, diagl[i] + diagr[j] - 2 * sim[i, j])
                 + offsets[groupsq[i]][groupsr[j]] + 400 for j in range(len(right))]
                for i in range(len(left))]
    scores = distances(q, r)
    if paired:
        other = distances(qb, rb)
        scores = [[x + y for x, y in zip(a, b)] for a, b in zip(scores, other)]
    expected, sizes, exact = [], [], []
    for i in range(len(q)):
        channels = [[] for _ in range(7)]; exactrow = [0, 0]
        for j in range(len(r)):
            surviving = set(full[j]) - {folds[i]}
            if not surviving: continue
            remaining = [t for t in labels[j] if t != folds[i]]
            identical = q[i] == r[j] and (not paired or qb[i] == rb[j])
            if identical:
                for t in remaining: exactrow[t] += 1
                continue
            channels[0].append(scores[i][j])
            for t in remaining:
                channels[1 + t].append(scores[i][j])
                if qgene[i] == rgene[j]: channels[3 + t].append(scores[i][j])
            pep = next(iter(surviving)) if len(surviving) == 1 else -1
            if pep in (0, 1): channels[5 + pep].append(scores[i][j])
        expected.append([[sum(x <= cut for x in row) for cut in cuts] for row in channels])
        sizes.append([len(row) for row in channels]); exact.append(exactrow)
    chunks = list(batch)
    assert [row for chunk in chunks for row in buf(chunk.counts)] == expected
    assert [row for chunk in chunks for row in buf(chunk.sizes)] == sizes
    assert [row for chunk in chunks for row in buf(chunk.exact_counts)] == exact


def test_diagonal_sums_are_not_self_alignment_scores():
    matrix = SubstitutionMatrix.from_similarity([[-1, -3, -3, -3], [-3, 7, -3, -3],
        [-3, -3, 2, -3], [-3, -3, -3, 3]])
    # The NT codec is A,C,G,T. Local self-AA alignment is0, raw diagonal sum is−2.
    assert pairwise.score('AA', 'AA', matrix, alphabet='nt', mode='local') == 0
    batch = gapblock.grouped_count_batch(['AA'], ['CC'], [12, 13, 14], alphabet='nt',
        matrix=matrix, alignment_mode='local', gap_open=13, gap_extend=13,
        reference_target_offsets=[0, 1], reference_target_ids=[0], target_peptide_ids=[0],
        reference_exclusive_peptide_ids=[0], query_gene_ids=[0], reference_gene_ids=[0])
    assert buf(next(batch).counts) == [[[1, 1, 1]] * 4]


def one_population(q=('CAAF',), r=('CRRF',), **options):
    return gapblock.grouped_count_batch(q, r, [0, 1, 600, 18000, 100000],
        matrix=SubstitutionMatrix.blosum62(), gap_open=13, gap_extend=13,
        reference_target_offsets=list(range(len(r)+1)), reference_target_ids=[0]*len(r),
        target_peptide_ids=[0], reference_exclusive_peptide_ids=[0]*len(r),
        query_gene_ids=[0]*len(q), reference_gene_ids=[0]*len(r), **options)


@pytest.mark.parametrize('options', [
    {'alignment_mode': 'nw'}, {'alignment_mode': None},
    {'alignment_scale': 2}, {'trim_left': 1},
    {'alignment_mode': 'global', 'alignment_scale': 0},
    {'alignment_mode': 'global', 'alignment_scale': True},
    {'alignment_mode': 'global', 'trim_left': -1},
    {'alignment_mode': 'global', 'trim_right': 1.5},
    {'alignment_mode': 'local', 'trim_left': 5},
    {'alignment_mode': 'global', 'query_position_weights': [[1]*4]},
    {'alignment_mode': 'local', 'position_weights_by_length': {4: [1]*4}},
    {'alignment_mode': 'global', 'gap_prior': lambda i,d,m: 0},
    {'alignment_mode': 'global', 'alignment_scale': 2**31-1},
])
def test_affine_option_validation(options):
    with pytest.raises((ValueError, OverflowError)):
        one_population(**options)


@pytest.mark.parametrize('mode', ['global', 'local'])
def test_empty_axes_and_fully_trimmed_sequences(mode):
    assert list(one_population(q=(), alignment_mode=mode)) == []
    empty = next(one_population(r=(), alignment_mode=mode))
    assert buf(empty.sizes) == [[0]*4]
    assert buf(empty.counts) == [[[0]*5]*4]
    trimmed = next(one_population(alignment_mode=mode, trim_left=2, trim_right=2))
    assert buf(trimmed.counts) == [[[1]*5]*4]
    assert buf(trimmed.exact_counts) == [[0]]


def test_affine_quantum_floors_scaled_junction_before_offset():
    matrix = SubstitutionMatrix.blosum62()
    q, r = 'CAAF', 'CRRF'
    raw = max(0, sum(matrix.similarity(a,a) for a in q) +
        sum(matrix.similarity(a,a) for a in r) -
        2*pairwise.score(q, r, matrix, mode='global', gap_open=13, gap_extend=13))
    expected = (raw*3)//100*100+7
    batch = gapblock.grouped_count_batch([q], [r], [expected-1, expected],
        matrix=matrix, gap_open=13, gap_extend=13, alignment_mode='global', alignment_scale=3,
        score_quantum=100, query_group_ids=[0], reference_group_ids=[0], group_distances=[[7]],
        reference_target_offsets=[0, 1], reference_target_ids=[0], target_peptide_ids=[0],
        reference_exclusive_peptide_ids=[0], query_gene_ids=[0], reference_gene_ids=[0])
    assert buf(next(batch).counts) == [[[0, 1]]*4]


@pytest.mark.parametrize('mode', ['global', 'local'])
def test_exact_evidence_is_retained_and_no_hit_keeps_sizes(mode):
    retained = next(one_population(r=('CAAF',), alignment_mode=mode))
    assert buf(retained.exact_counts) == [[1]]
    assert buf(retained.counts) == [[[1]*5]*4]
    nohit = gapblock.grouped_count_batch(['CAAF'], ['CRRF'], [-1],
        matrix=SubstitutionMatrix.blosum62(), gap_open=13, gap_extend=13,
        alignment_mode=mode, reference_target_offsets=[0,1], reference_target_ids=[0],
        target_peptide_ids=[0], reference_exclusive_peptide_ids=[0],
        query_gene_ids=[0], reference_gene_ids=[0])
    chunk = next(nohit)
    assert buf(chunk.counts) == [[[0]]*4]
    assert buf(chunk.sizes) == [[1]*4]


@pytest.mark.parametrize('field,value', [('alignment_mode', 'invalid'),
    ('alignment_scale', 0), ('trim_left', 2**32-1)])
def test_native_affine_validation_rechecks_python_snapshot(monkeypatch, field, value):
    original = gapblock._alignment_options
    def corrupt(geometry, *args):
        original(geometry, *args)
        setattr(geometry, field, value)
    monkeypatch.setattr(gapblock, '_alignment_options', corrupt)
    with pytest.raises(ValueError):
        one_population(alignment_mode='local')
