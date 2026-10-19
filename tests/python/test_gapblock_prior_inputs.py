"""Prior preparation must preserve exact scores without evaluating absent lengths."""
import pytest
from seqtree import gapblock as gb


@pytest.mark.parametrize('mode', ['matrix', 'count', 'paired', 'lane-priors'])
def test_prior_preparation_uses_input_lengths_and_matches_scalar_scores(mode):
    qa, ra = ['AC', 'A' * 40], ['ACA', 'AC']
    qb, rb = ['A', 'AC'], ['AC', 'A']
    seen, seen_beta = set(), set()
    def value(i, d, m):
        return (i + 2*d + m) % 7
    def prior(i, d, m):
        seen.add((m, m-d, i))
        return value(i, d, m)
    def beta(i, d, m):
        seen_beta.add((m, m-d, i))
        return value(i, d, m) + 1
    options = dict(gap_open=3, gap_extend=2, threads=1)
    def scores(q, r, p):
        return [[gb.gapblock_score(a, b, gap_open=3, gap_extend=2,
                    gap_prior=p)[0] for b in r] for a in q]
    expected = scores(qa, ra, value)
    cuts = [0, 10, 100]
    lengths = {len(s) for s in qa + ra}
    if mode == 'matrix':
        actual = memoryview(gb.score_matrix(qa, ra, gap_prior=prior, **options)).tolist()
    elif mode == 'count':
        actual = gb.count_batch(qa, ra, [cuts]*len(qa), gap_prior=prior, **options)
        expected = [[sum(v <= c for v in row) for c in cuts] for row in expected]
    else:
        lanes = mode == 'lane-priors'
        other = scores(qb, rb, lambda i, d, m: value(i, d, m) + lanes)
        expected = [[a+b for a,b in zip(x,y)] for x,y in zip(expected,other)]
        expected = [[[sum(v <= c for v in row) for c in cuts]]*4 for row in expected]
        priors = dict(gap_prior_alpha=prior, gap_prior_beta=beta) if lanes else dict(gap_prior=prior)
        chunks = list(gb.grouped_paired_sum_count_batch(qa, qb, ra, rb, cuts,
            reference_target_offsets=[0, 1, 2], reference_target_ids=[0, 0],
            target_peptide_ids=[0], reference_exclusive_peptide_ids=[0, 0],
            query_gene_ids=[0, 0], reference_gene_ids=[0, 0], **options, **priors))
        actual = [row for c in chunks for row in memoryview(c.counts).tolist()]
        if not lanes:
            lengths.update(len(s) for s in qb + rb)
    assert actual == expected
    assert seen and all(m in lengths and n in lengths and 0 <= i <= n for m, n, i in seen)
    if seen_beta:
        beta_lengths = {len(s) for s in qb + rb}
        assert all(m in beta_lengths and n in beta_lengths and 0 <= i <= n for m, n, i in seen_beta)


@pytest.mark.parametrize('queries,refs', [([], []), (['AC'], ['AC']), ([], ['AC']), (['AC'], [])])
def test_prior_is_not_evaluated_when_no_length_difference_is_used(queries, refs):
    def unused(i, d, m):
        raise AssertionError('no gap-prior cell can be read')
    matrix = gb.score_matrix(queries, refs, gap_prior=unused, threads=1)
    assert matrix.shape == (len(queries), len(refs))
    assert gb.count_batch(queries, refs, [[0, 1]]*len(queries), gap_prior=unused,
                          threads=1) == [[len(refs)]*2 for _ in queries]
