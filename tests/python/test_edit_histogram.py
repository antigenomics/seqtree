"""Native histogram must preserve the existing search's edit counts and multiplicity."""
import pytest
from seqtree import Index, SearchParams


def test_edit_histogram_matches_search():
    index = Index.build(['ACDE', 'ACDE', 'ACDF', 'ACD', 'ACGDE', 'GGGG'])
    queries = ['ACDE', 'ACD', 'GGGG', '']
    params = SearchParams(max_subs=2, max_ins=2, max_dels=2, max_total_edits=2)
    for exclude in (False, True):
        expected = []
        for row in index.search_batch(queries, params, threads=1):
            bins = [0, 0, 0]
            for hit in row:
                total = hit.n_subs + hit.n_ins + hit.n_dels
                if not (exclude and total == 0):
                    bins[total] += 1
            expected.append(bins)
        assert index.edit_histogram_batch(queries, params, threads=1, exclude_exact=exclude) == expected
        for threads in (4, 0, -1):
            assert index.edit_histogram_batch(queries, params, threads=threads, exclude_exact=exclude) == expected
    assert Index.build([]).edit_histogram_batch(['ACDE'], params, threads=1) == [[0, 0, 0]]
    assert index.edit_histogram_batch([], params, threads=1) == []


@pytest.mark.parametrize('params', [SearchParams(max_subs=2),
    SearchParams(max_total_edits=2, engine='seqtrie'),
    SearchParams(max_subs=2, max_total_edits=2, mode='top')])
def test_edit_histogram_rejects_incomplete_counts(params):
    with pytest.raises(ValueError):
        Index.build(['ACDE']).edit_histogram_batch(['ACDE'], params, threads=1)
