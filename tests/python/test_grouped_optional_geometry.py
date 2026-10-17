"""Draft regression tests for grouped optional values under nanobind's supported floor."""
import pytest
from seqtree import gapblock as gb, SubstitutionMatrix


def metadata():
    return dict(reference_target_offsets=[0, 1], reference_target_ids=[0],
        target_peptide_ids=[0], reference_exclusive_peptide_ids=[0],
        query_gene_ids=[0], reference_gene_ids=[0], threads=1)


def feature():
    return gb.PairedFeatureOffsets([[0, 0]], [[1, 1]], [[0, 3], [3, 0]])


@pytest.mark.parametrize('matrix', [None, SubstitutionMatrix.blosum62()])
@pytest.mark.parametrize('paired,alpha,beta', [
    (False, False, False), (False, True, False),
    (True, False, False), (True, True, False),
    (True, False, True), (True, True, True)])
def test_optional_features_default_and_explicit_lanes(matrix, paired, alpha, beta):
    opts = dict(matrix=matrix, **metadata())
    if paired:
        opts.update(feature_offsets_alpha=feature() if alpha else None,
                    feature_offsets_beta=feature() if beta else None)
        batch = gb.grouped_paired_sum_count_batch(['AA'], ['CC'], ['AA'], ['CC'],
                                                [0, 2, 3, 5, 6], **opts)
    else:
        opts['feature_offsets'] = feature() if alpha else None
        batch = gb.grouped_count_batch(['AA'], ['AA'], [0, 2, 3, 5, 6], **opts)
    chunk = next(batch)
    offset = 3 * (alpha + beta)
    expected = [int(offset <= cut) for cut in [0, 2, 3, 5, 6]]
    assert memoryview(chunk.counts).tolist() == [[expected] * 4]
    assert memoryview(chunk.sizes).tolist() == [[1] * 4]
    assert memoryview(chunk.exact_counts).tolist() == [[1]]


@pytest.mark.parametrize('self_ids', [None, [-1]])
def test_sparse_optional_reference_ids_default_and_explicit(self_ids):
    chunk = next(gb.grouped_count_batch(['AA'], ['AA'], [0],
        return_hits=True, query_target_ids=[0], query_reference_ids=self_ids,
        matrix=SubstitutionMatrix.blosum62(), **metadata()))
    assert memoryview(chunk.hits).tolist() == [[0, 0, 0]]
