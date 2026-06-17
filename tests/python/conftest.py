import pytest
import seqtree


@pytest.fixture
def cdr3_db():
    return [
        "CASSLAPGATNEKLFF",
        "CASSLELGATNEKLFF",
        "CASSPQGATNEKLFF",
        "CASSLAPGATNEKLF",
        "WASSLAPGATNEKLFF",
    ]


@pytest.fixture
def idx(cdr3_db):
    return seqtree.Index.build(cdr3_db, alphabet="aa")
