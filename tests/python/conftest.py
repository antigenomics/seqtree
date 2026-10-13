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


@pytest.fixture
def require_thread_budget(threads):
    """Run explicit parallel cases only when the native allocation permits them."""
    index = seqtree.Index.build(["A"])
    try:
        index.search_batch(["A"] * threads, seqtree.SearchParams(), threads=threads)
    except ValueError as error:
        if "available CPU budget" not in str(error):
            raise
        pytest.skip(f"requested {threads} threads exceed this runner's native CPU budget")
