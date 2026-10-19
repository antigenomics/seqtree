"""Bounded primitive search output retains the ordinary batch's exact hits."""
import pytest
from seqtree import Index, SearchParams


@pytest.mark.parametrize('threads', [1, 2])
@pytest.mark.parametrize('weighted', [False, True])
def test_consumer_equals_full_batch_and_keeps_owned_complete_buffers(threads, weighted):
    refs = ['AAAA', 'AAAC', 'AACC', 'CCCC', 'AAA', 'AAAAA', 'AAAA']
    queries = ['AAAA', 'CCCC', 'GGGG', 'AAAC'] * 9
    index = Index.build(refs, alphabet='nt')
    params = SearchParams(max_subs=2, max_ins=1, max_dels=1, max_total_edits=2)
    if weighted:
        from seqtree import SubstitutionMatrix
        params.matrix = SubstitutionMatrix.from_similarity([[2 if i == j else -1 for j in range(4)] for i in range(4)])
        params.gap_open = params.gap_extend = 3
    expected = [[(h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels) for h in hits]
        for hits in index.search_batch(queries, params, threads)]
    blocks = []
    def consume(start, stop, buffer):
        view = memoryview(buffer)
        assert view.readonly and view.shape[1] == 6
        assert view.nbytes <= 24 * len(refs) * 2
        assert 0 < stop-start <= 2
        blocks.append((start, stop, buffer))  # owns data after callback return
    index.search_batch_consume(queries, params, consume, threads=threads,
        chunk_size=2, max_hit_bytes=24*len(refs)*2)
    observed = [[] for _ in queries]
    covered = []
    for start, stop, buffer in sorted(blocks, key=lambda b: b[0]):
        covered.extend(range(start, stop))
        for row in memoryview(buffer).tolist():
            assert start <= row[0] < stop
            observed[row[0]].append(tuple(row[1:]))
    assert covered == list(range(len(queries)))
    assert observed == expected


@pytest.mark.parametrize('base_exception', [False, True])
def test_consumer_overflow_and_callback_errors_join_workers_and_preserve_index(base_exception):
    index = Index.build(['AAAA'] * 8, alphabet='nt')
    params = SearchParams(max_subs=1)
    with pytest.raises(OverflowError, match='max_hit_bytes'):
        index.search_batch_consume(['AAAA']*20, params, lambda *_: None,
            threads=2, max_hit_bytes=24)
    class CallbackError(BaseException if base_exception else Exception):
        pass
    calls = []
    def fail(*args):
        calls.append(1)
        raise CallbackError('consumer stopped')
    with pytest.raises(CallbackError, match='consumer stopped'):
        index.search_batch_consume(['AAAA']*20, params, fail, threads=2, chunk_size=1)
    assert calls == [1]
    assert len(index.search('AAAA', params)) == 8


def test_consumer_snapshots_params_and_bounds_dense_output():
    import threading
    from seqtree import PositionalMatrix, SubstitutionMatrix
    index = Index.build(['AAAA']*100, alphabet='nt')
    params = SearchParams()
    params.pos_matrix = PositionalMatrix.from_weights(SubstitutionMatrix.unit(4), [1, 1, 1, 0])
    calls = []
    workers = set()
    def consume(start, stop, buffer):
        assert stop-start == 1 and memoryview(buffer).nbytes == 2400
        params.max_subs = 0  # mutable Python parameters cannot change this batch
        params.pos_matrix = PositionalMatrix.from_weights(SubstitutionMatrix.unit(4), [1, 1, 1, 1])
        assert memoryview(buffer).shape == (100, 6)
        calls.append(start)
        workers.add(threading.get_native_id())
    index.search_batch_consume(['AAAC']*1000, params, consume, threads=2,
        chunk_size=32, max_hit_bytes=2400)
    assert sorted(calls) == list(range(1000))
    assert len(workers) <= 2  # one native pool, not one new pool per output block
    assert index.search('AAAC', params) == []


def test_consumer_emits_empty_queries_and_native_errors_leave_index_usable():
    index = Index.build(['AAAA'], alphabet='nt')
    blocks = []
    index.search_batch_consume(['CCCC'], SearchParams(),
        lambda start, stop, buffer: blocks.append((start, stop, memoryview(buffer).shape)))
    assert blocks == [(0, 1, (0, 6))]
    with pytest.raises(ValueError, match='invalid symbol'):
        index.search_batch_consume(['AAAA', '!!!!'], SearchParams(), lambda *_: None, threads=2)
    assert len(index.search('AAAA', SearchParams())) == 1


def test_callbacks_remain_serial_when_python_releases_the_gil():
    import time
    active = [0]
    calls = []
    def consume(start, stop, buffer):
        active[0] += 1
        assert active[0] == 1
        time.sleep(.001)
        calls.append(start)
        active[0] -= 1
    Index.build(['AAAA'], alphabet='nt').search_batch_consume(
        ['AAAA']*32, SearchParams(), consume, threads=2, chunk_size=1)
    assert sorted(calls) == list(range(32)) and active == [0]


def test_consumer_peak_memory_does_not_accumulate_dense_cohort_hits():
    import json
    import os
    import subprocess
    import sys
    if sys.platform == 'win32':
        pytest.skip('resource.ru_maxrss is unavailable on Windows; buffer/pool bounds are tested above')
    code = '''import json,resource,sys
from seqtree import Index,SearchParams
index=Index.build(['AAAA']*1000,alphabet='nt')
count=[0]
def consume(start,stop,buffer):
    count[0]+=memoryview(buffer).shape[0]
index.search_batch_consume(['AAAA']*int(sys.argv[1]),SearchParams(),consume,
    threads=2,chunk_size=32,max_hit_bytes=96000)
peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
print(json.dumps({'hits':count[0],'peak':peak*(1 if sys.platform=='darwin' else 1024)}))
'''
    small, large = [json.loads(subprocess.check_output([sys.executable, '-c', code, str(n)],
        text=True, env=os.environ, timeout=30)) for n in (200, 20000)]
    assert small['hits'] == 200000 and large['hits'] == 20000000
    assert large['peak'] <= small['peak'] + 32*1024**2


@pytest.mark.parametrize('kwargs', [{'chunk_size':0}, {'max_hit_bytes':0}])
def test_consumer_rejects_invalid_output_limits(kwargs):
    with pytest.raises(ValueError):
        Index.build(['AAAA'], alphabet='nt').search_batch_consume(
            ['AAAA'], SearchParams(), lambda *_: None, **kwargs)
