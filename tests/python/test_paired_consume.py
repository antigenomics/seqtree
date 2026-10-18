"""Coordinated full-population search preserves each chain's ordinary hits."""
import pytest
from seqtree import Index, SearchParams, PositionalMatrix, SubstitutionMatrix

@pytest.mark.parametrize('threads',[1,2])
@pytest.mark.parametrize('weighted',[False,True])
def test_paired_buffers_equal_legacy(threads,weighted):
    a=Index.build(['AAAA','AAAC','AAAA','CCCC'],alphabet='nt')
    b=Index.build(['CCCC','CCCA','AAAA'],alphabet='nt')
    qa=['AAAA','GGGG','AAAC']*25;qb=['CCCC','GGGG','AAAA']*25
    p=SearchParams(max_subs=2,max_ins=1,max_dels=1,max_total_edits=2)
    if weighted:
        p.matrix=SubstitutionMatrix.unit(4);p.gap_open=p.gap_extend=2
    expected=[[[tuple(h) for h in row] for row in index.search_batch(q,p,threads)] for index,q in [(a,qa),(b,qb)]]
    blocks=[]
    a.search_paired_batch_consume(b,qa,qb,p,lambda *args:blocks.append(args),threads=threads,chunk_size=2,max_hit_bytes=24*7*2)
    actual=[[[] for _ in qa] for _ in range(2)];covered=[]
    for start,stop,*buffers in sorted(blocks):
        covered.extend(range(start,stop))
        assert 0<stop-start<=2
        assert sum(memoryview(buf).nbytes for buf in buffers)<=24*7*2
        for chain,buf in enumerate(buffers):
            assert memoryview(buf).readonly
            for qi,*hit in memoryview(buf).tolist():actual[chain][qi].append(tuple(hit))
    assert covered==list(range(len(qa))) and actual==expected

@pytest.mark.parametrize('base',[False,True])
def test_paired_failure_cleanup_and_snapshot(base):
    a=Index.build(['AAAA']*10,alphabet='nt');b=Index.build(['AAAC']*10,alphabet='nt')
    p=SearchParams(max_subs=1)
    p.pos_matrix=PositionalMatrix.from_weights(SubstitutionMatrix.unit(4),[1,1,1,0])
    calls=[]
    def callback(start,stop,x,y):
        p.max_subs=0;p.pos_matrix=PositionalMatrix.from_weights(SubstitutionMatrix.unit(4),[1,1,1,1])
        assert memoryview(x).shape[0]==10 and memoryview(y).shape[0]==10
        calls.append(start)
    a.search_paired_batch_consume(b,['AAAA']*20,['AAAA']*20,p,callback,threads=2,chunk_size=1)
    assert sorted(calls)==list(range(20))
    class Stop(BaseException if base else Exception):pass
    calls.clear()
    def fail(*_):calls.append(1);raise Stop('stop')
    with pytest.raises(Stop,match='stop'):
        a.search_paired_batch_consume(b,['AAAA']*20,['AAAA']*20,p,fail,threads=2,chunk_size=1)
    assert calls==[1] and len(a.search('AAAA',p))==10
    with pytest.raises(OverflowError,match='max_hit_bytes'):
        a.search_paired_batch_consume(b,['AAAA'],['AAAA'],p,lambda *_:None,max_hit_bytes=24)
    with pytest.raises(ValueError,match='equal lengths'):
        a.search_paired_batch_consume(b,['AAAA'],[],p,lambda *_:None)
    with pytest.raises(ValueError,match='same alphabet'):
        a.search_paired_batch_consume(Index.build(['AAAA']),[],[],p,lambda *_:None)


def test_paired_callbacks_serial_and_invalid_query_cleanup():
    import threading,time
    a=Index.build(['AAAA']*3,alphabet='nt');b=Index.build(['CCCC']*4,alphabet='nt')
    p=SearchParams(max_subs=4);active=0;workers=set()
    def consume(start,stop,x,y):
        nonlocal active
        active+=1;assert active==1
        workers.add(threading.get_native_id())
        time.sleep(.001)
        assert active==1;active-=1
    a.search_paired_batch_consume(b,['AAAA']*40,['CCCC']*40,p,consume,threads=2,chunk_size=1)
    assert len(workers)<=2
    with pytest.raises(ValueError):
        a.search_paired_batch_consume(b,['AAAA']*20,['***']*20,p,consume,threads=2,chunk_size=1)
    assert len(b.search('CCCC',p))==4


def test_paired_peak_memory_does_not_retain_first_chain_population():
    import json,os,subprocess,sys
    if sys.platform=='win32':
        pytest.skip('resource.ru_maxrss unavailable; primitive ownership and caps tested above')
    code='''import json,resource,sys
from seqtree import Index,SearchParams
a=Index.build(['AAAA']*1000,alphabet='nt');b=Index.build(['CCCC']*1000,alphabet='nt');count=[0]
def consume(start,stop,x,y):count[0]+=memoryview(x).shape[0]+memoryview(y).shape[0]
a.search_paired_batch_consume(b,['AAAA']*int(sys.argv[1]),['CCCC']*int(sys.argv[1]),SearchParams(),consume,threads=2,chunk_size=32,max_hit_bytes=96000)
print(json.dumps({'hits':count[0],'peak':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform=='darwin' else 1024)}))'''
    small,large=[json.loads(subprocess.check_output([sys.executable,'-c',code,str(n)],
        env=os.environ,text=True,timeout=30)) for n in (200,20000)]
    assert small['hits']==400000 and large['hits']==40000000
    assert large['peak']<=small['peak']+32*1024**2
