"""A cold cache shared by concurrent processes must never hand back a half-written index.

This is the pytest-xdist / Snakemake / Nextflow shape: N workers start at once against an empty
``~/.cache``, and one is still serializing a 45 MB index when the next already sees the path and
calls ``Index.load``. Before ``Index::save`` wrote through a temporary and renamed it into place,
a reader landing in that ~55 ms window got ``RuntimeError: truncated or corrupt index`` -- 10 times
out of 10.

Once warm the caches are read-only and safe; this only bites on first use.

These tests spawn real processes, because the bug lives between them: threads would share the
writer's buffers and never reproduce it.
"""
import multiprocessing as mp
import os
import time

import pytest

import seqtree

# spawn, not fork: a forked child inherits the parent's already-imported extension module and its
# open handles, which is not how a Snakemake worker or an xdist process actually starts.
CTX = mp.get_context("spawn")

AA = "ACDEFGHIKLMNPQRSTVWY"

#: The payload has to be big enough to keep the write window open. A few hundred short sequences
#: serialize in well under a millisecond, and a reader will essentially never land inside that --
#: an earlier version of this file used 400 refs and passed happily against the *broken* code. At
#: ~180k refs the index is tens of MB and takes ~50 ms to write, which is the real window a
#: cold-cache fan-out races against.
N_REFS = 180_000


def _refs():
    out = []
    for i in range(N_REFS):
        x, y, z, w = i % 20, (i // 20) % 20, (i // 400) % 20, (i // 8000) % 20
        out.append(f"CASS{AA[x]}{AA[y]}GQ{AA[z]}YEQ{AA[w]}F{AA[(i // 160000) % 20]}")
    return list(dict.fromkeys(out))


def _save(path):
    """Serialize a fresh index into `path` (the writer)."""
    idx = seqtree.Index.build(_refs(), "aa")
    idx.save(path)
    return "wrote"


def _poll_then_load(path):
    """Poll for the file the way a fan-out worker does, then load the instant it appears."""
    for _ in range(4000):
        if os.path.exists(path):
            try:
                return ("ok", len(seqtree.Index.load(path)))
            except Exception as e:  # noqa: BLE001 -- any failure at all is the bug
                return ("torn", f"{type(e).__name__}: {e}")
        time.sleep(0.0005)
    return ("timeout", None)


def _load_control(args):
    cache_dir, stagger = args
    time.sleep(stagger)
    try:
        return ("ok", len(seqtree.load_control("human_trb_aa", cache_dir=cache_dir)))
    except Exception as e:  # noqa: BLE001
        return ("fail", f"{type(e).__name__}: {e}")


def test_a_reader_never_observes_a_half_written_index(tmp_path):
    """The regression. A reader racing a writer sees the complete file or no file, never a stub."""
    path = str(tmp_path / "idx.sqtree")
    for _ in range(5):
        if os.path.exists(path):
            os.remove(path)
        with CTX.Pool(2) as pool:
            reader = pool.apply_async(_poll_then_load, (path,))
            writer = pool.apply_async(_save, (path,))
            writer.get(timeout=60)
            status, payload = reader.get(timeout=60)
        assert status != "torn", f"reader observed a partially written index: {payload}"
        if status == "ok":
            assert payload == len(_refs())


def test_save_leaves_no_temporary_behind(tmp_path):
    path = str(tmp_path / "idx.sqtree")
    small = ["CASSLGQAYEQYF", "CASSPGQAYEQYF"]
    seqtree.Index.build(small, "aa").save(path)
    seqtree.Index.build(small, "aa").save(path)  # overwrite an existing index
    assert len(seqtree.Index.load(path)) == 2
    assert [f.name for f in tmp_path.iterdir()] == ["idx.sqtree"]


def test_a_failed_save_leaves_no_temporary_behind(tmp_path):
    """A save into a non-existent directory must clean up its temporary, not litter."""
    bad = str(tmp_path / "no-such-dir" / "idx.sqtree")
    with pytest.raises(Exception):  # noqa: B017 -- RuntimeError from the C++ layer
        seqtree.Index.build(["CASSLGQAYEQYF"], "aa").save(bad)
    assert list(tmp_path.iterdir()) == []


def test_concurrent_writers_leave_one_valid_index(tmp_path):
    """Last writer wins, and what lands on disk is always a complete, loadable index."""
    path = str(tmp_path / "idx.sqtree")
    with CTX.Pool(6) as pool:
        assert pool.map(_save, [path] * 6) == ["wrote"] * 6
    assert len(seqtree.Index.load(path)) == len(_refs())
    assert [f.name for f in tmp_path.iterdir()] == ["idx.sqtree"]


def test_cold_cache_fanout_through_load_control(tmp_path):
    """The end-to-end shape: staggered workers, one empty cache dir, the real 250k control."""
    n = 8
    args = [(str(tmp_path), i * 0.05) for i in range(n)]
    with CTX.Pool(n) as pool:
        results = pool.map(_load_control, args)

    failures = [r for r in results if r[0] != "ok"]
    assert not failures, f"cold-cache fan-out failed: {failures}"
    assert {r[1] for r in results} == {250_000}
    assert not [f for f in os.listdir(tmp_path) if ".tmp" in f], "stray temporary left behind"


def test_a_corrupt_cache_is_rebuilt_rather_than_raised(tmp_path):
    """A cache truncated by a full disk, or written by an older seqtree, must self-heal."""
    cache = tmp_path / "control_human_trb_aa_bundled.sqtree"
    cache.write_bytes(b"SQTR" + b"\x00" * 32)  # right magic, garbage body
    idx = seqtree.load_control("human_trb_aa", cache_dir=str(tmp_path))
    assert len(idx) == 250_000
    assert len(seqtree.Index.load(str(cache))) == 250_000  # and the cache was replaced


def test_the_lock_is_optional_and_correctness_does_not_depend_on_it(monkeypatch):
    """A bare ``pip install seqtree`` has no filelock -- it arrives only with huggingface_hub.

    The lock saves redundant work; the atomic save is what makes concurrent use *correct*. If this
    ever regresses into a hard import, a default install starts raising ImportError on first use.
    """
    import builtins

    from seqtree import control

    assert control._build_lock("/tmp/x").__class__.__name__ in {"FileLock", "UnixFileLock",
                                                                "WindowsFileLock", "_NoLock"}

    real_import = builtins.__import__

    def no_filelock(name, *args, **kwargs):
        if name == "filelock":
            raise ImportError("simulated: not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", no_filelock)
    lock = control._build_lock("/tmp/x")
    assert lock.__class__.__name__ == "_NoLock"
    with lock:  # must be a usable context manager, not a stub that explodes
        pass
