import seqtree


def test_bundled_control_loads_and_caches(tmp_path):
    idx = seqtree.load_control(size=5000, cache_dir=str(tmp_path))
    assert len(idx) == 5000
    # cache file written; second call reloads it
    assert any(p.suffix == ".sqtree" for p in tmp_path.iterdir())
    idx2 = seqtree.load_control(size=5000, cache_dir=str(tmp_path))
    assert len(idx2) == 5000


def test_control_members_are_findable(tmp_path):
    idx = seqtree.load_control(size=2000, cache_dir=str(tmp_path))
    p = seqtree.SearchParams(max_subs=1, engine="seqtm")
    member = idx.ref_seq(0)
    assert any(h.ref_id == 0 and h.score == 0 for h in idx.search(member, p))


def test_unknown_control_raises():
    import pytest

    with pytest.raises(ValueError):
        seqtree.load_control(name="klingon_trb_aa", size=10)


def test_sanitize_keeps_only_productive_clonotypes():
    """The control IS the E-value null and must share the target's background law (ass:match).

    VDJtools marks out-of-frame rearrangements with '_' and in-frame stops with '*'. Both are
    non-coding. '_' cannot be repaired at the amino-acid level: VDJtools collapses a RUN of
    untranslatable positions into one character, so the residue count is already gone.
    Out-of-frame junctions escape thymic selection and estimate Pgen, not P0.
    """
    from seqtree.control import sanitize

    # the three real markers from the mouse TRB table, plus a clean junction
    kept, dropped = sanitize(["CASSLYEQYF", "C*A_FF", "CASS*GTGGYEQYF", "CASISRTV_NTGQLYF"])
    assert kept == ["CASSLYEQYF"] and dropped == 3

    # ambiguous residues go too: pen(X, a) is barely half a mismatch and pen(X, X) is zero,
    # so an X column is a cheap wildcard that would inflate the control ball mass.
    kept, dropped = sanitize(["CASXF", "CASZF", "CASBF", "CASSF"])
    assert kept == ["CASSF"] and dropped == 3

    kept, dropped = sanitize(["ACGT", "ACG_"], "nt")     # non-aa alphabets keep their symbol set
    assert kept == ["ACGT"] and dropped == 1


def test_bundled_control_is_productive():
    """A regression guard on the shipped asset itself."""
    from seqtree.control import _read_bundled, sanitize

    seqs = _read_bundled("human_trb_aa")
    kept, dropped = sanitize(seqs)
    assert dropped == 0, f"{dropped} non-productive sequences in the bundled control"
    assert len(kept) == 250_000


def test_a_prefix_of_the_bundled_control_is_a_uniform_sample():
    """load_control(size=k) takes bundled[:k], so the asset's ORDER is part of its contract.

    It used to be abundance-sorted (a public-clone head). Sorting it alphabetically would have
    been just as bad. gen_control.py shuffles with a fixed seed, so any prefix is representative.
    """
    import statistics as stt

    from seqtree.control import _read_bundled

    seqs = _read_bundled("human_trb_aa")
    whole_cass = sum(s.startswith("CASS") for s in seqs) / len(seqs)
    whole_len = stt.mean(map(len, seqs))
    for k in (2_000, 20_000, 50_000):
        pre = seqs[:k]
        cass = sum(s.startswith("CASS") for s in pre) / k
        assert abs(cass - whole_cass) < 0.03, f"prefix[:{k}] CASS share {cass:.3f} vs {whole_cass:.3f}"
        assert abs(stt.mean(map(len, pre)) - whole_len) < 0.2


def test_download_filters_productive_and_samples_uniformly(tmp_path, monkeypatch):
    """Pins both control defects at once, with a fake table and no network.

    The upstream tables are sorted by clonotype abundance, so `head -n size` returns the most
    expanded public clones -- not a sample of P0. And they carry '_' (out of frame) and '*' (stop),
    which are non-coding and belong to Pgen, not P0.
    """
    import gzip
    import sys
    import types

    import seqtree.control as control

    aa = "ACDEFGHIKLMNPQRSTVWY"
    rows = ["count\tcdr3aa"]
    productive = [f"CASS{aa[i % 20]}{aa[(i // 20) % 20]}{aa[(i // 400) % 20]}EQYF" for i in range(1000)]
    productive = list(dict.fromkeys(productive))
    for n, s in enumerate(productive):                     # descending "abundance", as upstream
        rows.append(f"{len(productive) - n}\t{s}")
    rows += ["5\tCASS_EQYF", "4\tCASS*EQYF", "3\tCASSXEQYF"]   # out of frame, stop, ambiguous
    path = tmp_path / "fake.tsv.gz"
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(rows) + "\n")

    monkeypatch.setitem(control._HF, "fake", ("repo", "fake.tsv.gz", "cdr3aa"))
    monkeypatch.setitem(sys.modules, "huggingface_hub",
                        types.SimpleNamespace(hf_hub_download=lambda **kw: str(path)))

    everything = control._download("fake", None)
    assert everything == productive
    assert not any(c in s for s in everything for c in "_*X")

    k = 50
    sample = control._download("fake", k, seed=0)
    assert len(sample) == k and set(sample) <= set(everything)
    assert sample != everything[:k], "reservoir returned the abundance head"

    assert control._download("fake", k, seed=0) == sample     # reproducible from the seed
    assert control._download("fake", k, seed=1) != sample
