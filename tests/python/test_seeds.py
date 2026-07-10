"""Seed E-values: the flank/core asymmetry, and the union rule for overlapping seeds."""
import gzip
import statistics as stt
from importlib import resources

import pytest

from seqtree.seeds import SeedIndex, core_kmers


def _control(limit=None):
    with resources.files("seqtree.data").joinpath("control_human_trb_aa.txt.gz").open("rb") as fh:
        seqs = [ln.strip() for ln in gzip.open(fh, "rt") if ln.strip()]
    seqs = [s for s in seqs if s and s[0].isalpha()]
    return seqs[:limit] if limit else seqs


def test_core_kmers_excludes_the_flanks():
    assert core_kmers("CASSLGQAYEQYF", 4) == {"LGQA", "GQAY"}      # core is 'LGQAY'
    assert core_kmers("CASSLGQAYEQYF", 4, flank=0) == {
        s for s in ("CASS", "ASSL", "SSLG", "SLGQ", "LGQA", "GQAY", "QAYE", "AYEQ", "YEQY", "EQYF")
    }
    assert core_kmers("CASSF", 4) == set()   # core shorter than k
    with pytest.raises(ValueError):
        core_kmers("CASSLGQAYEQYF", 0)


def test_counts_sequences_not_occurrences():
    # 'GGGG' appears at two core offsets in the first sequence but the event is "contains".
    si = SeedIndex(["CASSGGGGGEQYF", "CASSGGGGAEQYF", "CASSLLLLLEQYF"], k=4, flank=4)
    assert si.count("GGGG") == 2
    assert len(si) == 3


def test_evalue_is_n_times_the_background_rate():
    si = SeedIndex(["CASS" + core + "EQYF" for core in ("AAAAA", "AAAAA", "BBBBB", "CCCCC")],
                   k=5, flank=4)
    assert len(si) == 4
    assert si.count("AAAAA") == 2
    assert si.evalue("AAAAA", n_target=100) == pytest.approx(100 * 2 / 4)
    # rule of three for an unseen seed, matching seqtree.evalue
    assert si.evalue("WWWWW", n_target=100) == pytest.approx(100 * 3 / 4)


def test_significant_filters_and_orders_by_rarity():
    seqs = ["CASS" + c + "EQYF" for c in ["AAAAA"] * 50 + ["BBBBB"] * 5 + ["CCCCC"]]
    si = SeedIndex(seqs, k=5, flank=4)
    q = "CASSCCCCCEQYF"
    assert si.significant(q, n_target=10, alpha=1.0) == ["CCCCC"]   # E = 10*1/56 << 1
    assert si.significant("CASSAAAAAEQYF", n_target=10, alpha=1.0) == []  # E = 10*50/56 ~ 8.9


def test_union_evalue_counts_the_union_not_the_sum():
    # Two overlapping core 4-mers of the query both hit the SAME control sequence, so the
    # union count is 1, while summing the per-seed E-values would double count it.
    si = SeedIndex(["CASSLGQAYEQYF", "CASSWWWWWEQYF"], k=4, flank=4)
    q = "CASSLGQAYEQYF"
    assert core_kmers(q, 4) == {"LGQA", "GQAY"}
    assert si.gather(q) == {0}
    assert si.union_evalue(q, n_target=10) == pytest.approx(10 * 1 / 2)
    boole = sum(si.evalue(w, 10) for w in core_kmers(q, 4))
    assert boole > si.union_evalue(q, n_target=10)   # Boole is an upper bound, and loose here


def test_gather_on_significant_seeds_only():
    seqs = ["CASS" + c + "EQYF" for c in ["AAAAA"] * 50 + ["CCCCC"]]
    si = SeedIndex(seqs, k=5, flank=4)
    q = "CASSCCCCCEQYF"
    assert len(si.gather(q)) == 1
    # a query carrying the common core gathers nearly everything unless filtered
    common = "CASSAAAAAEQYF"
    assert len(si.gather(common)) == 50
    assert si.gather(common, seeds=si.significant(common, n_target=10, alpha=1.0)) == set()


def test_from_index_round_trips():
    import seqtree as st
    seqs = ["CASSLGQAYEQYF", "CASSIRSSYEQYF"]
    idx = st.Index.build(list(seqs), "aa")
    a, b = SeedIndex.from_index(idx, k=4), SeedIndex(seqs, k=4)
    assert a.postings == b.postings and len(a) == len(b)


# ---------------------------------------------------------------- the measured science

def test_flank_seeds_are_useless_and_core_seeds_are_not():
    """The whole reason this module exists: an N-terminal 4-mer selects ~half the repertoire,
    a central 4-mer ~0.1%."""
    ctrl = _control()
    assert len(ctrl) == 250_000

    nterm = SeedIndex(ctrl, k=4, flank=0)
    core = SeedIndex(ctrl, k=4, flank=4)
    M = len(ctrl)

    # 'CASS' is the single most common N-terminal 4-mer.
    assert nterm.count("CASS") / M > 0.6
    # ...whereas the most common central 4-mer is far rarer.
    top_core = max(core.postings, key=lambda w: len(core.postings[w]))
    assert core.count(top_core) / M < 0.05


@pytest.mark.parametrize("k,lo,hi", [(4, 10.0, 100.0), (6, 0.0, 1.0)])
def test_median_seed_evalue_crosses_one_at_k6(k, lo, hi):
    """Occurrence-weighted median E_seed of a central k-mer, N = 1e5: 32.4 at k=4, 0.80 at
    k=6. This is the quantitative form of 'the same 4 residues from the NDN region is
    sometimes hardly by chance'."""
    ctrl = _control()
    si = SeedIndex(ctrl, k=k, flank=4)
    N = 100_000
    # weight each seed by how often a repertoire sequence carries it
    evs = [si.evalue(w, N) for w, post in si.postings.items() for _ in range(len(post))]
    med = stt.median(evs)
    assert lo <= med < hi, f"k={k}: median E_seed = {med:.2f}, expected in [{lo}, {hi})"
