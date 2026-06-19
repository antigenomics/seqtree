import seqtree
from seqtree import layout, pmhc

# The Dolton et al. (Cell 2023) HLA-A*02:01 cross-reactive trio (shared central motif
# x-x-x-A/G-I/L-G-I-x-x-x). The biological positive control for TCR-facing homology.
TRIPLE = ["EAAGIGILTV", "LLLGIGILVL", "NLSALGIFST"]


def _triple_store(k=4):
    recs = [{"epitope": e, "mhc": "HLA-A*02:01", "mhc_class": "MHCI", "gene": g}
            for e, g in zip(TRIPLE, ["MLANA", "BST2", "IGF2BP2"])]
    return pmhc.PMHCStore.from_records(recs, k=k)


def test_anchor_masking_and_kmers():
    spec = layout.DEFAULTS["mhc1"]  # P2 + C-terminus
    assert spec.resolve(10) == {1, 9}
    masked = layout.mask_anchors("EAAGIGILTV", spec)
    assert masked[1] == layout.MASK and masked[9] == layout.MASK
    assert masked[0] == "E" and masked[3] == "G"
    km = layout.kmers("EAAGIGILTV", 4, spec)
    assert "GIGI" in km  # central TCR-facing motif survives masking


def test_cell_triple_mutual_homologs():
    store = _triple_store()
    for q in TRIPLE:
        hits = store.search_homologs(q, "mhc1", mhc="HLA-A*02:01", max_subs=2, min_shared=1)
        found = {h.epitope for h in hits}
        assert found == set(TRIPLE) - {q}, f"{q}: {found}"


def test_mhc_filter():
    recs = [{"epitope": "EAAGIGILTV", "mhc": "HLA-A*02:01", "mhc_class": "MHCI"},
            {"epitope": "LLLGIGILVL", "mhc": "HLA-B*07:02", "mhc_class": "MHCI"}]
    store = pmhc.PMHCStore.from_records(recs, k=4)
    # same-allele restriction yields no A*02:01 homolog (the only homolog is B*07:02)
    assert store.search_homologs("EAAGIGILTV", "mhc1", mhc="HLA-A*02:01", max_subs=2) == []
    # querying B*07:02 finds the B*07:02 peptide
    assert any(h.epitope == "LLLGIGILVL"
               for h in store.search_homologs("EAAGIGILTV", "mhc1", mhc="HLA-B*07:02", max_subs=2))


def test_class_ii_register_agnostic_local():
    # class II: a shared 9-mer core inside longer peptides with different flanks (trim/shift)
    recs = [{"epitope": "AAA" + "PKYVKQNTLKLAT" + "AAA", "mhc": "HLA-DRB1*01:01", "mhc_class": "MHCII"},
            {"epitope": "WWWW" + "PKYVKQNTLKLAT", "mhc": "HLA-DRB1*01:01", "mhc_class": "MHCII"}]
    store = pmhc.PMHCStore.from_records(recs, k=5)
    hits = store.search_homologs("GG" + "PKYVKQNTLKLAT" + "GG", "mhc2", max_subs=0, min_shared=2)
    assert len(hits) == 2  # both found via the shared core regardless of flanks


def test_find_mimics_evalue():
    neo = "EAAGIGILTV"
    self_set = ["LLLGIGILVL"] + ["".join(c) for c in
                [("CASSLEAFF"), ("KLGGALQAK"), ("GILGFVFTL"), ("SLYNTVATL")]]
    bact = {"ecoli": ["NLSALGIFST", "MMMMMMMMM", "PQRSTVWYA"]}
    res = pmhc.find_mimics(neo, self_set, bacterial_sets=bact, max_subs=2, min_shared=1)
    assert "self" in res and "ecoli" in res
    self_hits = {h.epitope for h in res["self"]["hits"]}
    assert "LLLGIGILVL" in self_hits
    assert any(h.epitope == "NLSALGIFST" for h in res["ecoli"]["hits"])
    for r in res.values():
        assert 0.0 <= r["p_enrichment"] <= 1.0 and r["E"] >= 0.0


def test_assign_allele():
    # peptides with A*02:01-like anchors (P2=L, PΩ=V) vs B*07-like (P2=P, PΩ=L)
    recs = ([{"epitope": "ALAAAAAAV", "mhc": "HLA-A*02:01", "mhc_class": "MHCI"} for _ in range(5)]
            + [{"epitope": "GLDDDDDDV", "mhc": "HLA-A*02:01", "mhc_class": "MHCI"} for _ in range(5)]
            + [{"epitope": "APRRRRRRL", "mhc": "HLA-B*07:02", "mhc_class": "MHCI"} for _ in range(5)])
    store = pmhc.PMHCStore.from_records(recs, k=4)
    ranked = store.assign_allele("KLEEEEEEV", "mhc1")  # P2=L, PΩ=V -> A*02:01-like
    assert ranked[0][0] == "HLA-A*02:01"
