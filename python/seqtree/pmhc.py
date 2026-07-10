"""Peptide-MHC epitope homology search and neoantigen molecular-mimicry discovery.

Homology is a shared *TCR-facing* k-mer (anchors masked); see :mod:`seqtree.layout`.
Built on the C++ :class:`KmerIndex` (seed-and-gather). MHC restriction is a payload
filter; significance (E-values) lives in :mod:`seqtree.pmhc_evalue`.
"""
import csv
import gzip
from collections import Counter, defaultdict
from dataclasses import dataclass

from . import layout
from ._core import KmerIndex, SearchParams

# vdjdb/IEDB-style class labels -> our class keys
_CLASS = {"MHCI": "mhc1", "I": "mhc1", "mhc1": "mhc1",
          "MHCII": "mhc2", "II": "mhc2", "mhc2": "mhc2"}


@dataclass
class EpitopeHit:
    epitope: str
    mhc: str
    cls: str
    shared_kmers: int
    score: int
    gene: str = ""
    species: str = ""

    def __iter__(self):
        return iter((self.epitope, self.mhc, self.shared_kmers, self.score))


class _ClassIndex:
    """One class's peptides + payload + KmerIndex."""

    def __init__(self, cls, k, spec):
        self.cls = cls
        self.k = k
        self.spec = spec
        self.epitopes = []
        self.mhc = []
        self.gene = []
        self.species = []
        self.allele_to_id = {}
        self.kmer = None  # KmerIndex

    def add(self, epitope, mhc, gene, species):
        self.epitopes.append(epitope)
        self.mhc.append(mhc)
        self.gene.append(gene)
        self.species.append(species)
        self.allele_to_id.setdefault(mhc, len(self.allele_to_id))

    def build(self):
        allele_ids = [self.allele_to_id[m] for m in self.mhc]
        kmer_lists = [layout.kmers(p, self.k, self.spec) for p in self.epitopes]
        self.kmer = KmerIndex.build(kmer_lists, alphabet="aa", allele_ids=allele_ids)


class PMHCStore:
    """Searchable epitope store, partitioned by MHC class."""

    def __init__(self, k=4, anchor_overrides=None):
        self.k = k
        ov = anchor_overrides or {}
        self._cls = {c: _ClassIndex(c, k, layout.spec_for(c, ov.get(c)))
                     for c in ("mhc1", "mhc2")}

    # -- construction ---------------------------------------------------------
    @classmethod
    def from_records(cls, records, k=4, anchor_overrides=None):
        """records: iterable of dicts with epitope, mhc (or mhc_a), mhc_class[, gene, species]."""
        store = cls(k=k, anchor_overrides=anchor_overrides)
        AA = set("ACDEFGHIKLMNPQRSTVWY")
        for r in records:
            cl = _CLASS.get(str(r.get("mhc_class", "")).strip())
            ep = str(r.get("epitope", "")).strip().upper()
            if cl is None or not ep or not all(c in AA for c in ep):
                continue
            mhc = str(r.get("mhc") or r.get("mhc_a") or "").strip()
            store._cls[cl].add(ep, mhc, str(r.get("gene", "")), str(r.get("species", "")))
        store._build()
        return store

    @classmethod
    def from_pmhc(cls, path, classes=("mhc1", "mhc2"), species=None, k=4, anchor_overrides=None):
        """Stream the isalgo/pmhc_data TSV(.gz)."""
        csv.field_size_limit(10**7)
        op = gzip.open if str(path).endswith(".gz") else open
        keep = set(classes)
        with op(path, "rt") as fh:
            rows = csv.DictReader(fh, delimiter="\t")
            recs = []
            for row in rows:
                cl = _CLASS.get(str(row.get("mhc_class", "")).strip())
                if cl is None or cl not in keep:
                    continue
                if species and row.get("species") and row["species"] != species:
                    continue
                recs.append(row)
        return cls.from_records(recs, k=k, anchor_overrides=anchor_overrides)

    def _build(self):
        for ci in self._cls.values():
            ci.build()

    def __len__(self):
        return sum(len(ci.epitopes) for ci in self._cls.values())

    def size(self, cls):
        """Number of indexed epitopes for a presentation class (``"mhc1"`` / ``"mhc2"``)."""
        return len(self._cls[cls].epitopes)

    # -- search ---------------------------------------------------------------
    def search_homologs(self, query, cls, mhc=None, max_subs=1, matrix="", min_shared=1,
                        exclude_self=True, threads=0):
        """TCR-facing homologs of `query` in class `cls`, optionally restricted to `mhc`."""
        ci = self._cls[cls]
        if ci.kmer is None or ci.kmer.num_peptides() == 0:
            return []
        allele_filter = -1
        if mhc is not None:
            if mhc not in ci.allele_to_id:
                return []
            allele_filter = ci.allele_to_id[mhc]
        qk = layout.kmers(query.strip().upper(), ci.k, ci.spec)
        if not qk:
            return []
        p = SearchParams(max_subs=max_subs, engine="seqtm", matrix=matrix)
        cands = ci.kmer.seed_and_gather([qk], p, min_shared, allele_filter, threads)[0]
        out = []
        for c in cands:
            pid = c.peptide_id
            if exclude_self and ci.epitopes[pid] == query.strip().upper():
                continue
            out.append(EpitopeHit(ci.epitopes[pid], ci.mhc[pid], cls, c.shared_kmers, c.best_score,
                                  ci.gene[pid], ci.species[pid]))
        return out

    # -- reverse problem: which allele(s) present this peptide ----------------
    def assign_allele(self, query, cls, top=5):
        """Rank alleles by how typical the query's ANCHOR signature is among each
        allele's presented peptides (a lightweight presentation prior, not NetMHCpan).
        Returns [(allele, score, n_match, n_allele)] sorted by score desc.

        ``score`` is the log-odds of presentation by the allele vs the marginal
        background, so it doubles as a non-binder filter: ``score <= 0`` means the
        anchors are not enriched for that allele (a non-binder for it), and if every
        allele scores ``<= 0`` the peptide binds nothing in the panel. Class II is
        promiscuous, so expect several alleles with ``score > 0`` (multi-label)."""
        ci = self._cls[cls]
        spec = layout.DEFAULTS[cls] if not spec_has_anchors(ci.spec) else ci.spec
        q = query.strip().upper()
        qsig = _anchor_sig(q, spec)
        per_allele = defaultdict(Counter)
        totals = Counter()
        bg = Counter()
        for ep, m in zip(ci.epitopes, ci.mhc):
            sig = _anchor_sig(ep, spec)
            per_allele[m][sig] += 1
            totals[m] += 1
            bg[sig] += 1
        n_total = sum(totals.values()) or 1
        bg_p = bg[qsig] / n_total
        ranked = []
        for allele, n_allele in totals.items():
            n_match = per_allele[allele][qsig]
            p_a = n_match / n_allele if n_allele else 0.0
            # log-odds of presentation by this allele vs the marginal background
            import math
            score = math.log((p_a + 1e-9) / (bg_p + 1e-9))
            ranked.append((allele, round(score, 3), n_match, n_allele))
        ranked.sort(key=lambda x: x[1], reverse=True)
        return ranked[:top]


def spec_has_anchors(spec):
    """Does this layout pin any anchor positions? Class II specs do not."""
    return bool(spec.anchors)


def _anchor_sig(pep, spec):
    idx = sorted(spec.resolve(len(pep)))
    return "".join(pep[i] for i in idx) if idx else ""


def build_kmer_index(peptides, cls="mhc1", k=4, anchor_overrides=None):
    """Build a bare KmerIndex over a peptide set (for find_mimics / custom corpora)."""
    spec = layout.spec_for(cls, (anchor_overrides or {}).get(cls))
    return KmerIndex.build([layout.kmers(p.strip().upper(), k, spec) for p in peptides],
                           alphabet="aa"), spec


def find_mimics(neoantigen, self_set, bacterial_sets=None, control=None, cls="mhc1",
                k=4, max_subs=1, matrix="", min_shared=1, anchor_overrides=None, threads=0):
    """Discover TCR-facing molecular mimics of a neoantigen, with presentation-aware E-values.

    All inputs are assumed presented by the same (compatible) MHC -- cross-allele groove
    similarity is out of scope. ``self_set`` is the host presented peptidome (same allele);
    ``bacterial_sets`` maps organism -> presented peptides; ``control`` is the per-allele
    presented background (defaults to ``self_set``). Returns one entry per source with the
    homolog hits, the expected count E, and the enrichment p-value.
    """
    from .pmhc_evalue import homolog_evalue

    neoantigen = neoantigen.strip().upper()
    spec = layout.spec_for(cls, (anchor_overrides or {}).get(cls))
    qk = layout.kmers(neoantigen, k, spec)
    params = SearchParams(max_subs=max_subs, engine="seqtm", matrix=matrix)

    def homologs(peptides):
        seqs = [p.strip().upper() for p in peptides]
        idx = KmerIndex.build([layout.kmers(p, k, spec) for p in seqs], alphabet="aa")
        cands = idx.seed_and_gather([qk], params, min_shared, -1, threads)[0]
        return seqs, cands

    ctrl_seqs, ctrl_cands = homologs(control if control is not None else self_set)
    n_control = sum(1 for c in ctrl_cands if ctrl_seqs[c.peptide_id] != neoantigen)
    M = len(ctrl_seqs)

    out = {}
    sources = {"self": self_set}
    sources.update(bacterial_sets or {})
    for name, peptides in sources.items():
        seqs, cands = homologs(peptides)
        hits = [EpitopeHit(seqs[c.peptide_id], "", cls, c.shared_kmers, c.best_score)
                for c in cands if seqs[c.peptide_id] != neoantigen]
        ev = homolog_evalue(len(hits), n_control, len(seqs), M)
        out[name] = {"hits": sorted(hits, key=lambda h: (-h.shared_kmers, h.score)), **ev}
    return out
