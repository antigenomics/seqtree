"""E-values for shared k-mer seeds in the variable core of an anchored loop.

The germline flanks of a CDR3 junction carry almost no evidence: an exact N-terminal
4-mer is shared by 31.0% of the 250k human TRB control repertoire (``CASS`` alone by 56.5%),
and a C-terminal 4-mer by 14.1%. A *central* 4-mer is shared by 0.080% -- about 386x more
selective. Same four residues, ~2e4-fold different evidence.

So the significance of a shared seed has to be **computed, not assumed**. For a seed ``w``
and a target set of ``N`` sequences, the expected number of chance targets sharing ``w``
with the query is

.. math::  E_{\\mathrm{seed}}(w) = N \\cdot n_C(w) / M

where ``n_C(w)`` counts *control sequences containing* ``w`` and ``M = |C|``. Under the null
the target draw is independent of the query, so conditioning on "the query contains ``w``"
is vacuous -- do not square the probability. Occurrence-weighted over the bundled control
with ``N = 1e5``, the median ``E_seed`` of a central k-mer is 20.8 (k=4), 2.0 (k=5), 0.40
(k=6). A *typical* shared central 4-mer is therefore not significant, but 4.9% of them are;
the median crosses ``E_seed < 1`` at k=6.

This cannot be modelled. The residual KL divergence of the empirical central-k-mer
distribution from a fitted background is 0.85 / 2.28 / 5.46 bits (independent per-position),
0.49 / 1.77 / 4.79 (Markov-1) and 0.43 / 1.60 / 4.48 (Markov-2) at k = 4 / 5 / 6 -- growing
with k, because D-gene germline runs (``GGG``, ``LAGG``, ``SGGG``) correlate. Count directly.

**Scope.** Seeds buy *precision*, not recall. Among same-epitope VDJdb pairs that lie in
different sequence islands (i.e. beyond the reach of any anchored alignment), only 0.5%
share a central 4-mer at all and 0.0% share a 6-mer -- a real 4x enrichment over
cross-epitope pairs, but negligible coverage. Use this to decide whether a seed you already
found is meaningful, and to prune uninformative seeds before gathering. Do not expect it to
connect distant relatives.
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence

from ._core import Index

__all__ = ["core_kmers", "SeedIndex"]


def core_kmers(seq: str, k: int, flank: int = 4) -> set[str]:
    """The distinct ``k``-mers of ``seq``'s variable core, excluding ``flank`` residues at
    each end.

    ``flank=4`` is the junction default: it drops the conserved Cys plus the first three
    germline residues at the 5' end, and the conserved Phe/Trp plus three at the 3' end.
    Sequences with a core shorter than ``k`` yield nothing.

    Example:
        >>> sorted(core_kmers("CASSLGQAYEQYF", 4))   # core is 'LGQAY'
        ['GQAY', 'LGQA']
    """
    if k < 1:
        raise ValueError("k must be >= 1")
    if flank < 0:
        raise ValueError("flank must be >= 0")
    core = seq[flank:len(seq) - flank] if flank else seq
    if len(core) < k:
        return set()
    return {core[i:i + k] for i in range(len(core) - k + 1)}


class SeedIndex:
    """Inverted index over the core k-mers of a sequence set.

    Built over a *control* repertoire it calibrates seed significance
    (:meth:`evalue`, :meth:`significant`). Built over a *target* set it gathers candidates
    (:meth:`gather`). Counts are per sequence, not per occurrence: the event is "sequence
    ``x`` contains ``w``".
    """

    def __init__(self, seqs: Iterable[str], k: int = 5, flank: int = 4):
        self.k = k
        self.flank = flank
        self.postings: dict[str, list[int]] = {}
        n = 0
        for i, s in enumerate(seqs):
            n += 1
            for w in core_kmers(s, k, flank):
                self.postings.setdefault(w, []).append(i)
        self.n = n

    @classmethod
    def from_index(cls, index: Index, k: int = 5, flank: int = 4) -> "SeedIndex":
        """Build from an existing :class:`Index`, e.g. the one from :func:`load_control`."""
        return cls((index.ref_seq(i) for i in range(len(index))), k, flank)

    def __len__(self) -> int:
        return self.n

    def count(self, seed: str) -> int:
        """Number of indexed sequences whose core contains ``seed``."""
        return len(self.postings.get(seed, ()))

    def evalue(self, seed: str, n_target: int) -> float:
        """``E_seed = n_target * count(seed) / len(self)``.

        The expected number of sequences in a target set of size ``n_target`` that share
        ``seed`` with the query by chance. Below 1, the shared seed is itself evidence.
        Uses the rule of three (``3/M``) when the seed is absent from the control, matching
        :mod:`seqtree.evalue`.
        """
        if self.n == 0:
            raise ValueError("seed index is empty")
        c = self.count(seed)
        return n_target * (3.0 if c == 0 else float(c)) / self.n

    def seed_evalues(self, query: str, n_target: int) -> dict[str, float]:
        """``E_seed`` for every core k-mer of ``query``, keyed by seed."""
        return {w: self.evalue(w, n_target) for w in core_kmers(query, self.k, self.flank)}

    def significant(self, query: str, n_target: int, alpha: float = 1.0) -> list[str]:
        """Core k-mers of ``query`` with ``E_seed < alpha``, rarest first.

        ``alpha=1.0`` keeps seeds expected to arise fewer than once by chance.
        """
        ev = self.seed_evalues(query, n_target)
        return sorted((w for w, e in ev.items() if e < alpha), key=lambda w: ev[w])

    def union_evalue(self, query: str, n_target: int, seeds: Sequence[str] | None = None) -> float:
        """E-value for sharing *at least one* of ``query``'s seeds with a target sequence.

        The query carries several overlapping core k-mers, so the per-seed E-values form a
        family. Rather than correct for that multiplicity, count the union directly: it is a
        single measurable set, so ``E = n_target * |union of postings| / M`` is exact and
        needs no correction constant. (A Boole bound over the seeds is only ~7% loose here,
        because overlapping string seeds gather nearly disjoint background sets -- but it is
        unnecessary.)
        """
        if self.n == 0:
            raise ValueError("seed index is empty")
        ws = core_kmers(query, self.k, self.flank) if seeds is None else set(seeds)
        hit: set[int] = set()
        for w in ws:
            hit.update(self.postings.get(w, ()))
        return n_target * (3.0 if not hit else float(len(hit))) / self.n

    def gather(self, query: str, seeds: Sequence[str] | None = None) -> set[int]:
        """Indices of sequences sharing at least one of ``query``'s core k-mers.

        Pass ``seeds=control.significant(query, ...)`` to gather on informative seeds only:
        an N-terminal 4-mer would otherwise pull in nearly half the repertoire.
        """
        ws = core_kmers(query, self.k, self.flank) if seeds is None else set(seeds)
        hit: set[int] = set()
        for w in ws:
            hit.update(self.postings.get(w, ()))
        return hit
