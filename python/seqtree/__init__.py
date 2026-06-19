"""seqtree: fast fuzzy search over biological sequences (amino acid / nucleotide).

Build an immutable index once, then search single queries or massive batches in
parallel. Two engines: ``seqtm`` (branch-and-bound, exact per-type edit caps,
fast Hamming path) and ``seqtrie`` (banded DP, matrix-weighted score budgets).
Payload-agnostic: results are ``(ref_id, score, n_subs, n_ins, n_dels)``.
"""
from ._core import (
    Index,
    SearchParams,
    Hit,
    Alignment,
    SubstitutionMatrix,
    PositionalMatrix,
    KmerIndex,
    Candidate,
    pairwise_batch,
    alphabet_symbols,
    amino_acids,
)
from .control import load_control
from .evalue import evalues
from . import layout, pmhc
from .pmhc import PMHCStore, find_mimics

__all__ = [
    "layout",
    "pmhc",
    "PMHCStore",
    "find_mimics",
    "Index",
    "SearchParams",
    "Hit",
    "Alignment",
    "SubstitutionMatrix",
    "PositionalMatrix",
    "KmerIndex",
    "Candidate",
    "pairwise_batch",
    "alphabet_symbols",
    "amino_acids",
    "load_control",
    "evalues",
]
