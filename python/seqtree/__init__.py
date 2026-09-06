"""seqtree: fast fuzzy search over biological sequences (amino acid / nucleotide).

Build an immutable index once, then search single queries or massive batches in
parallel. Two engines: ``seqtm`` (branch-and-bound, exact per-type edit caps, fast
Hamming path) and ``seqtrie`` (full-width DP carried down the trie, budget-only --
it ignores the per-type caps). Payload-agnostic: results are
``(ref_id, score, n_subs, n_ins, n_dels)``.

For anchored loops (CDR3 / junction), :mod:`seqtree.gapblock` restricts the alignment
to one contiguous indel and picks its position with a gap prior.

For the other shape of the problem -- a short query against a long *text* (a proteome, a
genome) rather than against a set of equal-length references -- :class:`TextIndex` does
exact k-mismatch search from a single index that serves every query length.
"""
from importlib.metadata import version as _distribution_version

from ._core import (
    Index,
    SearchParams,
    Hit,
    Alignment,
    SubstitutionMatrix,
    PositionalMatrix,
    KmerIndex,
    Candidate,
    TextIndex,
    TextResult,
    TextHit,
    ArrayView,
    pairwise_batch,
    alphabet_symbols,
    amino_acids,
)
from .control import load_control
from .evalue import evalues, thetas_from_scores, threshold_for_evalue
from . import distance, gapblock, layout, pairwise, pmhc, seeds
from .gapblock import (
    GapBlockIndex, IslandProfile, ScoreMatrix, central_prior, embed_in_frame, frame_prior,
    gapblock_score, positions_prior, profile_prior, score_matrix,
)
from .seeds import SeedIndex, core_kmers
from .pmhc import PMHCStore, find_mimics

#: Read from the installed distribution metadata, so ``pyproject.toml``'s ``project.version``
#: is the single source. A hand-copied literal here would silently drift from the version
#: actually published to PyPI.
__version__ = _distribution_version("seqtree")

__all__ = [
    "distance",
    "gapblock",
    "layout",
    "pairwise",
    "pmhc",
    "seeds",
    "GapBlockIndex",
    "IslandProfile",
    "ScoreMatrix",
    "gapblock_score",
    "score_matrix",
    "central_prior",
    "profile_prior",
    "frame_prior",
    "positions_prior",
    "embed_in_frame",
    "SeedIndex",
    "core_kmers",
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
    "TextIndex",
    "TextResult",
    "TextHit",
    "ArrayView",
    "pairwise_batch",
    "alphabet_symbols",
    "amino_acids",
    "load_control",
    "evalues",
    "thetas_from_scores",
    "threshold_for_evalue",
]
