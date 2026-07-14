"""Needleman-Wunsch and Smith-Waterman: ordinary protein alignment, without BioPython.

Everything else in seqtree *minimises a non-negative penalty* -- that is what a search ball and
an E-value need. This module does the opposite: it **maximises a raw log-odds similarity**, the
way BLAST and BioPython do, because that is what an ordinary pairwise alignment means and what
downstream code expects to get back.

The two views live on the same :class:`~seqtree.SubstitutionMatrix`::

    mat.penalty(a, b)     >= 0, zero on the diagonal   -- search, E-values, gap blocks
    mat.similarity(a, b)  signed log-odds              -- the aligners here

The penalty is the Gram transform of the similarity, ``pen = s(a,a) + s(b,b) - 2·s(a,b)``, which
is **lossy**: it forces the diagonal to zero and destroys ``s(a,a)``. So the raw grid is kept
rather than reconstructed, and a similarity score cannot be recovered from a penalty.

**Conventions, all verified against BioPython** (``tests/python/test_pairwise.py`` runs 6,720
comparisons across three matrices, ten gap settings and both modes; zero disagreements):

* a gap run of length ``L`` costs ``gap_open + (L-1)·gap_extend`` -- ``gap_open`` is the cost of
  the *first* gap column, not a surcharge on top of it;
* ``gap_open == gap_extend`` gives **linear** gaps. There is no separate mode for it;
* ``mode="global"`` charges end gaps like any other (true Needleman-Wunsch, not semi-global);
* ``mode="local"`` never lets the score fall below zero and takes the best cell anywhere
  (Smith-Waterman).

Gap costs are **positive magnitudes** and are subtracted. BLAST's protein defaults are
``gap_open=11, gap_extend=1``; BioPython's ``PairwiseAligner("blastp")`` preset uses ``12, 1``.

Example:
    >>> import seqtree
    >>> from seqtree.pairwise import align, score
    >>> mat = seqtree.SubstitutionMatrix.blosum62()
    >>> score("CASSLGQAYEQYF", "CASSPGQAYEQF", mat)          # global, BLAST defaults
    45
    >>> score("CASSLGQAYEQYF", "CASSPGQAYEQF", mat, gap_open=12)   # BioPython's 'blastp' preset
    44
    >>> aln = align("WWWAAAWWW", "KKKAAAKKK", mat, mode="local")   # Smith-Waterman
    >>> aln.score, aln.aligned_query, aln.aligned_ref
    (12, 'AAA', 'AAA')
"""
from __future__ import annotations

from collections.abc import Sequence

from ._core import Alignment, ScoreMatrix, SubstitutionMatrix
from ._core import align_dist_matrix as _dist_matrix
from ._core import align_pair as _align_pair
from ._core import align_score as _align_score
from ._core import align_score_matrix as _score_matrix

__all__ = ["score", "align", "score_matrix", "dist_matrix"]


def score(
    query: str,
    ref: str,
    matrix: SubstitutionMatrix,
    mode: str = "global",
    gap_open: int = 11,
    gap_extend: int = 1,
    alphabet: str = "aa",
) -> int:
    """Optimal alignment score of ``query`` against ``ref``.

    Args:
        query: First sequence.
        ref: Second sequence.
        matrix: Scoring matrix; its ``similarity`` view is used, not its penalty.
        mode: ``"global"`` for Needleman-Wunsch, ``"local"`` for Smith-Waterman. ``"nw"`` and
            ``"sw"`` are accepted too.
        gap_open: Cost of the first column of a gap. Positive; it is subtracted.
        gap_extend: Cost of each further column. Equal to ``gap_open`` means linear gaps.
        alphabet: ``"aa"``, ``"nt"`` or ``"iupac"``.

    Returns:
        The score, signed. Higher is more similar -- the opposite sense to the rest of seqtree.

    Raises:
        ValueError: On a negative gap cost, an unknown mode, or a symbol outside the alphabet.

    Example:
        >>> m = SubstitutionMatrix.blosum62()
        >>> score("AAA", "AAA", m)
        12
        >>> score("AAA", "AAAAA", m)          # a length-2 gap: 11 + 1*1 = 12
        0
    """
    return _align_score(query, ref, matrix, mode=mode, gap_open=gap_open,
                        gap_extend=gap_extend, alphabet=alphabet)


def align(
    query: str,
    ref: str,
    matrix: SubstitutionMatrix,
    mode: str = "global",
    gap_open: int = 11,
    gap_extend: int = 1,
    alphabet: str = "aa",
) -> Alignment:
    """As :func:`score`, but also returns the aligned strings and the edit ops.

    Returns:
        An ``Alignment``. Note ``Alignment.score`` here is a **similarity** (signed, higher is
        better), whereas the same field from :meth:`seqtree.Index.align` is a penalty. In local
        mode the aligned strings are the matched sub-sequences only.

    Example:
        >>> m = SubstitutionMatrix.blosum62()
        >>> a = align("CASSLGQAYEQYF", "CASSPGQAYEQF", m)
        >>> a.aligned_query, a.aligned_ref
        ('CASSLGQAYEQYF', 'CASSPGQAYEQ-F')
    """
    return _align_pair(query, ref, matrix, mode=mode, gap_open=gap_open,
                       gap_extend=gap_extend, alphabet=alphabet)


def score_matrix(
    queries: Sequence[str],
    refs: Sequence[str],
    matrix: SubstitutionMatrix,
    mode: str = "global",
    gap_open: int = 11,
    gap_extend: int = 1,
    alphabet: str = "aa",
    threads: int = 0,
) -> ScoreMatrix:
    """Every query against every reference, in C++ with the GIL released.

    Returns:
        A :class:`~seqtree.ScoreMatrix` of shape ``(len(queries), len(refs))`` holding signed
        similarity scores. ``numpy.asarray`` wraps it without copying.
    """
    return _score_matrix(list(queries), list(refs), matrix, mode=mode, gap_open=gap_open,
                         gap_extend=gap_extend, alphabet=alphabet, threads=threads)


def dist_matrix(
    queries: Sequence[str],
    refs: Sequence[str],
    matrix: SubstitutionMatrix,
    mode: str = "global",
    gap_open: int = 11,
    gap_extend: int = 1,
    alphabet: str = "aa",
    threads: int = 0,
) -> ScoreMatrix:
    """Alignment **distances**: ``d(a, b) = s(a,a) + s(b,b) - 2·s(a,b)``.

    The Gram transform, applied at the *sequence* level to the alignment scores rather than
    per residue. Non-negative, symmetric, zero on the diagonal -- so it is a distance, and it is
    what a prototype-distance embedding actually wants. This is the quantity users of BioPython
    hand-roll, and it is computed here without a Python loop: the self-scores are taken once per
    sequence, not once per pair.

    Returns:
        A :class:`~seqtree.ScoreMatrix` of shape ``(len(queries), len(refs))``.

    Example:
        >>> m = SubstitutionMatrix.blosum62()
        >>> d = dist_matrix(["CASSLGQAYEQYF"], ["CASSLGQAYEQYF", "CASSPGQAYEQF"], m)
        >>> d[0, 0], d[0, 1] > 0
        (0, True)
    """
    return _dist_matrix(list(queries), list(refs), matrix, mode=mode, gap_open=gap_open,
                        gap_extend=gap_extend, alphabet=alphabet, threads=threads)
