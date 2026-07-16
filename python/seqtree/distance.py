"""Plain string edit distances: Hamming and Levenshtein, in C++, without a dependency.

These are the *unweighted* distances -- unit costs, no substitution matrix, no gap model, no
alphabet. That is the whole point: when all you need is "how many edits apart are these two
strings", you should not have to build an :class:`~seqtree.SubstitutionMatrix` or reach for
``python-Levenshtein`` / ``rapidfuzz``. seqtree still needs nothing at runtime.

* :func:`hamming` -- number of differing positions; defined only for **equal-length** sequences,
  and raises :class:`ValueError` otherwise;
* :func:`levenshtein` -- the classic insertion / deletion / substitution edit distance, each edit
  costing 1;
* :func:`hamming_matrix` / :func:`levenshtein_matrix` -- every ``a`` against every ``b`` in one
  GIL-released, multi-threaded C++ call, returned as a zero-copy :class:`~seqtree.ScoreMatrix`.

Comparison is **case-sensitive**, byte for byte -- unlike the search engines, which fold case.
For a *weighted* alignment (a substitution matrix, affine gaps, local mode), use
:mod:`seqtree.pairwise` instead.

Example:
    >>> from seqtree.distance import hamming, levenshtein
    >>> hamming("CASSLGQYF", "CASSPGQYF")
    1
    >>> levenshtein("kitten", "sitting")
    3
"""
from __future__ import annotations

from collections.abc import Sequence

from ._core import ScoreMatrix
from ._core import hamming as _hamming
from ._core import hamming_matrix as _hamming_matrix
from ._core import levenshtein as _levenshtein
from ._core import levenshtein_matrix as _levenshtein_matrix

__all__ = ["hamming", "levenshtein", "hamming_matrix", "levenshtein_matrix"]


def hamming(a: str, b: str) -> int:
    """Number of positions at which ``a`` and ``b`` differ.

    Args:
        a: First sequence.
        b: Second sequence, of the **same length** as ``a``.

    Returns:
        The count of differing positions (0 when identical).

    Raises:
        ValueError: If ``a`` and ``b`` have different lengths -- Hamming distance is undefined
            for unequal lengths; use :func:`levenshtein` for that.

    Example:
        >>> hamming("AAAA", "AAAA")
        0
        >>> hamming("AAAA", "ATAT")
        2
    """
    return _hamming(a, b)


def levenshtein(a: str, b: str) -> int:
    """Edit distance: fewest single-character insert / delete / substitute steps from ``a`` to ``b``.

    Args:
        a: First sequence.
        b: Second sequence; may be any length.

    Returns:
        The edit distance (0 when identical, ``max(len(a), len(b))`` at most).

    Example:
        >>> levenshtein("flaw", "lawn")
        2
        >>> levenshtein("CASSLGQAYEQYF", "CASSPGQAYEQF")
        2
    """
    return _levenshtein(a, b)


def hamming_matrix(a: Sequence[str], b: Sequence[str], threads: int = 0) -> ScoreMatrix:
    """Hamming distance of every ``a`` against every ``b``, in parallel C++.

    Args:
        a: Query sequences (the rows).
        b: Reference sequences (the columns).
        threads: Worker threads; ``0`` uses all cores.

    Returns:
        A :class:`~seqtree.ScoreMatrix` of shape ``(len(a), len(b))`` of int32 distances.
        ``numpy.asarray`` wraps it without copying.

    Raises:
        ValueError: If any ``(a[i], b[k])`` pair has mismatched lengths.

    Example:
        >>> import numpy as np
        >>> d = np.asarray(hamming_matrix(["AAAA", "AAAT"], ["AAAA", "TTTT"]))
        >>> d.tolist()
        [[0, 4], [1, 3]]
    """
    return _hamming_matrix(list(a), list(b), threads=threads)


def levenshtein_matrix(a: Sequence[str], b: Sequence[str], threads: int = 0) -> ScoreMatrix:
    """Levenshtein distance of every ``a`` against every ``b``, in parallel C++.

    Unlike :func:`hamming_matrix`, sequences may differ in length freely.

    Args:
        a: Query sequences (the rows).
        b: Reference sequences (the columns).
        threads: Worker threads; ``0`` uses all cores.

    Returns:
        A :class:`~seqtree.ScoreMatrix` of shape ``(len(a), len(b))`` of int32 distances,
        zero-copy through ``numpy.asarray``.
    """
    return _levenshtein_matrix(list(a), list(b), threads=threads)
