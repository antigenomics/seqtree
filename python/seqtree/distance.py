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

The same module also *enumerates* a Hamming ball rather than scoring a pair you already hold:

* :func:`neighbourhood` -- every sequence within ``r`` substitutions of one centre;
* :func:`neighbourhood_union` -- the same over **many** centres, each distinct sequence emitted
  once. Deduplication happens during generation, never as a pass over the ``sum(19*L_i)``
  multiset;
* :func:`union_size` -- the cardinality alone, for sizing a job before running it.

Comparison is **case-sensitive**, byte for byte -- unlike the search engines, which fold case.
For a *weighted* alignment (a substitution matrix, affine gaps, local mode), use
:mod:`seqtree.pairwise` instead.

Example:
    >>> from seqtree.distance import hamming, levenshtein, union_size
    >>> hamming("CASSLGQYF", "CASSPGQYF")
    1
    >>> levenshtein("kitten", "sitting")
    3
    >>> union_size(["CASSLGQYF", "CASSPGQYF"])      # 2 * 172 balls, 20 sequences shared
    324
"""
from __future__ import annotations

from collections.abc import Iterable, Sequence

from ._core import ScoreMatrix
from ._core import amino_acids as _amino_acids
from ._core import hamming as _hamming
from ._core import hamming_matrix as _hamming_matrix
from ._core import levenshtein as _levenshtein
from ._core import levenshtein_matrix as _levenshtein_matrix

__all__ = ["hamming", "levenshtein", "hamming_matrix", "levenshtein_matrix",
           "neighbourhood", "neighbourhood_union", "union_size"]

#: The 20 standard residues: :func:`~seqtree.amino_acids` minus the ambiguity codes ``B``/``Z``/
#: ``X`` and the stop ``*``, none of which is a substitution anybody means.
_STANDARD_AA = "".join(c for c in _amino_acids() if c not in "BZX*")


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


def _shells(seqs: Iterable[str], r: int, alphabet: str | None,
            include_self: bool, shell: bool) -> list[set[str]]:
    """The selected shells of the union, as sets: shell ``r`` alone, or shells ``0..r``.

    Multi-source breadth-first walk with one shared ``seen`` set, so a sequence lands in the
    shell of its *nearest* centre and is generated at most once -- the ``sum(19*L_i)`` multiset
    is never materialised. Substitution only, so centres of different lengths never meet.
    """
    if r < 0:
        raise ValueError(f"radius must be non-negative, got r={r}")
    alpha = _STANDARD_AA if alphabet is None else alphabet
    frontier = set(seqs)
    seen = set(frontier)
    levels = [frontier]
    for _ in range(r):
        nxt: set[str] = set()
        for s in frontier:
            for i, c in enumerate(s):
                head, tail = s[:i], s[i + 1:]
                for a in alpha:
                    if a != c:
                        t = head + a + tail
                        if t not in seen:
                            seen.add(t)
                            nxt.add(t)
        levels.append(nxt)
        frontier = nxt
    if not include_self:
        levels[0] = set()       # the centres are exactly the r = 0 shell
    return [levels[r]] if shell else levels


def neighbourhood_union(seqs: Iterable[str], r: int = 1, alphabet: str | None = None,
                        include_self: bool = True, shell: bool = False) -> list[str]:
    """Every sequence within ``r`` substitutions of **any** of ``seqs``, each listed once.

    This is the union of the per-sequence Hamming balls, not their concatenation. Near-duplicate
    centres -- co-specific TCR junctions, an error family around one UMI -- have balls that
    overlap heavily, so the difference is the dominant term rather than a correction. The
    dedup happens as the ball is walked, so the ``sum(19*L_i)`` multiset never exists.

    Args:
        seqs: The centres. Duplicates and mixed lengths are fine; substitution-only means
            balls of different lengths cannot overlap.
        r: Radius in substitutions.
        alphabet: Substituting alphabet; ``None`` is the 20 standard amino acids. Characters of
            ``seqs`` outside it are still substitutable, so a centre containing ``X`` has
            ``|alphabet|`` neighbours at that position rather than ``|alphabet| - 1``.
        include_self: Keep the centres themselves. ``False`` drops **every** centre, including
            one that happens to be a neighbour of another -- the ``r = 0`` shell of a union is
            exactly ``set(seqs)``.
        shell: Return only the members at distance exactly ``r`` from the *nearest* centre,
            instead of the whole closed ball.

    Returns:
        The distinct members, ordered by distance and sorted within each shell (so the centres
        come first). Empty input gives ``[]``.

    Raises:
        ValueError: If ``r`` is negative.

    Example:
        >>> neighbourhood_union(["AA", "AC"], 1, alphabet="AC")
        ['AA', 'AC', 'CA', 'CC']
        >>> a = neighbourhood("CASSLGQYF")
        >>> b = neighbourhood("CASSPGQYF")                  # one substitution away
        >>> len(a) + len(b), len(neighbourhood_union(["CASSLGQYF", "CASSPGQYF"]))
        (344, 324)
    """
    return [s for level in _shells(seqs, r, alphabet, include_self, shell) for s in sorted(level)]


def neighbourhood(seq: str, r: int = 1, alphabet: str | None = None,
                  include_self: bool = True, shell: bool = False) -> list[str]:
    """Every sequence within ``r`` substitutions of ``seq``: its closed Hamming ball.

    Substitution only, so every member has ``len(seq)`` characters -- Hamming distance is
    undefined across lengths (see :func:`hamming`), and an indel is a different question. Over a
    ``k``-letter alphabet the ball holds ``sum((k-1)^d * C(L, d) for d in 0..r)`` sequences;
    at ``r = 1`` that is ``19*L + 1``.

    Args:
        seq: The centre.
        r: Radius in substitutions.
        alphabet: Substituting alphabet; ``None`` is the 20 standard amino acids.
        include_self: Keep ``seq`` itself (the ``d = 0`` member).
        shell: Return only the members at distance exactly ``r``, instead of the closed ball.
            The cognacy-retention profile is estimated per shell, not per ball.

    Returns:
        The distinct members, ordered by distance and sorted within each shell, so ``seq``
        itself comes first when ``include_self``.

    Raises:
        ValueError: If ``r`` is negative.

    Example:
        >>> len(neighbourhood("CASSLGQYF"))                 # 19 * 9 + 1
        172
        >>> len(neighbourhood("CASSLGQYF", include_self=False))
        171
        >>> neighbourhood("A", 1, alphabet="ACGT")
        ['A', 'C', 'G', 'T']
        >>> len(neighbourhood("CASSLGQYF", 2, shell=True))   # 19^2 * C(9, 2)
        12996
    """
    return neighbourhood_union([seq], r, alphabet, include_self, shell)


def union_size(seqs: Iterable[str], r: int = 1, alphabet: str | None = None,
               include_self: bool = True, shell: bool = False) -> int:
    """How many distinct sequences :func:`neighbourhood_union` would return.

    Same walk, but the result list is never built or sorted -- use it to size a job (an OLGA
    ``P_gen`` pass over the union, say) before committing to it. The dedup set is still held,
    so this bounds the output, not the peak memory.

    Args:
        seqs: The centres.
        r: Radius in substitutions.
        alphabet: Substituting alphabet; ``None`` is the 20 standard amino acids.
        include_self: Count the centres themselves.
        shell: Count only distance-exactly-``r`` members.

    Returns:
        ``len(neighbourhood_union(seqs, r, alphabet, include_self, shell))``.

    Raises:
        ValueError: If ``r`` is negative.

    Example:
        >>> union_size(["CASSLGQYF", "CASSPGQYF"])          # 344 with double-counting
        324
        >>> union_size(["CASSLGQYF", "CASSLGQYF"])          # a repeat adds nothing
        172
    """
    return sum(len(level) for level in _shells(seqs, r, alphabet, include_self, shell))
