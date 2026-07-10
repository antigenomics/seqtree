"""Single-gap-block alignment for anchored loops (CDR3 / junction).

CDR3 length variation comes from V/J trimming plus N-addition -- **one** contiguous
indel event, not scattered indels. So we restrict the alignment to exactly one gap block
of length ``d = abs(len(q) - len(r))``, placed anywhere, and pick its position by score.

Two facts make this cheap and exact:

* The optimum over all block positions is a prefix/suffix sum, so :func:`gapblock_score`
  is ``O(min(m, n))`` rather than the ``O(m*n)`` of a full DP.
* Each contiguous-``d``-deletion variant of the query *is* one block position. Hamming-
  matching those variants against the ordinary trie therefore reproduces the same optimum
  with no new engine -- see :class:`GapBlockIndex`.

The score is a non-negative penalty with ``s(q, q) == 0``, so it defines a ball in the
sense of ``appendix/evalue.tex`` and flows through :func:`seqtree.evalues` unchanged.

**Gap costs must be on the matrix's scale.** The Gram transform puts a typical BLOSUM62
mismatch at ``SubstitutionMatrix.scale() == 14``; a ``gap_open`` of 1 would make gaps ~14x
cheaper than substitutions and every alignment would degenerate to gaps.

**A sequence score alone does not locate the gap.** Measured against 58 pairs of TCR-pMHC
crystal structures sharing peptide and MHC, the minimum-BLOSUM62 block position agreed with
the structurally correct one 8.6% of the time -- indistinguishable from picking at random
(8.6%). A central prior lifts that to 25.9% and cuts loop CA-RMSD from 2.15 A to 1.62 A
(oracle: 1.52 A). Pass ``gap_prior=central_prior(...)`` unless you have a better one.

**A prior is also what makes a column frame possible.** Pairwise-optimal gap placement is not
transitive: align A to B and B to C independently and the two column assignments do not
compose, so a set of unequal-length sequences has no consistent column index -- and hence no
profile. A *rule* that maps length to block position supplies one. Only a rule whose block
start is **constant in** ``d`` is transitive, i.e. one that pins the block to a fixed frame
column ``c``; see :func:`frame_prior` and :func:`embed_in_frame`. :func:`central_prior` is not
such a rule -- its block start drifts with ``d`` and the correspondence it induces between two
shorter members splits into two blocks.
"""
from __future__ import annotations

import collections
import math
from collections.abc import Callable, Iterable, Sequence

from ._core import Index, ScoreMatrix, SearchParams, SubstitutionMatrix
from ._core import gapblock_matrix as _gapblock_matrix

__all__ = [
    "gapblock_score", "score_matrix", "deletion_variants", "central_prior", "profile_prior",
    "frame_prior", "positions_prior", "embed_in_frame", "gap_cost", "GapBlockIndex", "ScoreMatrix",
    "IslandProfile",
]

#: A prior value large enough that no layout carrying it can ever win, yet small enough that
#: adding a real substitution cost cannot overflow int32.
UNREACHABLE = 1 << 20

#: ``prior(block_start, block_length, longer_length) -> int``. The block occupies columns
#: ``[i, i + d)`` of the longer sequence. Two requirements, both load-bearing: the result must
#: be ``>= 0`` (admissible trie pruning) and must be 0 when ``d == 0``, or ``s(q, q)`` stops
#: being zero and the score no longer defines a ball (``appendix/evalue.tex``).
#:
#: Monotonicity in ``d`` is *not* required and does not hold in general -- growing a leading
#: block drags its midpoint toward the centre, so :func:`central_prior` decreases.
GapPrior = Callable[[int, int, int], int]


def gap_cost(d: int, gap_open: int, gap_extend: int) -> int:
    """Affine cost of one gap block of length ``d``. ``d == 0`` costs nothing.

    The ``d == 0`` guard is not cosmetic: ``gap_open + (d-1)*gap_extend`` evaluates to
    ``gap_open - gap_extend`` at ``d == 0``, which is negative whenever
    ``gap_open < gap_extend``. A negative score would break ``s >= 0``, scope monotonicity
    and admissible trie pruning.
    """
    if d < 0:
        raise ValueError("gap length must be >= 0")
    return 0 if d == 0 else gap_open + (d - 1) * gap_extend


def central_prior(lam: int) -> GapPrior:
    """Penalise blocks whose midpoint sits away from the centre of the longer sequence.

    ``lam * abs(block_midpoint - m/2)``, where the block spans ``[i, i + d)``. ``lam ~ 1.5 *
    matrix.scale()`` reproduces the structurally-fitted optimum. Returns an integer so the
    total score stays an exact non-negative penalty.

    Not a transitive frame rule: its block start ``(m - d) // 2`` moves with ``d``.
    """
    if lam < 0:
        raise ValueError("lam must be >= 0")

    def prior(i: int, d: int, m: int) -> int:
        return 0 if d == 0 else int(lam * abs(2 * i + d - m) // 2)

    return prior


def profile_prior(lam: int, w: Callable[[int, int], float] | Sequence[float]) -> GapPrior:
    """Charge ``lam`` per unit of positional weight the block deletes.

    ``lam * sum(w(j, m) for j in block)``. With ``w(j, m)`` the probability that position ``j``
    of a length-``m`` sequence is germline-templated, this reads as *lam times the expected
    number of templated residues the gap had to remove* -- deleting conserved framework is
    implausible, deleting a non-templated insert is free.

    Args:
        lam: Cost per unit weight. Must be ``>= 0``.
        w: ``w(j, m) -> float`` in ``[0, 1]``, or a fixed sequence indexed by ``j`` when the
            frame has a fixed width. Prefer the callable: loop length varies.

    Returns:
        A :data:`GapPrior`. Zero at ``d == 0`` (the sum is empty) and non-negative. Unlike
        :func:`central_prior` it is also monotone non-decreasing in ``d``: a longer block can
        only delete more weight.

    Raises:
        ValueError: If ``lam`` is negative.
    """
    if lam < 0:
        raise ValueError("lam must be >= 0")
    wfn = w if callable(w) else (lambda j, m: w[j])

    def prior(i: int, d: int, m: int) -> int:
        return int(lam * sum(wfn(j, m) for j in range(i, i + d)))

    return prior


def frame_prior(lam: int, c: int) -> GapPrior:
    """Pin the block to frame column ``c``: ``lam * abs(i - c)``.

    The block start does not depend on ``d``, which makes this the **only** kind of rule under
    which embedding two sequences into a common frame reproduces their pairwise single-block
    alignment (see :func:`embed_in_frame`). Equivalent to left-anchoring the first ``c``
    residues and right-anchoring the rest.

    A large ``lam`` makes the pin hard: exactly one layout survives.
    """
    if lam < 0:
        raise ValueError("lam must be >= 0")
    if c < 0:
        raise ValueError("c must be >= 0")

    def prior(i: int, d: int, m: int) -> int:
        return 0 if d == 0 else lam * abs(i - c)

    return prior


def positions_prior(starts: Iterable[int]) -> GapPrior:
    """Allow the block to open only at ``starts``; let the score choose among them.

    A non-negative start counts from the sequence's beginning, a negative one from the end of
    the *shorter* sequence -- so ``(3, 4, -4, -3)`` reproduces the fixed gap set that
    ``mir.distances.aligner.JunctionAligner`` hardcodes for every locus. Starts outside
    ``[0, shorter]`` clamp into range, so at least one layout always survives.

    This is the "score several candidate placements and keep the best" rule. It is weaker than
    it looks: measured on human TRB retrieval at a matched false-positive rate, candidates
    ``(3, 4, mid)`` reached precision 0.156 against 0.414 for a single hard-pinned centre. The
    freer the placement, the more readily an unrelated reference manufactures a low score.

    Args:
        starts: Permitted block starts. Negative values index from the end.

    Returns:
        A :data:`GapPrior` returning 0 at a permitted start and :data:`UNREACHABLE` elsewhere.

    Raises:
        ValueError: If ``starts`` is empty.
    """
    ss = tuple(starts)
    if not ss:
        raise ValueError("need at least one start position")

    def prior(i: int, d: int, m: int) -> int:
        if d == 0:
            return 0
        shorter = m - d
        for p in ss:
            allowed = min(p, shorter) if p >= 0 else max(0, shorter + p)
            if i == allowed:
                return 0
        return UNREACHABLE

    return prior


def embed_in_frame(seq: str, width: int, c: int, gap: str = "-") -> str:
    """Place ``seq`` into a ``width``-column frame, gaps blocked at column ``c``.

    Columns ``[0, c)`` hold the sequence's own prefix (left-anchored) and columns
    ``[c + d, width)`` its suffix (right-anchored), with ``d = width - len(seq)`` gap columns
    between them. Applying this to every member of a set yields a multiple alignment whose
    columns are consistent -- which is what a position weight matrix needs.

    Example:
        >>> # c = 4: the V-templated CASS stays left, the J-templated EQYF stays right.
        >>> for s in ("CASSLGQGAYEQYF", "CASSLGQAYEQYF", "CASSGQAYEQYF"):
        ...     print(embed_in_frame(s, 14, 4))
        CASSLGQGAYEQYF
        CASS-LGQAYEQYF
        CASS--GQAYEQYF
    """
    d = width - len(seq)
    if d < 0:
        raise ValueError(f"sequence of length {len(seq)} does not fit a width-{width} frame")
    if not 0 <= c <= len(seq):
        raise ValueError(f"frame column {c} outside [0, {len(seq)}]")
    return seq[:c] + gap * d + seq[c:]


_AA = "ACDEFGHIKLMNPQRSTVWY"


class IslandProfile:
    """A position weight matrix over one island, scored as a non-negative penalty.

    Once a set of related sequences has been embedded into a common frame (see
    :func:`embed_in_frame`), each column has a residue distribution and a query can be scored
    column by column instead of against every member. The column penalty is measured **against the
    column's own consensus**::

        pen(j, a) = round(lam * log(p_max_j / p_j(a)))

    which is what keeps the score usable. A textbook PWM log-odds score is signed; this one is
    ``>= 0`` and exactly ``0`` on the consensus sequence, so it still defines a ball -- centred on
    the consensus rather than on any one member -- and still flows through
    :func:`seqtree.thetas_from_scores`.

    The gap is a column symbol like any other, so a column never gapped in the training members
    charges heavily for a gap there. There is no separate affine gap term: the island's own members
    say where a gap is tolerated.

    **When this is worth it depends on your cutoff, and there are two regimes.** The E-value's
    ``k = floor(e_target * M / N)`` is the number of control neighbours the cutoff may admit, so
    the false-positive rate is ``k / M`` and it moves with ``N``, the size of the set you annotate.

    * Building islands *within* one epitope group puts ``N`` at the group size (median 88), so
      ``k`` has median 142 of ``M = 250,000``: FPR ~ 5.7e-4.
    * Annotating a whole repertoire against known islands puts ``N`` at ~20,000. Then
      ``e_target = 0.05`` gives ``k = 0``, which :func:`seqtree.thetas_from_scores` reports as
      ``-1``: the rule of three certifies no ``E`` below ``3N/M = 0.236``. At that smallest
      certifiable ``E``, ``k = 3``: FPR ~ 1.2e-5.

    Recall on held-out members of 108 calibrated VDJdb islands of >= 10 (human TRB, three splits
    each, paired bootstrap over islands, 250,000 control negatives):

    ===============  ==========  ==================  ==============  ======================
    regime           FPR         min-over-members    IslandProfile   difference [95% CI]
    ===============  ==========  ==================  ==============  ======================
    loose reference  1%          **99.5 %**          99.1 %          -0.40 [-1.09, +0.14]
    per-epitope      0.0568%     88.3 %              **89.3 %**      +0.93 [-0.80, +2.79]
    repertoire       0.0012%     37.6 %              **48.5 %**      +10.90 [+7.69, +14.21]
    ===============  ==========  ==================  ==============  ======================

    So: **no significant difference when you are building the islands**, and a large one when you
    use them to annotate a repertoire. On islands of >= 50 members the repertoire-regime gap is
    9.8 % against 22.6 %.

    It does **not** generalise. Junctions specific to the same epitope that fall in a *different*
    island are recovered by neither this nor min-over-members (3.5 % vs 3.7 % at a 1 % FPR, by
    neither at either operating point). Distinct islands share no motif either representation finds.

    Nor is it a compression: 14 columns x 21 symbols x 4 B is 1,176 B against 182 B of member
    strings. An island needs 84 members before the profile is the smaller of the two, which 3.7 %
    of real islands reach.

    Args:
        penalties: One dict per frame column, mapping symbol to a non-negative integer penalty.
        width: Frame width; equals ``len(penalties)``.
        c: Frame column where the gap block opens.

    Example:
        >>> members = ["CASSLGQAYEQYF", "CASSLGQGYEQYF", "CASSLGQAYEQYF"]
        >>> p = IslandProfile.fit(members)
        >>> p.score(p.consensus())
        0
        >>> p.score("CASSLGQAYEQYF") <= p.score("CASSLGQGYEQYF")
        True
    """

    def __init__(self, penalties: list[dict[str, int]], width: int, c: int) -> None:
        self.penalties = penalties
        self.width = width
        self.c = c

    @classmethod
    def fit(
        cls,
        members: Sequence[str],
        c: int | None = None,
        lam: int = 1000,
        pseudocount: float = 0.5,
        gap: str = "-",
    ) -> IslandProfile:
        """Fit a profile to an island's members.

        Args:
            members: The island. Must be non-empty. The frame width is the longest member.
            c: Frame column for the gap block. ``None`` picks the column minimising summed column
                entropy -- the frame the members themselves prefer. On real islands the mode lands
                at 6, where crystal structures put the block.
            lam: Score resolution. Scores are compared against a control-calibrated cutoff, so any
                monotone rescaling cancels; ``lam`` only controls integer rounding.
            pseudocount: Added to every symbol count, so an unseen residue is expensive but finite.
            gap: The gap symbol used in the frame.

        Returns:
            A fitted :class:`IslandProfile`.

        Raises:
            ValueError: If ``members`` is empty, ``lam`` is negative, ``pseudocount`` is not
                positive, or an explicit ``c`` exceeds the shortest member's length.
        """
        members = list(members)
        if not members:
            raise ValueError("cannot fit a profile to an empty island")
        if lam < 0:
            raise ValueError("lam must be >= 0")
        if pseudocount <= 0:
            raise ValueError("pseudocount must be > 0")

        width = max(len(s) for s in members)
        shortest = min(len(s) for s in members)
        if c is None:
            c = min(range(shortest + 1), key=lambda k: cls._entropy(members, width, k, gap))
        elif not 0 <= c <= shortest:
            raise ValueError(f"frame column {c} outside [0, {shortest}]")

        alpha = _AA + gap
        # Cap a single column so a full-width score can never reach the UNREACHABLE sentinel.
        cap = UNREACHABLE // (width + 1)
        penalties = []
        for j in range(width):
            counts = collections.Counter(embed_in_frame(s, width, c, gap)[j] for s in members)
            denom = sum(counts.values()) + pseudocount * len(alpha)
            probs = {a: (counts.get(a, 0) + pseudocount) / denom for a in alpha}
            p_max = max(probs.values())
            penalties.append(
                {a: min(cap, int(round(lam * math.log(p_max / p)))) for a, p in probs.items()}
            )
        return cls(penalties, width, c)

    @staticmethod
    def _entropy(members: Sequence[str], width: int, c: int, gap: str) -> float:
        total = 0.0
        for j in range(width):
            counts = collections.Counter(embed_in_frame(s, width, c, gap)[j] for s in members)
            n = sum(counts.values())
            total -= sum((v / n) * math.log(v / n) for v in counts.values() if v)
        return total

    def consensus(self, gap: str = "-") -> str:
        """The zero-penalty sequence: each column's most frequent symbol, gaps stripped.

        This is the centre of the ball the profile defines. ``score(consensus()) == 0``.
        """
        out = "".join(min(col, key=col.get) for col in self.penalties)
        return out.replace(gap, "")

    def score(self, seq: str, gap: str = "-") -> int:
        """Penalty of ``seq`` in this island's frame; ``0`` on the consensus.

        A sequence that does not fit the frame -- longer than ``width``, or shorter than ``c`` --
        cannot be embedded and scores :data:`UNREACHABLE`. That is a rejection, not an error: it
        must still count as a scored sequence when a cutoff is calibrated against a control.
        """
        if len(seq) > self.width or self.c > len(seq):
            return UNREACHABLE
        emb = embed_in_frame(seq, self.width, self.c, gap)
        cap = UNREACHABLE // (self.width + 1)
        return sum(self.penalties[j].get(ch, cap) for j, ch in enumerate(emb))

    def score_batch(self, seqs: Iterable[str], gap: str = "-") -> list[int]:
        """:meth:`score` for many sequences. Score the whole control this way, then hand the
        result to :func:`seqtree.thetas_from_scores` for a calibrated cutoff."""
        return [self.score(s, gap) for s in seqs]

    def __repr__(self) -> str:
        return f"IslandProfile(width={self.width}, c={self.c})"


def deletion_variants(q: str, d: int) -> list[tuple[int, str]]:
    """Every contiguous-``d``-deletion variant of ``q``, as ``(block_position, variant)``.

    Variant ``i`` is ``q[:i] + q[i+d:]``; it is exactly the query as seen through a gap block
    opened at position ``i``. ``d == 0`` has no block, so it yields a single identity variant
    (not ``len(q) + 1`` copies of it).

    Example:
        >>> deletion_variants("CAST", 1)
        [(0, 'AST'), (1, 'CST'), (2, 'CAT'), (3, 'CAS')]
        >>> deletion_variants("CAST", 0)
        [(0, 'CAST')]
    """
    if d < 0 or d > len(q):
        raise ValueError(f"cannot delete {d} residues from a length-{len(q)} sequence")
    if d == 0:
        return [(0, q)]
    return [(i, q[:i] + q[i + d:]) for i in range(len(q) - d + 1)]


def _pen_table(matrix: SubstitutionMatrix | None, symbols: str) -> dict[tuple[str, str], int]:
    if matrix is None:
        return {(a, b): int(a != b) for a in symbols for b in symbols}
    return {(a, b): matrix.penalty(a, b) for a in symbols for b in symbols}


def _default_gap_open(matrix: SubstitutionMatrix | None) -> int:
    return 2 * matrix.scale() if matrix is not None else 1


def gapblock_score(
    q: str,
    r: str,
    matrix: SubstitutionMatrix | None = None,
    gap_open: int | None = None,
    gap_extend: int = 1,
    gap_prior: GapPrior | None = None,
    _pen: dict[tuple[str, str], int] | None = None,
) -> tuple[int, int]:
    """Optimal single-gap-block alignment score of ``q`` against ``r``.

    Returns ``(score, block_position)``. ``score`` is a non-negative penalty, zero iff the
    sequences are identical. ``block_position`` indexes the shorter sequence and is
    inclusive at both ends: 0 is a leading block, ``min(m, n)`` a trailing one.

    The gap block sits in whichever sequence is shorter, so this is symmetric in ``q``/``r``.
    Being a strict restriction of affine alignment, the score is always ``>=`` the
    unrestricted Gotoh optimum; on pairs differing by a pure length change the two agree
    exactly, and at a calibrated ``gap_open`` they agree on 90% of pairs that also carry
    substitutions.

    Args:
        q: Query sequence.
        r: Reference sequence.
        matrix: Substitution penalties. ``None`` means unit cost (1 per mismatch).
        gap_open: Cost of opening the block. Defaults to ``2 * matrix.scale()``, or 1 for
            unit cost. Must be ``>= 0``.
        gap_extend: Cost of each additional gap column. Must be ``>= 0``.
        gap_prior: :data:`GapPrior`, added to each candidate position. See
            :func:`central_prior`, :func:`profile_prior`, :func:`frame_prior`. ``None``
            disables it. It applies only when there is a block to place (``d > 0``); otherwise
            ``s(q, q)`` would be non-zero and the score would no longer define a ball.
        _pen: Precomputed penalty lookup, for hot loops.

    Returns:
        ``(score, block_position)``.

    Raises:
        ValueError: If a gap cost is negative.

    Example:
        >>> m = SubstitutionMatrix.blosum62()
        >>> gapblock_score("CASSLGQAYEQYF", "CASSLGQAYEQYF", m)
        (0, 0)
    """
    if gap_open is None:
        gap_open = _default_gap_open(matrix)
    if gap_open < 0 or gap_extend < 0:
        raise ValueError("gap_open and gap_extend must be >= 0")

    m, n = len(q), len(r)
    d = abs(m - n)
    shorter = min(m, n)
    pen = _pen if _pen is not None else _pen_table(matrix, "".join(sorted(set(q + r))))

    # A[i] = penalty of q[:i] vs r[:i] on diagonal 0 (before the block).
    prefix = [0] * (shorter + 1)
    for j in range(shorter):
        prefix[j + 1] = prefix[j] + pen[q[j], r[j]]

    # B[i] = penalty of the tails on diagonal d (after the block). The offset sits on the
    # LONGER sequence's index -- getting this backwards silently misaligns the suffix.
    suffix = [0] * (shorter + 1)
    for j in range(shorter - 1, -1, -1):
        step = pen[q[j + d], r[j]] if m >= n else pen[q[j], r[j + d]]
        suffix[j] = suffix[j + 1] + step

    # No block to place when the lengths match, so no positional cost -- otherwise s(q, q)
    # would pick up the prior's minimum and stop being zero.
    longer = max(m, n)
    prior = gap_prior if (gap_prior is not None and d > 0) else (lambda i, d_, m_: 0)
    best_i, best = 0, None
    for i in range(shorter + 1):
        cand = prefix[i] + suffix[i] + prior(i, d, longer)
        if best is None or cand < best:
            best, best_i = cand, i
    return best + gap_cost(d, gap_open, gap_extend), best_i


def _prior_cube(prior: GapPrior, width: int) -> list[int]:
    """Flatten ``prior`` to ``[m][d][i]`` with stride ``width + 1``, for the C++ kernel.

    Only reachable cells are filled: ``d >= 1`` (the prior never applies to equal lengths) and
    ``i <= m - d`` (the block start indexes the shorter sequence). The rest stay zero and are
    never read.
    """
    w1 = width + 1
    cube = [0] * (w1 * w1 * w1)
    for m in range(w1):
        for d in range(1, m + 1):
            base = (m * w1 + d) * w1
            for i in range(m - d + 1):
                v = int(prior(i, d, m))
                if v < 0:
                    raise ValueError(f"gap_prior({i}, {d}, {m}) = {v} is negative")
                cube[base + i] = v
    return cube


def score_matrix(
    queries: Sequence[str],
    refs: Sequence[str],
    matrix: SubstitutionMatrix | None = None,
    gap_open: int | None = None,
    gap_extend: int = 1,
    gap_prior: GapPrior | None = None,
    alphabet: str = "aa",
    threads: int = 0,
) -> ScoreMatrix:
    """Gap-block penalty of every query against every reference, in C++ with the GIL released.

    The exhaustive counterpart of :meth:`GapBlockIndex.search`: no budget, no trie, every cell
    scored. This is the shape a prototype-distance embedding wants -- ``n`` clonotypes against a
    few thousand fixed references -- and it is where :func:`gapblock_score` stops being fast
    enough, at roughly 0.4 M pairs/s in Python against ~50 M in the kernel.

    Args:
        queries: Query sequences (the matrix rows).
        refs: Reference sequences (the matrix columns).
        matrix: Substitution penalties; ``None`` means unit cost.
        gap_open: Block-opening cost. Defaults to ``2 * matrix.scale()``. See the module note:
            leaving this at 1 with a real matrix makes gaps ~14x cheaper than substitutions.
        gap_extend: Cost per additional gap column.
        gap_prior: :data:`GapPrior`, materialized once into a lookup cube and then read from C++.
            ``None`` lets the score alone choose the block position.
        alphabet: ``"aa"``, ``"nt"``, or ``"iupac"``. Symbols outside it raise.
        threads: Worker threads; ``0`` means one per core. Rows are disjoint, so this scales.

    Returns:
        A :class:`ScoreMatrix` of shape ``(len(queries), len(refs))``. It holds
        ``4 * len(queries) * len(refs)`` bytes -- 1.2 GB at 100k x 3000 -- so chunk the queries
        if that does not fit. ``numpy.asarray`` wraps it without copying.

    Raises:
        ValueError: If a gap cost is negative, or the prior returns a negative value.

    Example:
        >>> m = SubstitutionMatrix.blosum62()
        >>> sm = score_matrix(["CASSLGQAYEQYF"], ["CASSLGQAYEQYF", "CASSLGAYEQYF"], m)
        >>> sm.shape
        (1, 2)
        >>> sm[0, 0]
        0
    """
    q, r = list(queries), list(refs)
    if gap_open is None:
        gap_open = _default_gap_open(matrix)
    if gap_open < 0 or gap_extend < 0:
        raise ValueError("gap_open and gap_extend must be >= 0")
    width = max((len(s) for s in q + r), default=0)
    cube = _prior_cube(gap_prior, width) if gap_prior is not None else []
    return _gapblock_matrix(q, r, alphabet, matrix, gap_open, gap_extend, cube, width, threads)


class GapBlockIndex:
    """Search a reference set under the single-gap-block model, reusing the Hamming engine.

    Refs *shorter* than the query are reached by Hamming-matching the query's deletion
    variants against the ordinary index (the Hamming path only ever terminates on refs of
    the query's length, so no length partitioning is needed). Refs *longer* than the query
    are reached by pre-indexing the refs' own deletion variants, one auxiliary index per
    block length.

    Building costs ``O(d_max * total_residues)`` extra index entries -- roughly 14x the base
    index for CDR3 at ``d_max=1``. Build once, query many.

    Profiled over the bundled 250k control at ``d_max=3``: 91% of query time is the *first*
    branch (the query's own deletion variants against the base index) and only 9% the auxiliary
    indices, despite those holding 9.8M entries. Netting the prior out of each variant's budget
    already cuts that branch from ~15 sub-searches per query to 2.5. Deduplicating variants
    would touch 7-10% of them before pruning and fewer after, and bucketing the auxiliary
    indices by reference length saves no memory -- neither is worth the code.
    """

    def __init__(self, refs: Iterable[str], alphabet: str = "aa", d_max: int = 1):
        if d_max < 0:
            raise ValueError("d_max must be >= 0")
        self.refs = list(refs)
        self.d_max = d_max
        self.alphabet = alphabet
        self._base = Index.build(list(self.refs), alphabet)
        # For each d >= 1: an index over every d-deletion variant of every ref, plus the
        # (ref_id, block_position) each variant came from.
        self._var: dict[int, tuple[Index, list[tuple[int, int]]]] = {}
        for d in range(1, d_max + 1):
            variants: list[str] = []
            owner: list[tuple[int, int]] = []
            for rid, r in enumerate(self.refs):
                if len(r) < d:
                    continue
                for i, v in deletion_variants(r, d):
                    variants.append(v)
                    owner.append((rid, i))
            self._var[d] = (Index.build(variants, alphabet), owner) if variants else (None, [])

    def __len__(self) -> int:
        return len(self.refs)

    def search(
        self,
        query: str,
        max_penalty: int,
        matrix: SubstitutionMatrix | None = None,
        gap_open: int | None = None,
        gap_extend: int = 1,
        gap_prior: GapPrior | None = None,
    ) -> list[tuple[int, int, int, int]]:
        """All refs within ``max_penalty`` under the gap-block score.

        Returns ``(ref_id, score, block_length, block_position)`` per ref, best score kept,
        sorted by ascending score. ``block_position`` indexes the shorter of the two.
        """
        if gap_open is None:
            gap_open = _default_gap_open(matrix)
        if max_penalty < 0:
            raise ValueError("max_penalty must be >= 0")
        zero: GapPrior = lambda i, d, m: 0  # noqa: E731
        mat = matrix if matrix is not None else ""
        best: dict[int, tuple[int, int, int]] = {}

        def offer(rid: int, score: int, d: int, i: int) -> None:
            if score <= max_penalty and (rid not in best or score < best[rid][0]):
                best[rid] = (score, d, i)

        for d in range(0, self.d_max + 1):
            g = gap_cost(d, gap_open, gap_extend)
            budget = max_penalty - g
            if budget < 0:
                continue
            prior = gap_prior if (gap_prior is not None and d > 0) else zero

            # Refs shorter than the query (the block sits in the query, the longer sequence):
            # delete d from the query.
            if d <= len(query):
                for i, qv in deletion_variants(query, d):
                    extra = prior(i, d, len(query))
                    if budget - extra < 0:
                        continue
                    for h in self._base.search(qv, self._params(budget - extra, matrix=mat, qlen=len(qv))):
                        offer(h.ref_id, h.score + g + extra, d, i)

            # Refs longer than the query (the block sits in the ref, of length len(query) + d):
            # the refs' variants were pre-deleted, so match the full query against them.
            if d >= 1:
                idx, owner = self._var[d]
                if idx is None:
                    continue
                for h in idx.search(query, self._params(budget, matrix=mat, qlen=len(query))):
                    rid, i = owner[h.ref_id]
                    offer(rid, h.score + g + prior(i, d, len(query) + d), d, i)

        return sorted(((rid, s, d, i) for rid, (s, d, i) in best.items()), key=lambda t: (t[1], t[0]))

    @staticmethod
    def _params(budget: int, matrix, qlen: int) -> SearchParams:
        # budget == 0 must mean "exact only". Passing max_penalty=0 would instead be read as
        # "no explicit budget" and, in matrix mode, resolve to an unbounded ball -- so cap the
        # substitution count to 0 in that case and let the caps do the work.
        if budget <= 0:
            return SearchParams(max_subs=0, matrix=matrix, engine="seqtm")
        return SearchParams(max_subs=qlen, max_penalty=budget, matrix=matrix, engine="seqtm")
