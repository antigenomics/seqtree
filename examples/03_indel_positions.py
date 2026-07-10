#!/usr/bin/env python
"""Where does the gap go? Preferred indel positions, and what a prior actually assumes.

Two sequences of different length can be aligned with the gap block anywhere. ``gapblock_score``
scores every position and keeps the best -- but "best" can mean best *score*, or best score plus a
*prior*, and those disagree on about 90% of pairs.

This example plants a deletion at a known position and asks each rule to find it again. It then
shows the cost of the prior's assumption by planting the deletion somewhere the prior does not
expect. A prior is not free: it is a claim about where indels occur, and you should be able to
state it and to measure what it costs when it is wrong.

Self-contained: only the 250k control repertoire bundled with seqtree.

Run:
    python examples/03_indel_positions.py
"""
import collections
import random

import seqtree
from seqtree.gapblock import central_prior, gapblock_score

HUGE = 10 ** 6


def plant(seq, pos, d):
    """Delete ``d`` residues of ``seq`` starting at ``pos``. The block position is then ``pos``."""
    return seq[:pos] + seq[pos + d:]


def rules(lam):
    return [
        ("fixed centre", central_prior(HUGE)),   # one layout; the score never votes
        (f"central prior lam={lam}", central_prior(lam)),
        ("flat (score alone)", None),            # every layout, lowest score wins
    ]


AA = "ACDEFGHIKLMNPQRSTVWY"


def make_pairs(pool, rng, n, where, subs=0):
    """``(long, short, planted_position)``: one residue deleted at a known place, then ``subs``
    substitutions, because real relatives are never a pure indel."""
    out = []
    while len(out) < n:
        s = pool[rng.randrange(len(pool))]
        if where == "core":
            p = rng.randrange(4, len(s) - 5)
        else:
            p = rng.choice([1, 2, len(s) - 3, len(s) - 2])
        # Deleting from inside a homopolymer run is positionally ambiguous: several p give the
        # same shorter string, so no rule can be graded on it. Skip those.
        r = plant(s, p, 1)
        if any(plant(s, k, 1) == r for k in range(len(s)) if k != p):
            continue
        if subs:
            t = list(r)
            for pos in rng.sample(range(len(t)), subs):
                t[pos] = rng.choice(AA)
            r = "".join(t)
        out.append((s, r, p))
    return out


def recover(pairs, mat, gap_open, lam):
    """For each rule: how often does its chosen block position equal the planted one?"""
    out = {}
    for name, prior in rules(lam):
        exact = near = 0
        for q, r, true_pos in pairs:
            _, i = gapblock_score(q, r, mat, gap_open, 1, prior)
            exact += i == true_pos
            near += abs(i - true_pos) <= 1
        out[name] = (exact / len(pairs), near / len(pairs))
    return out


def main():
    rng = random.Random(0)
    control = seqtree.load_control("human_trb_aa")
    pool = [control.ref_seq(i) for i in range(len(control))]
    pool = [s for s in pool if 13 <= len(s) <= 16]

    mat = seqtree.SubstitutionMatrix.blosum62()
    gap_open, lam = 2 * mat.scale(), int(1.5 * mat.scale())
    n = 3000

    print(f"{n:,} pairs per cell: a control junction with ONE residue deleted at a known place,")
    print("then k substitutions -- because a real relative is never a pure indel. Homopolymer")
    print("deletions are skipped; several positions give the same string, so no rule can be graded.")
    print(f"BLOSUM62 Gram scale = {mat.scale()}, gap_open = {gap_open}, lam = {lam}\n")

    # ---- 1. the score can only find the gap when nothing else is happening -------------------
    print("=== deletion planted in the CORE (positions 4 .. L-5): exact recovery ===")
    print(f"  {'rule':<26}" + "".join(f"{f'k={k} subs':>11}" for k in (0, 1, 2, 3)))
    for name, prior in rules(lam):
        cells = []
        for k in (0, 1, 2, 3):
            pairs = make_pairs(pool, random.Random(100 + k), n, "core", subs=k)
            exact = sum(gapblock_score(q, r, mat, gap_open, 1, prior)[1] == p for q, r, p in pairs)
            cells.append(f"{100*exact/len(pairs):>10.1f}%")
        print(f"  {name:<26}" + "".join(cells))

    print("\n  Read this carefully, because it is easy to draw the wrong conclusion. When two")
    print("  sequences REALLY ARE one deletion plus a little noise, the score finds the gap:")
    print("  100% at k = 0, still 76% at k = 3. The prior, which never reads the sequence, sits")
    print("  near 26% throughout -- it is right exactly as often as the planted position happens")
    print("  to be the centre.")
    print("\n  So why constrain the block at all? Because the ground truth above is a CONSTRUCTION.")
    print("  Against TCR-pMHC crystal structures -- where 'correct' means 'the residues that")
    print("  actually superpose' -- the score agrees about a tenth of the time and a central prior")
    print("  30-42%. And in retrieval at a matched false-positive rate, letting the score choose")
    print("  halves precision. The score finds the alignment you planted. It does not find the")
    print("  alignment biology used, and on unrelated sequences it invents one (see below).\n")

    # ---- 2. the prior's assumption, priced ---------------------------------------------------
    print("=== deletion planted at the EDGE (within 2 of an anchor), k = 1 sub ===")
    print(f"  {'rule':<26}{'exact':>9}{'within 1':>11}")
    for name, (ex, nr) in recover(make_pairs(pool, random.Random(7), n, "edge", subs=1),
                                  mat, gap_open, lam).items():
        print(f"  {name:<26}{100*ex:>8.1f}%{100*nr:>10.1f}%")

    print("\n  The core is where V/J trimming and N-addition put the indel, and where a central")
    print("  prior is right. Move the deletion to the anchors and the prior is confidently wrong.")
    print("  That gap is the prior's assumption, priced. State it; do not inherit it silently.\n")

    # ---- 2. the marginal: where does each rule LAND on ordinary pairs? -----------------------
    print("=== marginal distribution of the chosen block position, on unplanted pairs ===")
    print("  (random control junctions differing by one residue -- no ground truth, just the pile-up)")
    unplanted = []
    while len(unplanted) < 4000:
        a, b = pool[rng.randrange(len(pool))], pool[rng.randrange(len(pool))]
        if abs(len(a) - len(b)) == 1:
            unplanted.append((a, b))

    print(f"\n  {'rule':<26}{'i/L mean':>10}{'sd':>7}   histogram of i/L (0 = Cys end, 1 = Phe end)")
    for name, prior in rules(lam):
        rel = []
        for a, b in unplanted:
            _, i = gapblock_score(a, b, mat, gap_open, 1, prior)
            rel.append(i / min(len(a), len(b)))
        hist = collections.Counter(min(int(x * 10), 9) for x in rel)
        peak = max(hist.values())
        bars = "".join("#" if hist[k] > peak * 0.5 else ("+" if hist[k] > peak * 0.1 else ".")
                       for k in range(10))
        mean = sum(rel) / len(rel)
        sd = (sum((x - mean) ** 2 for x in rel) / len(rel)) ** 0.5
        print(f"  {name:<26}{mean:>10.3f}{sd:>7.3f}   {bars}")

    edges = 0
    for a, b in unplanted:
        _, i = gapblock_score(a, b, mat, gap_open, 1, None)
        L = min(len(a), len(b))
        edges += i <= 2 or i >= L - 2
    print(f"\n  Left free, the score puts the block within 2 residues of an anchor on "
          f"{100*edges/len(unplanted):.1f}% of these")
    print("  unrelated pairs. Those are the layouts where a spurious low score hides: shift the")
    print("  whole sequence by one and the conserved CASS / EQYF frame still lines up. That is")
    print("  the same mechanism by which unrestricted affine alignment undercuts the block score")
    print("  by a median of 106 penalty units on unrelated pairs -- it manufactures similarity.")
    print("\n  Measured against TCR-pMHC crystal structures, the true block sits at the loop apex")
    print("  (Cys-offset 6, both chains) and does not drift with block length. The central prior")
    print("  recovers it 42.4% (TRA) / 30.1% (TRB) of the time; the score alone, about a tenth.")


if __name__ == "__main__":
    main()
