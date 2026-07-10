#!/usr/bin/env python
"""Gapped search: find neighbours that differ in length, and know when to believe them.

An ordinary Hamming ball can only ever reach references of the query's own length. Half of the
co-specific neighbours of a TCR junction are not that length. ``GapBlockIndex`` reaches them by
allowing exactly one contiguous indel, and ``threshold_for_evalue`` decides how close is close
enough -- per query, because a fixed cutoff is not a calibrated one.

Self-contained: the only data used is the 250k control repertoire bundled with seqtree, plus a
small family planted inside it.

Run:
    python examples/01_gapped_search.py
"""
import random

import seqtree
from seqtree.gapblock import GapBlockIndex, central_prior

HUB = "CASSLGQAYEQYF"
AA = "ACDEFGHIKLMNPQRSTVWY"


def plant_family(hub, mat, rng, n=12):
    """Relatives of ``hub``: some substituted, some one or two residues shorter or longer.

    Half the substitutions are conservative (a residue the matrix considers cheap), because real
    relatives are not uniform random mutants and the distinction is what a matrix is for.
    """
    cheap = {a: [b for b in AA if b != a and mat.penalty(a, b) <= 6] for a in AA}
    fam = []
    for k in range(n):
        s = list(hub)
        for pos in rng.sample(range(4, len(hub) - 4), rng.choice([1, 1, 2])):
            pool = cheap[s[pos]] if (k % 2 == 0 and cheap[s[pos]]) else AA
            s[pos] = rng.choice(pool)
        if k % 3 == 1:                                   # delete a residue from the core
            cut = rng.randrange(5, len(s) - 5)
            del s[cut]
        elif k % 3 == 2:                                 # insert one
            s.insert(rng.randrange(5, len(s) - 4), rng.choice(AA))
        fam.append("".join(s))
    return sorted(set(fam))


def main():
    rng = random.Random(0)
    control = seqtree.load_control("human_trb_aa")
    background = [control.ref_seq(i) for i in range(len(control))]

    mat = seqtree.SubstitutionMatrix.blosum62()
    family = plant_family(HUB, mat, rng)
    # N must stay well under M/3 or no cutoff can certify E* < 3N/M -- see the end of this script.
    refs = background[:2_000] + family                   # the haystack, plus the needles
    print(f"planted {len(family)} relatives of {HUB} among {len(refs):,} references")
    for f in family:
        d = len(f) - len(HUB)
        print(f"    {f:<16} {'same length' if d == 0 else f'{d:+d} residue'}")

    gap_open = 2 * mat.scale()                           # 28. NEVER leave this at its default of 1.
    prior = central_prior(int(1.5 * mat.scale()))

    # ---- 1. a Hamming ball cannot see the length-different relatives -----------------------
    hamming = seqtree.Index.build(refs, "aa")
    p = seqtree.SearchParams(max_subs=len(HUB), max_penalty=40, matrix=mat, engine="seqtm")
    ham_hits = {refs[h.ref_id] for h in hamming.search(HUB, p)}
    print(f"\nHamming ball (max_penalty 40): {len(ham_hits)} hits, "
          f"{sum(1 for h in ham_hits if len(h) != len(HUB))} of them length-different")

    # ---- 2. the gap-block ball reaches them ------------------------------------------------
    gbi = GapBlockIndex(refs, "aa", d_max=2)
    hits = gbi.search(HUB, 60, mat, gap_open=gap_open, gap_prior=prior)
    gapped = [(refs[i], s, d, pos) for i, s, d, pos in hits if d]
    print(f"gap-block ball  (max_penalty 60): {len(hits)} hits, "
          f"{len(gapped)} of them length-different")

    # An indel costs gap_open = 28, so gapped hits never top a score-sorted list. Show them.
    print(f"\n  best length-different hits (invisible to the Hamming ball):")
    print(f"  {'reference':<16}{'score':>7}{'block len':>11}{'block pos':>11}{'planted':>9}")
    for r, score, d, pos in gapped[:8]:
        print(f"  {r:<16}{score:>7}{d:>11}{pos:>11}{'yes' if r in family else '':>9}")
    print(f"\n  of the {len(family)} planted relatives, "
          f"{sum(1 for r, *_ in gapped if r in family)} length-different ones were recovered.")

    # ---- 3. how close is close enough? ------------------------------------------------------
    # theta is derived from the control, not chosen. One control scan serves every query.
    target = seqtree.Index.build(refs, "aa")
    N, M = len(target), len(control)
    ceiling = seqtree.SearchParams(max_subs=len(HUB), max_penalty=50, matrix=mat, engine="seqtm")
    queries = [HUB, background[len(background) // 2]]      # a hub, and an ordinary junction
    print(f"\nN = {N:,} references, M = {M:,} control. Rule of three: no E* below "
          f"3N/M = {3 * N / M:.3f} can be certified.")

    for e in (1.0, 0.05, 0.001):
        thetas = seqtree.threshold_for_evalue(target, control, queries, ceiling,
                                              e_target=e, exclude_exact=True)
        shown = "  ".join(f"{q[:6]}..={t if t >= 0 else 'unreachable'}"
                          for q, t in zip(queries, thetas))
        print(f"  E* = {e:<6} theta per query:  {shown}")

    print("\n  The two queries get DIFFERENT cutoffs. CASSLGQAYEQYF sits in a dense, near-germline")
    print("  neighbourhood, so the same score buys it far more chance neighbours; it is held to a")
    print("  stricter bar. That is the entire reason the cutoff is per query.")

    print(f"\n  {'E*':>7}{'theta':>7}{'hits':>7}{'planted':>9}{'background':>12}{'precision':>11}")
    for e in (1.0, 0.05):
        th = seqtree.threshold_for_evalue(target, control, [HUB], ceiling, e_target=e,
                                          exclude_exact=True)[0]
        kept = [refs[i] for i, s, _, _ in hits if s <= th and refs[i] != HUB]
        tp = sum(1 for r in kept if r in family)
        prec = tp / len(kept) if kept else float("nan")
        print(f"  {e:>7}{th:>7}{len(kept):>7}{tp:>9}{len(kept) - tp:>12}{prec:>11.2f}")
    print("\n  Tightening E* trades recall for precision, on a scale that means something:")
    print("  'at most one chance neighbour per twenty queries' rather than 'score below 40'.")


if __name__ == "__main__":
    main()
