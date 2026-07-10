#!/usr/bin/env python
"""A dendrogram over an island, built on the gap-block distance.

Two cautions before the code, both measured rather than assumed.

*Do not build one over a whole repertoire.* A dense condensed distance matrix over the ~90k CDR3
of a VDJdb shortlist is 32 GB. Cluster first (union-find on calibrated edges), then draw a tree
inside each island, where the member count is in the tens.

*Do not read anything into the merges above the calibrated cutoff.* Distinct sequence islands for
one epitope share no motif that similarity can find: only ~0.5% of cross-island co-specific pairs
share even a rare central 4-mer, and none share a 6-mer. The tree above theta is arithmetic, not
biology.

Average linkage over a handful of sequences is thirty lines, so this uses no scipy -- seqtree has
no runtime dependencies and this example keeps it that way. Output is Newick (paste into any tree
viewer) plus an ASCII rendering.

Run:
    python examples/02_sequence_dendrogram.py
"""
import random

import seqtree
from seqtree.gapblock import central_prior, gapblock_score

AA = "ACDEFGHIKLMNPQRSTVWY"


def island(hub, rng, n=14):
    """A family with internal structure: two sub-branches, one carrying a length change."""
    left = [hub]
    for _ in range(n // 2 - 1):
        s = list(rng.choice(left))
        s[rng.randrange(4, len(s) - 4)] = rng.choice(AA)
        left.append("".join(s))
    right = []
    trunk = list(hub)
    del trunk[rng.randrange(5, len(trunk) - 5)]           # the branch defined by one deletion
    right.append("".join(trunk))
    for _ in range(n - len(left) - 1):
        s = list(rng.choice(right))
        s[rng.randrange(4, len(s) - 4)] = rng.choice(AA)
        right.append("".join(s))
    return sorted(set(left + right))


def distances(seqs, mat, gap_open, prior, d_max=3):
    """Condensed gap-block distances. Pairs too far apart in length never merge."""
    pen = {(a, b): mat.penalty(a, b) for a in AA for b in AA}
    inf = 10 ** 6
    d = {}
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            a, b = seqs[i], seqs[j]
            if abs(len(a) - len(b)) > d_max:
                d[i, j] = inf
            else:
                d[i, j] = gapblock_score(a, b, gap_open=gap_open, gap_extend=1,
                                         gap_prior=prior, _pen=pen)[0]
    return d


def average_linkage(seqs, d):
    """Merge the closest two clusters, distance = mean over cross pairs. Returns a tree + heights."""
    def between(x, y):
        return sum(d[min(i, j), max(i, j)] for i in x for j in y) / (len(x) * len(y))

    clusters = {i: (i,) for i in range(len(seqs))}          # id -> member indices
    node = {i: seqs[i] for i in range(len(seqs))}           # id -> newick subtree
    height = {i: 0.0 for i in range(len(seqs))}
    nxt = len(seqs)
    merges = []
    while len(clusters) > 1:
        (a, b), h = min((((a, b), between(clusters[a], clusters[b]))
                         for a in clusters for b in clusters if a < b), key=lambda t: t[1])
        merges.append((node[a], node[b], h))
        node[nxt] = f"({node[a]}:{(h - height[a]) / 2:.1f},{node[b]}:{(h - height[b]) / 2:.1f})"
        clusters[nxt] = clusters.pop(a) + clusters.pop(b)
        height[nxt] = h
        del node[a], node[b]
        nxt += 1
    root = next(iter(node))
    return node[root] + ";", merges, height[root]


def ascii_tree(merges, width=52):
    """One line per merge, indented by height. Crude, dependency-free, and legible."""
    top = max(h for *_, h in merges) or 1.0
    print(f"  {'merge':<46}{'height':>8}")
    for a, b, h in merges:
        bar = "-" * int(width * h / top)
        la = a if len(a) < 18 else "(...)"
        lb = b if len(b) < 18 else "(...)"
        print(f"  {la:<20}+ {lb:<22}{h:>8.1f}  {bar}")


def main():
    rng = random.Random(1)
    hub = "CASSLGQAYEQYF"
    seqs = island(hub, rng)
    mat = seqtree.SubstitutionMatrix.blosum62()
    gap_open, lam = 2 * mat.scale(), int(1.5 * mat.scale())
    prior = central_prior(lam)

    print(f"island of {len(seqs)} members around {hub} "
          f"(lengths {min(map(len, seqs))}-{max(map(len, seqs))})")
    for s in seqs:
        print(f"    {s}")

    d = distances(seqs, mat, gap_open, prior)
    newick, merges, root_h = average_linkage(seqs, d)

    print(f"\naverage linkage on the gap-block distance (gap_open = {gap_open}, lam = {lam}):\n")
    ascii_tree(merges)

    print(f"\nNewick (root height {root_h:.1f}):\n{newick}")

    # Where would a calibrated cutoff cut this tree? The E-value's N is the size of the set you
    # are actually annotating, not of the island -- an island of 12 makes everything "significant".
    control = seqtree.load_control("human_trb_aa")
    annotation_set = [control.ref_seq(i) for i in range(20_000)] + seqs
    target = seqtree.Index.build(annotation_set, "aa")
    ceiling = seqtree.SearchParams(max_subs=20, max_penalty=60, matrix=mat, engine="seqtm")

    print(f"\nN = {len(target):,} (the annotation set), M = {len(control):,} (the control)")
    for e in (1.0, 0.1):
        thetas = seqtree.threshold_for_evalue(target, control, seqs, ceiling,
                                              e_target=e, exclude_exact=True)
        usable = [t for t in thetas if t >= 0]
        if not usable:
            print(f"  E* = {e}: unreachable (3N/M = {3 * len(target) / len(control):.3f})")
            continue
        cut = sum(usable) / len(usable)
        below = sum(1 for *_, h in merges if h <= cut)
        print(f"  E* = {e:<5} mean theta {cut:5.1f}   {below}/{len(merges)} merges survive the cut"
              f"   -> {len(seqs) - below} clusters")
    print("\n  Merges above the cut join sequences whose proximity the control already explains.")
    print("  Cut there. Do not interpret what lies above -- including the root, which here is")
    print(f"  the {gap_open}-unit gap cost separating the deletion branch from the rest.")


if __name__ == "__main__":
    main()
