#!/usr/bin/env python
"""A per-island profile, and the one number that decides whether it is worth building.

Given a set of related junctions -- an island -- you can score a new sequence two ways: against
every member, keeping the closest, or against a position weight matrix fitted to the island. The
second is a profile. Which one wins is not a matter of taste. It depends on how strict your cutoff
is, and the cutoff is not something you choose freely: the E-value fixes it.

``k = floor(e_target * M / N)`` is how many control neighbours a calibrated cutoff may admit, so
the false-positive rate is ``k / M``. It moves with ``N``, the size of the set you are annotating:

    building islands inside one epitope group   N = group size (~90)   ->  FPR ~ 6e-4
    annotating a whole repertoire against them  N ~ 20,000             ->  FPR ~ 1e-5

Those are two different questions and they get two different answers. At the first, a profile and
min-over-members are indistinguishable. At the second, on 108 real VDJdb islands, the profile
recovers 48.5% of held-out members against min-over-members' 37.6%.

This script shows the mechanism on one planted island. **One island cannot settle the question** --
the numbers above come from `vdjmatch/bench/island_pwm.py` over 108 of them. What it can show is
where the crossover lives and why the E-value, not taste, puts you on one side of it.

Self-contained: the 250k control repertoire bundled with seqtree, plus a family planted in it.

Run:
    python examples/04_island_profile.py
"""
import bisect
import random

import seqtree
from seqtree.gapblock import IslandProfile, central_prior, score_matrix

AA = "ACDEFGHIKLMNPQRSTVWY"
HUB = "CASSLGQAYEQYF"
CORE = (4, 9)               # the untemplated core; CASS and EQYF are the conserved flanks


def grow_island(hub, mat, rng, n=60):
    """An island: mutate the core, leave the flanks, and let one member lose a residue.

    Real islands look like this. Of the calibrated VDJdb islands with >= 10 members, 83% span
    exactly one length, and the variation sits away from the anchors.
    """
    cheap = {a: [b for b in AA if b != a and mat.penalty(a, b) <= 8] for a in AA}
    out = {hub}
    while len(out) < n:
        s = list(rng.choice(sorted(out)))
        for pos in rng.sample(range(*CORE), rng.choice([1, 1, 2])):
            s[pos] = rng.choice(cheap[s[pos]] or AA)
        out.add("".join(s))
    island = sorted(out)
    island[-1] = island[-1][:6] + island[-1][7:]       # one member is a residue short
    return island


def recall_at_fpr(pos, neg_sorted, fpr):
    """Recall at the largest cutoff admitting AT MOST ``fpr`` of the negatives.

    Not ``neg_sorted[k]``. Gap-block scores are small integers and tie heavily, so the k-th
    smallest can let hundreds of negatives through -- which would quietly hand min-over-members a
    looser effective FPR than the profile. This conservative rule is the one thetas_from_scores
    applies, for the same reason.
    """
    budget = int(fpr * len(neg_sorted))
    cut = None
    for v in sorted(set(neg_sorted)):
        if bisect.bisect_right(neg_sorted, v) > budget:
            break
        cut = v
    if cut is None:
        return sum(1 for x in pos if x < neg_sorted[0]) / len(pos)
    return sum(1 for x in pos if x <= cut) / len(pos)


def main():
    rng = random.Random(4)
    mat = seqtree.SubstitutionMatrix.blosum62()
    gap_open, lam = 2 * mat.scale(), int(1.5 * mat.scale())
    prior = central_prior(lam)

    control = seqtree.load_control("human_trb_aa")
    negatives = [control.ref_seq(i) for i in range(len(control))]
    M = len(negatives)

    island = grow_island(HUB, mat, rng)
    rng.shuffle(island)
    test, train = island[:20], island[20:]
    print(f"island of {len(island)} junctions around {HUB}: "
          f"{len(train)} for training, {len(test)} held out")
    print(f"negatives: {M:,} control junctions\n")

    profile = IslandProfile.fit(train)
    print(f"fitted {profile}")
    print(f"  consensus {profile.consensus()},  score(consensus) = "
          f"{profile.score(profile.consensus())}")
    print("  the frame column c minimises column entropy -- the members' own preference, not a")
    print("  constant we imposed. Across real islands its mode is 6, where the crystal structures")
    print("  put the block.\n")

    # ---- two scorers over the same negatives -------------------------------------------------
    # min-over-members is a row minimum of the (negatives x members) gap-block matrix. That matrix
    # is 250,000 x 40; scoring it in Python would take minutes, and in C++ it takes milliseconds.
    sm = score_matrix(negatives, train, mat, gap_open=gap_open, gap_prior=prior)
    mn_neg = sorted(min(sm.row(i)) for i in range(len(negatives)))
    st = score_matrix(test, train, mat, gap_open=gap_open, gap_prior=prior)
    mn_pos = [min(st.row(i)) for i in range(len(test))]

    pf_neg = sorted(profile.score_batch(negatives))
    pf_pos = profile.score_batch(test)

    # ---- where do the calibrated cutoffs actually sit? ---------------------------------------
    e_target = 0.05
    regimes = []
    for label, N in (("per-epitope island", 90), ("repertoire annotation", 20_000)):
        k = int(e_target * M / N)
        if k < 1:                       # rule of three: no E below 3N/M can be certified
            e_min = 3 * N / M
            k = int(e_min * M / N)
            note = f"E*={e_target} unreachable (3N/M={e_min:.3f}); using E*={e_min:.3f}"
        else:
            note = f"E*={e_target}"
        regimes.append((label, N, k, k / M, note))

    print("=== the E-value fixes the cutoff; the cutoff decides the winner ===")
    print(f"  {'regime':<24}{'N':>8}{'k':>6}{'FPR':>11}   note")
    for label, N, k, fpr, note in regimes:
        print(f"  {label:<24}{N:>8,}{k:>6}{100*fpr:>10.4f}%   {note}")

    print(f"\n=== recall on the {len(test)} held-out members, at a matched false-positive rate ===")
    print(f"  {'FPR':>11}{'min-over-members':>20}{'IslandProfile':>16}   regime")
    rows = [(1e-2, "loose reference (nobody uses this)")]
    rows += [(fpr, label) for label, _, _, fpr, _ in regimes]
    for fpr, label in sorted(rows, reverse=True):
        m = recall_at_fpr(mn_pos, mn_neg, fpr)
        p = recall_at_fpr(pf_pos, pf_neg, fpr)
        print(f"  {100*fpr:>10.4f}%{100*m:>19.1f}%{100*p:>15.1f}%   {label}")

    print(f"\n  Twenty held-out members is far too few to call a winner from -- one sequence moves")
    print("  a cell by 5 points. The direction is what to take away, and the measured verdict over")
    print("  108 real islands is: no significant difference at the per-epitope cutoff")
    print("  (+0.93 [-0.80, +2.79]), and +10.90 [+7.69, +14.21] at the repertoire cutoff.")

    print("\n  A profile is not a compression and does not generalise. It is 21 penalties per")
    print("  column against one string per member -- larger below 84 members -- and same-epitope")
    print("  junctions in a *different* island are recovered by neither representation. Fit one to")
    print("  recognise this island's members more sharply, not to discover new ones.")


if __name__ == "__main__":
    main()
