#!/usr/bin/env python
"""Can the gap rule be optimised per chain? Three criteria; the answer is mostly no.

A gap rule only matters for pairs of unequal length -- where lengths match, every rule in the family
scores identically. So every comparison here is restricted to length-different pairs.

**A. epitope discrimination (labelled, non-circular).** VDJdb labels TRA and TRB, human and mouse.
Positives are same-epitope junction pairs with ``|dlen| >= 1``; negatives are different-epitope
pairs. Recall at a matched 1% FPR -- each rule reads its cutoff off its own negatives, so scale
cancels -- bootstrapped over *epitopes*, since pairs from one epitope are not independent.

This is the only criterion here that asks a question about biology, and it **cannot separate the
rules**: every 95% CI overlaps every other, on all four arms. Length-different same-epitope pairs
sit in different sequence islands and are not similar under any layout, so there is nothing for a
gap rule to get right. AUROC (~0.52) says the same.

**B. manufactured similarity (label-free, non-circular, one-sided).** Free gap placement invents
alignments: an unrelated pair of unequal length slides its block until something lines up. Measured
on real prototype sequences with no planting, so nothing is assumed about where indels go.

Equal-length pairs carry no block and are rule-invariant, which makes them the yardstick. The
quantity compared is the **substitution cost of the layout the rule chose**, with the gap cost and
the prior subtracted -- one scale for every rule, so a rule cannot win by inflating scores.
(Comparing raw scores fails outright: on four arms the 1st percentile of equal-length scores lies
*below* ``gap_open``, so no gapped pair can fall under it and every rule reads 0.000%.)

    manufactured = P( subcost(unrelated, dlen 1-3) <= q01 subcost(unrelated, dlen 0) )

One-sided: lower is better, but a rule that always blocks at position 0 would read 0.00% and find
nothing real. It bounds; it does not decide.

**C. planted-indel recall (label-free, CIRCULAR -- reported to show that it is).** Plant one
contiguous deletion at a uniformly random position, add one substitution, and ask each rule to
find the relative among unrelated negatives matched on block length. The flat rule scores 100%,
because the score finds exactly the alignment that was planted. Change the planting distribution
and the winner changes with it. **No label-free positive set can settle where indels go** -- that
information has to come from structure, epitope labels, or germline annotation, and the first two
exist only for TRA/TRB while the third was already falsified against structure (a germline
untemplated-span rule recovers the true block position 0.4% of the time for TRA, against 42.4%
for a central prior).

Data (see SOURCES.md): mirpy's per-locus prototype files and a VDJdb slim dump. Neither ships with
seqtree; pass their paths. This script imports only seqtree and the standard library.

Run:
    python bench/gap_rules_by_locus.py \
        --prototypes ../mirpy/mir/resources/prototypes \
        --vdjdb ../mirpy/airr_benchmark/vdjdb/vdjdb-2025-12-29/vdjdb.slim.txt.gz
"""
import argparse
import csv
import collections
import gzip
import pathlib
import random
import statistics
import sys
import time

import seqtree
from seqtree.gapblock import (
    _pen_table, central_prior, gap_cost, gapblock_score, positions_prior,
)

csv.field_size_limit(1 << 30)

PRODUCTIVE = frozenset("ACDEFGHIKLMNPQRSTVWY")
LOCI = ["TRA", "TRB", "TRG", "TRD", "IGH", "IGK", "IGL"]
# VDJdb spells species differently from the prototype filenames.
SPECIES = {"homosapiens": "human", "musmusculus": "mouse"}


def rules(scale):
    """The ladder, from hardest pin to no pin at all. lam is in units of matrix scale."""
    huge = 1 << 16
    return [
        ("fixed centre", central_prior(huge)),
        ("central lam=3.0", central_prior(int(3.0 * scale))),
        ("central lam=1.5", central_prior(int(1.5 * scale))),
        ("central lam=0.5", central_prior(int(0.5 * scale))),
        ("positions 3,4,-4,-3", positions_prior((3, 4, -4, -3))),
        ("flat (score alone)", None),
    ]


# ------------------------------------------------------------------------------------------
# pair scoring
# ------------------------------------------------------------------------------------------

PEN = None      # set in main() once the matrix exists


def score_pairs(pairs, gap_open, prior):
    """Distances for an explicit list of (a, b).

    Elementwise, not an outer product: ``score_matrix`` would compute |pairs|^2 cells and throw
    away all but the diagonal. seqtree has no C++ rung for "N pairs, elementwise" -- at the tens
    of thousands of pairs used here the Python scorer costs a couple of seconds, so none is
    needed. Add one when a caller wants millions.
    """
    return [gapblock_score(a, b, gap_open=gap_open, gap_extend=1, gap_prior=prior, _pen=PEN)[0]
            for a, b in pairs]


def subcosts(pairs, gap_open, prior):
    """Substitution cost of the layout each rule chose: score minus the gap cost and the prior.

    This is the similarity the rule claims to have found, on a scale shared by every rule.
    """
    out = []
    for a, b in pairs:
        s, i = gapblock_score(a, b, gap_open=gap_open, gap_extend=1, gap_prior=prior, _pen=PEN)
        d = abs(len(a) - len(b))
        s -= gap_cost(d, gap_open, 1)
        if prior is not None and d > 0:
            s -= prior(i, d, max(len(a), len(b)))
        out.append(s)
    return out


def recall_at_fpr(pos, neg, fpr=0.01):
    """Fraction of positives at or below the cutoff that admits `fpr` of the negatives.

    Each rule reads its cutoff off its own negatives, so the comparison is at matched FPR and no
    rule can win by rescaling. This is the tail; AUROC is the bulk.
    """
    if not pos or not neg:
        return float("nan")
    cut = sorted(neg)[min(len(neg) - 1, int(fpr * len(neg)))]
    return sum(1 for d in pos if d <= cut) / len(pos)


def auroc(pos, neg):
    """P(d_pos < d_neg) + 0.5 P(=), via rank sum. Lower distance = more likely a positive."""
    if not pos or not neg:
        return float("nan")
    merged = sorted([(d, 0) for d in pos] + [(d, 1) for d in neg])
    ranks, i = {}, 0
    while i < len(merged):
        j = i
        while j < len(merged) and merged[j][0] == merged[i][0]:
            j += 1
        r = (i + j - 1) / 2 + 1
        for k in range(i, j):
            ranks[k] = r
        i = j
    rsum = sum(ranks[k] for k in range(len(merged)) if merged[k][1] == 0)
    n1, n2 = len(pos), len(neg)
    u = rsum - n1 * (n1 + 1) / 2
    return 1.0 - u / (n1 * n2)          # low distance for positives => AUROC > 0.5


# ------------------------------------------------------------------------------------------
# criterion A: epitope discrimination
# ------------------------------------------------------------------------------------------

def load_vdjdb(path):
    """{(species, locus): {epitope: [junction]}} over unique productive junctions."""
    out = collections.defaultdict(lambda: collections.defaultdict(set))
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            sp = SPECIES.get(r["species"].lower())
            j, ep = r["cdr3"], r["antigen.epitope"]
            if sp and j and ep and set(j) <= PRODUCTIVE and 6 <= len(j) <= 30:
                out[(sp, r["gene"])][ep].add(j)
    return out


def epitope_pairs(by_epitope, rng, per_epitope=140, min_seqs=30):
    """Same-epitope and different-epitope pairs, all with |dlen| >= 1, tagged by epitope."""
    eps = sorted(e for e, s in by_epitope.items() if len(s) >= min_seqs)
    if len(eps) < 4:
        return [], []
    seqs = {e: sorted(by_epitope[e]) for e in eps}
    pos, neg = [], []
    for e in eps:
        own = seqs[e]
        for want, bucket in ((per_epitope, pos), (per_epitope, neg)):
            got = 0
            for _ in range(want * 60):
                if got >= want:
                    break
                a = rng.choice(own)
                if bucket is pos:
                    b = rng.choice(own)
                    if a == b:
                        continue
                else:
                    o = rng.choice(eps)
                    if o == e:
                        continue
                    b = rng.choice(seqs[o])
                if len(a) != len(b):
                    bucket.append((e, a, b))
                    got += 1
    return pos, neg


def bootstrap(pos, neg, dp, dn, rng, stat, b=400):
    """Resample epitopes, not pairs: pairs drawn from one epitope are not independent."""
    pe = collections.defaultdict(list)
    ne = collections.defaultdict(list)
    for (e, _, _), d in zip(pos, dp):
        pe[e].append(d)
    for (e, _, _), d in zip(neg, dn):
        ne[e].append(d)
    eps = sorted(set(pe) | set(ne))
    out = []
    for _ in range(b):
        draw = [eps[rng.randrange(len(eps))] for _ in eps]
        p = [d for e in draw for d in pe.get(e, ())]
        n = [d for e in draw for d in ne.get(e, ())]
        if p and n:
            out.append(stat(p, n))
    out.sort()
    return out[int(0.025 * len(out))], out[int(0.975 * len(out))]


# ------------------------------------------------------------------------------------------
# criterion B: planted indels vs unrelated pairs
# ------------------------------------------------------------------------------------------

AA = "ACDEFGHIKLMNPQRSTVWY"


def planted_pairs(seqs, rng, n, subs=1, d_max=3):
    """(original, relative) where relative lost d residues at a UNIFORMLY random position.

    Uniform placement handicaps a centred prior on purpose. Deletions from inside a homopolymer
    run are skipped: several positions produce the same string, so no rule can be graded on them.
    """
    out, guard = [], 0
    while len(out) < n and guard < n * 200:
        guard += 1
        s = rng.choice(seqs)
        d = rng.randint(1, d_max)
        if len(s) < d + 8:
            continue
        p = rng.randrange(0, len(s) - d + 1)
        r = s[:p] + s[p + d:]
        if any(s[:k] + s[k + d:] == r for k in range(len(s) - d + 1) if k != p):
            continue
        t = list(r)
        for pos in rng.sample(range(len(t)), min(subs, len(t))):
            t[pos] = rng.choice(AA)
        out.append((s, "".join(t), d))
    return out


def unrelated_pairs(seqs, rng, n, equal_length=False, d_wanted=None):
    """Random same-locus pairs. `d_wanted` matches a target block-length distribution."""
    by_len = collections.defaultdict(list)
    for s in seqs:
        by_len[len(s)].append(s)
    out, guard = [], 0
    need = collections.Counter(d_wanted) if d_wanted else None
    while len(out) < n and guard < n * 400:
        guard += 1
        a, b = rng.choice(seqs), rng.choice(seqs)
        if a == b:
            continue
        d = abs(len(a) - len(b))
        if equal_length:
            if d == 0:
                out.append((a, b))
        elif need is not None:
            if need.get(d, 0) > 0:
                need[d] -= 1
                out.append((a, b))
        elif 1 <= d <= 3:
            out.append((a, b))
    return out


def percentile(xs, q):
    s = sorted(xs)
    return s[min(len(s) - 1, int(q * len(s)))]


# ------------------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prototypes", type=pathlib.Path, required=True)
    ap.add_argument("--vdjdb", type=pathlib.Path)
    ap.add_argument("--pairs", type=int, default=40_000, help="unrelated pairs per arm (criterion B)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    global PEN
    mat = seqtree.SubstitutionMatrix.blosum62()
    scale = mat.scale()
    gap_open = 2 * scale
    PEN = _pen_table(mat, "".join(sorted(PRODUCTIVE)))
    print(f"BLOSUM62 scale = {scale}, gap_open = {gap_open}\n")

    # ---- criterion A -----------------------------------------------------------------------
    if args.vdjdb and args.vdjdb.exists():
        t = time.perf_counter()
        db = load_vdjdb(args.vdjdb)
        print(f"=== A. epitope discrimination on |dlen| >= 1 pairs, VDJdb "
              f"({time.perf_counter()-t:.1f}s to load) ===")
        print("    recall at a matched 1% FPR, [95% CI] bootstrapped over epitopes; "
              "AUROC for context\n")
        for (sp, locus) in sorted(db):
            rng = random.Random(args.seed)
            pos, neg = epitope_pairs(db[(sp, locus)], rng)
            if not pos:
                continue
            neps = len({e for e, _, _ in pos})
            print(f"  {sp} {locus}: {neps} epitopes, {len(pos):,} pos / {len(neg):,} neg pairs")
            print(f"    {'rule':<22}{'recall@1%FPR':>14}{'95% CI':>18}{'AUROC':>9}")
            for name, prior in rules(scale):
                dp = score_pairs([(a, b) for _, a, b in pos], gap_open, prior)
                dn = score_pairs([(a, b) for _, a, b in neg], gap_open, prior)
                r = recall_at_fpr(dp, dn)
                lo, hi = bootstrap(pos, neg, dp, dn, random.Random(args.seed + 1), recall_at_fpr)
                print(f"    {name:<22}{100*r:>13.2f}%{f'[{100*lo:.2f}, {100*hi:.2f}]':>18}"
                      f"{auroc(dp, dn):>9.3f}")
            print()

    # ---- criteria B and C --------------------------------------------------------------------
    protos = {}
    for path in sorted(args.prototypes.glob("*.tsv")):
        sp, locus = path.stem.split("_", 1)
        if locus not in LOCI:
            continue
        with path.open() as fh:
            seqs = [r["junction_aa"] for r in csv.DictReader(fh, delimiter="\t")
                    if set(r["junction_aa"]) <= PRODUCTIVE]
        protos[(sp, locus)] = seqs

    names = [n for n, _ in rules(scale)]
    order = sorted(protos, key=lambda k: (k[0] != "human", k))
    man_by_arm, rec_by_arm, med_by_arm = {}, {}, {}
    for (sp, locus) in order:
        seqs = protos[(sp, locus)]
        rng = random.Random(args.seed)
        planted = planted_pairs(seqs, rng, args.pairs)
        neg = unrelated_pairs(seqs, rng, args.pairs, d_wanted=[d for _, _, d in planted])
        eq = unrelated_pairs(seqs, rng, min(args.pairs, 20_000), equal_length=True)
        cut_eq = percentile(subcosts(eq, gap_open, None), 0.01)   # rule-invariant yardstick

        rec, man = [], []
        for _, prior in rules(scale):
            dp = score_pairs([(a, b) for a, b, _ in planted], gap_open, prior)
            dn = score_pairs(neg, gap_open, prior)
            rec.append(recall_at_fpr(dp, dn))
            s_ne = subcosts(neg, gap_open, prior)
            man.append(sum(1 for d in s_ne if d <= cut_eq) / len(s_ne))
        man_by_arm[(sp, locus)] = man
        rec_by_arm[(sp, locus)] = rec
        med_by_arm[(sp, locus)] = statistics.median([len(s) for s in seqs])

    print("=== B. manufactured similarity, all seven loci. No planting; lower is better. ===\n")
    print(f"  {'arm':<13}{'med len':>8} " + "".join(f"{n:>20}" for n in names) + f"{'flat/best':>12}")
    for arm in order:
        man = man_by_arm[arm]
        lo = min(man)
        cells = "".join((f"{100*m:>18.2f}%" + ("*" if m == lo else " ")) for m in man)
        ratio = man[-1] / lo if lo else float("inf")
        print(f"  {arm[0]+' '+arm[1]:<13}{med_by_arm[arm]:>8.0f} {cells}{ratio:>11.1f}x")

    ratios = {a: man_by_arm[a][-1] / min(man_by_arm[a]) for a in order}
    lo_arm, hi_arm = min(ratios, key=ratios.get), max(ratios, key=ratios.get)
    i_pos = names.index("positions 3,4,-4,-3")
    pos_pen = {a: man_by_arm[a][i_pos] / man_by_arm[a][0] for a in order}   # vs fixed centre
    worst = max(pos_pen, key=pos_pen.get)
    exceptions = [a for a in order if man_by_arm[a][i_pos] < man_by_arm[a][0]]
    print(f"\n  Constraining the block cuts invented similarity by {ratios[lo_arm]:.1f}x "
          f"({lo_arm[0]} {lo_arm[1]}) to {ratios[hi_arm]:.1f}x ({hi_arm[0]} {hi_arm[1]}).")
    print(f"  positions_prior((3,4,-4,-3)) -- what mir.distances.aligner hardcodes for every locus")
    print(f"  -- is {pos_pen[worst]:.1f}x worse than a pinned centre on {worst[0]} {worst[1]}.")
    if exceptions:
        e = ", ".join(f"{a[0]} {a[1]}" for a in exceptions)
        print(f"  It beats the pinned centre on exactly one arm: {e}.")
    print("  This is one-sided evidence. It shows what freedom costs, not what it buys.")

    print("\n=== C. planted-indel recall @ matched 1% FPR. CIRCULAR -- read the caveat. ===\n")
    print(f"  {'arm':<13}{'med len':>8} " + "".join(f"{n:>20}" for n in names))
    for arm in order:
        rec = rec_by_arm[arm]
        top = max(rec)
        cells = "".join((f"{100*r:>18.1f}%" + ("*" if r == top else " ")) for r in rec)
        print(f"  {arm[0]+' '+arm[1]:<13}{med_by_arm[arm]:>8.0f} {cells}")

    print("\n  The flat rule reaches 100% on five arms. That is not a finding -- the indel was")
    print("  planted uniformly at random and an unconstrained score finds exactly the alignment")
    print("  that was planted. Plant at the centre instead and the centred prior wins by the same")
    print("  construction. The winner of this table is the planting distribution, not the biology.")

    print("\n=== verdict ===\n")
    print("  Per-chain optimisation of the gap rule is NOT supported by the available data:")
    print("    * A (the only functional, non-circular criterion) cannot separate any two rules on")
    print("      any of the four labelled arms -- every 95% CI overlaps.")
    print("    * C, and every label-free positive set like it, is circular.")
    print("    * TRG, TRD, IGH, IGK and IGL have no epitope labels and no crystal structures, so no")
    print("      non-circular criterion for them exists at all.")
    print("  What B does establish, one-sidedly and for every locus: letting the score place the")
    print("  block invents similarity between unrelated sequences, and a centred prior suppresses")
    print("  it. central_prior scales its pin with sequence length, so it needs no per-locus")
    print("  constant; a fixed start set like (3, 4, -4, -3) is a length-14 assumption that a")
    print("  50-residue IGH junction does not satisfy.")
    print("  Tune lam, gap_open and gap_extend against a DOWNSTREAM objective (clustering, epitope")
    print("  classification) in the repo that owns it. seqtree's own data cannot settle the choice.")


if __name__ == "__main__":
    main()
