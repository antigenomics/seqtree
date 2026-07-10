#!/usr/bin/env python
"""The gap-freedom ladder: what does constraining the gap actually cost?

Five rungs, ordered by how free the gap block is to move:

    (i)   fixed centre     central_prior(huge)  -- exactly one layout survives
    (ii)  central prior    central_prior(lam)   -- soft pull to the loop centre
    (iii) profile prior    profile_prior(lam, w)-- soft pull to a low-weight window
    (iv)  flat             gap_prior=None       -- any layout, chosen by score alone
    (v)   full affine      Index.align          -- any NUMBER of blocks, anywhere

lam -> infinity turns (ii) into (i); lam -> 0 turns it into (iv). Only (v) leaves the
single-gap-block family at all.

Two questions, both measured here on generic synthetic data (no domain assumptions):

  G1  how often, and by how much, does the single-block score exceed the unrestricted affine
      optimum? Gap-block is a strict restriction, so it can only ever be >=.
  perf  us/query, hits/query and peak RSS for the searchable rungs (i)-(iv) over a real index.

Rung (v) is not searchable -- ``Index.align`` scores one pair at a time and ``pairwise_batch``
is a *bounded trie search*, not affine alignment. It appears here only as the score oracle.

Usage:
    python bench/bench_gapblock.py                    # fast tier
    RUN_BENCHMARK=1 python bench/bench_gapblock.py    # adds the 250k-reference throughput table

2026-07-10
"""
from __future__ import annotations

import argparse
import gzip
import os
import random
import resource
import statistics as stt
import sys
import time
from importlib import resources

import seqtree
from seqtree.gapblock import GapBlockIndex, central_prior, gapblock_score, profile_prior

AA = "ACDEFGHIKLMNPQRSTVWY"
HUGE = 10 ** 6


def peak_rss_gb() -> float:
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return r / (1024 ** 3) if sys.platform == "darwin" else r / (1024 ** 2)


def u_shaped(j: int, m: int) -> float:
    """A generic 'framework at the ends, insert in the middle' weight, in [0, 1]."""
    return abs(2 * j - (m - 1)) / (m - 1) if m > 1 else 0.0


def rungs(lam: int):
    return [
        ("i   fixed centre", central_prior(HUGE)),
        ("ii  central prior", central_prior(lam)),
        ("iii profile prior", profile_prior(lam, u_shaped)),
        ("iv  flat (free)", None),
    ]


def synth(n, rng, lo=10, hi=18):
    return ["".join(rng.choice(AA) for _ in range(rng.randint(lo, hi))) for _ in range(n)]


def related_pairs(n, rng, d_choices=(0, 1, 2, 3, 4), subs=(0, 1, 2)):
    """``r`` is ``q`` with ONE contiguous block removed and a few substitutions. The single-block
    model is exactly the right model here, so any excess over affine is the price of the
    restriction on data that satisfies its assumption."""
    out = []
    for _ in range(n):
        q = "".join(rng.choice(AA) for _ in range(rng.randint(11, 18)))
        d = rng.choice(d_choices)
        i = rng.randrange(0, len(q) - d) if d else 0
        r = list(q[:i] + q[i + d:])
        for pos in rng.sample(range(len(r)), min(rng.choice(subs), len(r))):
            r[pos] = rng.choice(AA)
        out.append((q, "".join(r)))
    return out


def _affine(mat, q, r, go):
    idx = seqtree.Index.build([r], "aa")
    params = seqtree.SearchParams(max_subs=len(q) + len(r), max_penalty=HUGE,
                                  matrix=mat, gap_open=go, gap_extend=1, engine="seqtm")
    return idx.align(0, q, params).score


def g1(mat, related, unrelated, gap_opens):
    """Score cost of the single-block restriction, against Index.align's affine optimum.

    On *related* pairs the restriction should be nearly free. On *unrelated* pairs affine will
    always undercut it by inventing gaps -- which is not a win, it is manufactured similarity,
    and it is exactly why extra gap freedom hurts retrieval precision.
    """
    print("\n=== G1: does restricting to ONE gap block cost score? ===")
    print("  gap-block >= affine always (a strict restriction). Flat prior throughout, so the")
    print("  only difference between the two is how many blocks affine is allowed to open.\n")
    for label, pairs in (("RELATED (1 indel + 0-2 subs)", related),
                         ("UNRELATED (independent random)", unrelated)):
        print(f"  --- {label}, n = {len(pairs)} ---")
        print(f"  {'gap_open':>9}{'exact tie':>12}{'block > affine':>16}{'median excess':>15}"
              f"{'p90 excess':>12}")
        for go in gap_opens:
            ex = [gapblock_score(q, r, mat, go, 1, None)[0] - _affine(mat, q, r, go)
                  for q, r in pairs]
            assert min(ex) >= 0, "gap-block scored BELOW affine: it is not a restriction"
            worse = [e for e in ex if e > 0]
            med = stt.median(worse) if worse else 0
            p90 = sorted(worse)[int(0.9 * len(worse))] if worse else 0
            print(f"  {go:>9}{100*(1-len(worse)/len(ex)):>11.1f}%{100*len(worse)/len(ex):>15.1f}%"
                  f"{med:>15}{p90:>12}")
        print()


def layout_agreement(mat, pairs, lam, go):
    """Do the rungs even choose different layouts? If not, the ladder is moot."""
    print("\n=== do the rungs disagree about where the block goes? ===")
    print(f"  {'rung':<20}{'pairs':>8}{'agrees with flat':>20}{'agrees with centre':>21}")
    ref_flat, ref_cen = {}, {}
    for k, (q, r) in enumerate(pairs):
        if len(q) == len(r):
            continue
        ref_flat[k] = gapblock_score(q, r, mat, go, 1, None)[1]
        ref_cen[k] = gapblock_score(q, r, mat, go, 1, central_prior(HUGE))[1]
    for name, prior in rungs(lam):
        same_f = same_c = 0
        for k, (q, r) in enumerate(pairs):
            if k not in ref_flat:
                continue
            i = gapblock_score(q, r, mat, go, 1, prior)[1]
            same_f += i == ref_flat[k]
            same_c += i == ref_cen[k]
        n = len(ref_flat)
        print(f"  {name:<20}{n:>8}{100*same_f/n:>19.1f}%{100*same_c/n:>20.1f}%")


def throughput(mat, lam, budget, d_max, n_queries, seed):
    print(f"\n=== throughput: 250k references, budget {budget}, d_max {d_max} ===")
    with resources.files("seqtree").joinpath("data/control_human_trb_aa.txt.gz").open("rb") as fh:
        refs = [ln.strip() for ln in gzip.open(fh, "rt") if ln.strip()]
    rss0 = peak_rss_gb()
    t0 = time.perf_counter()
    gbi = GapBlockIndex(refs, "aa", d_max=d_max)
    build = time.perf_counter() - t0
    print(f"  build {build:.1f}s, peak RSS {peak_rss_gb():.2f} GB (+{peak_rss_gb()-rss0:.2f})")
    rng = random.Random(seed)
    qs = [refs[rng.randrange(len(refs))] for _ in range(n_queries)]
    go = 2 * mat.scale()
    print(f"\n  {'rung':<20}{'us/query':>11}{'hits/query':>13}")
    for name, prior in rungs(lam):
        t0 = time.perf_counter()
        hits = sum(len(gbi.search(q, budget, mat, gap_open=go, gap_prior=prior)) for q in qs)
        dt = 1e6 * (time.perf_counter() - t0) / len(qs)
        print(f"  {name:<20}{dt:>11.0f}{hits/len(qs):>13.1f}")
    print("\n  Hamming-only reference (d_max = 0, same budget):")
    t0 = time.perf_counter()
    base = seqtree.Index.build(refs, "aa")
    p = seqtree.SearchParams(max_subs=30, max_penalty=budget, matrix=mat, engine="seqtm")
    hits = sum(len(base.search(q, p)) for q in qs)
    print(f"  {'    Hamming ball':<20}{1e6*(time.perf_counter()-t0)/len(qs):>11.0f}"
          f"{hits/len(qs):>13.1f}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pairs", type=int, default=600)
    ap.add_argument("--budget", type=int, default=40)
    ap.add_argument("--d-max", type=int, default=2)
    ap.add_argument("--queries", type=int, default=200)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    mat = seqtree.SubstitutionMatrix.blosum62()
    lam = int(1.5 * mat.scale())
    print(f"BLOSUM62 Gram scale (median mismatch) = {mat.scale()}, lam = {lam}")

    rng = random.Random(args.seed)
    related = related_pairs(args.pairs, rng)
    pool = synth(2 * args.pairs, rng)
    unrelated = [(pool[2 * k], pool[2 * k + 1]) for k in range(args.pairs)]
    unrelated = [(q, r) for q, r in unrelated if abs(len(q) - len(r)) <= 4]

    t0 = time.perf_counter()
    g1(mat, related, unrelated, [1, 14, 28])
    layout_agreement(mat, related, lam, 2 * mat.scale())

    if os.getenv("RUN_BENCHMARK"):
        throughput(mat, lam, args.budget, args.d_max, args.queries, args.seed)
    else:
        print("\n  (set RUN_BENCHMARK=1 for the 250k-reference throughput table)")

    print(f"\nwall {time.perf_counter()-t0:.0f}s, peak RSS {peak_rss_gb():.2f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
