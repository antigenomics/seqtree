#!/usr/bin/env python
"""Throughput of the batched gap-block kernel, at the shape a prototype embedding asks for.

A distance-vector embedding (TCREMP-style) scores every clonotype against a few thousand fixed
prototypes. That is a dense N x K matrix, not a bounded search: no reference can be pruned,
because the embedding needs the distance to all of them. `GapBlockIndex.search` is the wrong
tool -- `gapblock.score_matrix` is the right one.

Three rungs, ordered by the work each pair does:

    scalar     seqtree.gapblock.gapblock_score      pure Python, all L+1 block positions
    matrix-1   seqtree.gapblock.score_matrix        C++, all L+1, one thread
    matrix-N   seqtree.gapblock.score_matrix        C++, all L+1, one thread per core

For reference, mir.distances.aligner (mirpy's C extension, which considers exactly 4 fixed
block positions) measured 44.4 M pairs/s single-threaded on this machine, and BioPython's
PairwiseAligner 0.24 M. Those are not run here -- seqtree has no dependency on either.

Run:
    python bench/bench_score_matrix.py                # 3000 refs, up to 20k queries
    RUN_BENCHMARK=1 python bench/bench_score_matrix.py --queries 100000
"""
import argparse
import os
import random
import resource
import time

import seqtree
from seqtree.gapblock import central_prior, gapblock_score, positions_prior, score_matrix, _pen_table

AA = "ACDEFGHIKLMNPQRSTVWY"


def peak_rss_gb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / 1e9 if os.uname().sysname == "Darwin" else rss / 1e6


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--refs", type=int, default=3000, help="prototype references (matrix columns)")
    ap.add_argument("--queries", type=int, nargs="+", default=[400, 5000, 20000])
    ap.add_argument("--scalar-cap", type=int, default=400, help="skip pure Python above this")
    args = ap.parse_args()

    rng = random.Random(0)
    control = seqtree.load_control("human_trb_aa")
    pool = [control.ref_seq(i) for i in range(len(control))]
    refs = rng.sample(pool, args.refs)

    mat = seqtree.SubstitutionMatrix.blosum62()
    gap_open, lam = 2 * mat.scale(), int(1.5 * mat.scale())
    pen = _pen_table(mat, AA)
    cores = os.cpu_count() or 1

    print(f"{args.refs:,} prototype refs, human TRB control, BLOSUM62, gap_open={gap_open}, "
          f"{cores} cores\n")
    print(f"  {'rung':<26}{'queries':>9}{'pairs':>14}{'wall s':>9}{'M pairs/s':>12}{'speedup':>9}")

    base = {}
    for nq in args.queries:
        q = rng.sample(pool, nq)
        n = nq * len(refs)

        if nq <= args.scalar_cap:
            t = time.perf_counter()
            for a in q:
                for b in refs:
                    gapblock_score(a, b, gap_open=gap_open, gap_extend=1, _pen=pen)
            dt = time.perf_counter() - t
            base[nq] = n / dt
            print(f"  {'scalar (python)':<26}{nq:>9}{n:>14,}{dt:>9.3f}{n/dt/1e6:>12.2f}{1.0:>9.1f}x")

        for label, threads, prior in (
            ("matrix-1 (C++)", 1, None),
            ("matrix-1 + central prior", 1, central_prior(lam)),
            (f"matrix-{cores} + central prior", 0, central_prior(lam)),
            (f"matrix-{cores} + positions", 0, positions_prior((3, 4, -4, -3))),
        ):
            t = time.perf_counter()
            score_matrix(q, refs, mat, gap_open=gap_open, gap_prior=prior, threads=threads)
            dt = time.perf_counter() - t
            ref_rate = base.get(nq) or base.get(args.scalar_cap)
            speed = f"{(n/dt)/ref_rate:>8.0f}x" if ref_rate else f"{'':>9}"
            print(f"  {label:<26}{nq:>9}{n:>14,}{dt:>9.3f}{n/dt/1e6:>12.2f}{speed}")
        print()

    print(f"peak RSS {peak_rss_gb():.2f} GB")
    print("\n  The prior is free: it is materialised once into a [m][d][i] lookup cube, so the")
    print("  kernel pays one array read per candidate position and never re-enters Python.")
    print("  The result is int32 and travels by buffer protocol -- numpy.asarray does not copy.")
    print(f"  A {args.refs:,}-column matrix costs 4 * N * {args.refs:,} bytes; chunk N if that")
    print("  does not fit.")


if __name__ == "__main__":
    main()
