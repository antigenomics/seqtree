#!/usr/bin/env python
"""Needleman-Wunsch / Smith-Waterman throughput, against BioPython.

seqtree.pairwise exists so that ordinary protein alignment does not need BioPython. It is a
drop-in: `tests/python/test_pairwise.py` runs 1,800 comparisons per (matrix, mode, gap) cell
against Bio.Align.PairwiseAligner and requires exact agreement. This script asks the other
question -- how much faster.

The shape that matters is the one a prototype embedding or a germline distance table has: a few
hundred to a few thousand sequences, all-against-all, where BioPython's per-pair Python call
overhead dominates.

BioPython is optional here; skip it with --no-oracle.

Run:
    python bench/bench_pairwise.py
    python bench/bench_pairwise.py --sizes 100 300 --lengths 90
"""
import argparse
import random
import time

import seqtree
from seqtree.pairwise import dist_matrix, score, score_matrix

AA = "ACDEFGHIKLMNPQRSTVWY"


def seqs(rng, n, length):
    return ["".join(rng.choice(AA) for _ in range(rng.randint(length - 8, length + 8)))
            for _ in range(n)]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="+", default=[100, 200, 400])
    ap.add_argument("--lengths", type=int, nargs="+", default=[15, 90],
                    help="15 ~ a CDR3 junction; 90 ~ a germline V gene")
    ap.add_argument("--no-oracle", action="store_true", help="skip the BioPython comparison")
    args = ap.parse_args()

    rng = random.Random(0)
    mat = seqtree.SubstitutionMatrix.blosum62()
    go, ge = 11, 1

    oracle = None
    if not args.no_oracle:
        try:
            from Bio import Align

            oracle = Align.PairwiseAligner()
            oracle.mode = "global"
            from Bio.Align import substitution_matrices

            oracle.substitution_matrix = substitution_matrices.load("BLOSUM62")
            oracle.open_gap_score = -go
            oracle.extend_gap_score = -ge
        except ImportError:
            print("BioPython not installed; running seqtree only.\n")

    print(f"BLOSUM62, global, gap_open={go} gap_extend={ge}\n")
    header = f"  {'length':>7}{'n':>6}{'pairs':>10}{'seqtree 1t':>13}{'seqtree Nt':>13}"
    if oracle:
        header += f"{'biopython':>13}{'speedup':>9}"
    print(header)

    for length in args.lengths:
        for n in args.sizes:
            xs = seqs(rng, n, length)
            pairs = n * n

            t = time.perf_counter()
            score_matrix(xs, xs, mat, gap_open=go, gap_extend=ge, threads=1)
            t1 = time.perf_counter() - t

            t = time.perf_counter()
            score_matrix(xs, xs, mat, gap_open=go, gap_extend=ge, threads=0)
            tn = time.perf_counter() - t

            row = (f"  {length:>7}{n:>6}{pairs:>10,}"
                   f"{pairs/t1/1e3:>10,.0f} k/s{pairs/tn/1e3:>10,.0f} k/s")

            if oracle:
                # BioPython on a subset, extrapolated -- the full matrix would take minutes
                k = min(n, 20)
                t = time.perf_counter()
                for a in xs[:k]:
                    for b in xs[:k]:
                        oracle.align(a, b).score
                tb = (time.perf_counter() - t) * (pairs / (k * k))
                row += f"{pairs/tb/1e3:>10,.1f} k/s{tb/tn:>8.0f}x"
            print(row)

    print("\n  seqtree's per-pair cost is a C++ Gotoh with O(min(m,n)) memory and no Python in the")
    print("  loop; BioPython pays a Python call and an Alignment object per pair. The N-thread")
    print("  column releases the GIL, so it scales with cores.")

    # the distance matrix is what a prototype embedding actually consumes
    xs = seqs(rng, 200, 90)
    t = time.perf_counter()
    d = dist_matrix(xs, xs, mat, gap_open=go, gap_extend=ge, threads=0)
    dt = time.perf_counter() - t
    print(f"\n  dist_matrix (200 x 200 germline-length): {dt*1000:.0f} ms")
    print(f"    d = s(a,a) + s(b,b) - 2*s(a,b); self-scores taken once per sequence, not per pair.")
    print(f"    diagonal is zero: {all(d[i, i] == 0 for i in range(len(xs)))}")


if __name__ == "__main__":
    main()
