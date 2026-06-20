#!/usr/bin/env python3
"""Produce the performance table (build/search time + peak RSS) on a fixed workload.

Runs the compiled seqtree from the current repo on a seeded synthetic reference set
and writes ``metric<TAB>value`` rows. Unlike the retrieval table this is *not* an
exact oracle — timings and memory vary by machine — so the regression test compares
against ``perf_baseline.tsv`` within a tolerance, not byte-for-byte.

  python bench/tables/gen_perf_table.py --out bench/tables/perf_baseline.tsv
"""
import argparse
import sys
import time

import seqtree as st

from gen_retrieval_table import lcg_pool  # shared deterministic pool

N_REFS, N_QUERIES, LENGTH, SEED = 50_000, 5_000, 14, 1


def peak_rss_mb():
    try:
        import resource

        peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        # ru_maxrss is bytes on macOS, kibibytes on Linux.
        return peak / (1024 * 1024) if sys.platform == "darwin" else peak / 1024
    except Exception:
        return 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="-")
    ap.add_argument("--repeats", type=int, default=3, help="take the best (min) of N runs")
    args = ap.parse_args()

    refs = lcg_pool(N_REFS, LENGTH, SEED)
    queries = [refs[(i * 7919) % N_REFS] for i in range(N_QUERIES)]
    p = st.SearchParams(max_subs=2, max_total_edits=2, engine="seqtm")

    build_ms = search_ms = float("inf")
    for _ in range(args.repeats):
        t0 = time.perf_counter()
        idx = st.Index.build(refs, alphabet="aa")
        build_ms = min(build_ms, (time.perf_counter() - t0) * 1000)
        t0 = time.perf_counter()
        idx.search_batch(queries, p, threads=1)
        search_ms = min(search_ms, (time.perf_counter() - t0) * 1000)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    out.write(f"# perf: n_refs={N_REFS} n_queries={N_QUERIES} length={LENGTH} "
              f"scope=2subs threads=1 best_of={args.repeats}\n")
    out.write("metric\tvalue\n")
    out.write(f"build_ms\t{build_ms:.1f}\n")
    out.write(f"search_ms\t{search_ms:.1f}\n")
    out.write(f"peak_rss_mb\t{peak_rss_mb():.1f}\n")
    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
