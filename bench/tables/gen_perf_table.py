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
from pathlib import Path

import seqtree as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # bench/, for _common

from _common import peak_rss_mb
from gen_retrieval_table import lcg_pool  # shared deterministic pool

N_REFS, N_QUERIES, LENGTH, SEED = 50_000, 5_000, 14, 1

# TextIndex reads the same pool as a concatenated text (~700 k residues). Two search rows, not
# one: L=14 takes an all-exact scheme (b = max_subs + 1) and L=9 takes a multi-block scheme
# carrying a per-block error budget. Only the second moves if the dispatch regresses, which is
# the whole reason it is gated -- restoring the pre-1.0 two-path dispatch made it 23x slower.
TEXT_K, TEXT_SHORT, TEXT_LONG, TEXT_SUBS = 4, 9, 14, 2


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

    text_build_ms = text_short_ms = text_long_ms = float("inf")
    for _ in range(args.repeats):
        t0 = time.perf_counter()
        tix = st.TextIndex.build(refs, alphabet="aa", k=TEXT_K)
        text_build_ms = min(text_build_ms, (time.perf_counter() - t0) * 1000)
        for length, slot in ((TEXT_SHORT, "short"), (TEXT_LONG, "long")):
            tq = [refs[(i * 7919) % N_REFS][:length] for i in range(N_QUERIES)]
            t0 = time.perf_counter()
            tix.search_batch(tq, max_subs=TEXT_SUBS, threads=1)
            dt = (time.perf_counter() - t0) * 1000
            if slot == "short":
                text_short_ms = min(text_short_ms, dt)
            else:
                text_long_ms = min(text_long_ms, dt)

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    out.write(f"# perf: n_refs={N_REFS} n_queries={N_QUERIES} length={LENGTH} "
              f"scope=2subs threads=1 best_of={args.repeats} "
              f"text_k={TEXT_K} text_L={TEXT_SHORT},{TEXT_LONG}\n")
    out.write("metric\tvalue\n")
    out.write(f"build_ms\t{build_ms:.1f}\n")
    out.write(f"search_ms\t{search_ms:.1f}\n")
    out.write(f"text_build_ms\t{text_build_ms:.1f}\n")
    out.write(f"text_search_short_ms\t{text_short_ms:.1f}\n")
    out.write(f"text_search_long_ms\t{text_long_ms:.1f}\n")
    out.write(f"peak_rss_mb\t{peak_rss_mb():.1f}\n")
    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
