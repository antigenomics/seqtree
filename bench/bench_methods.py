#!/usr/bin/env python3
"""Compare the two engines across sizes, scopes, and budgets, and measure the
cost of fetching alignments.

  python bench/bench_methods.py                 # fast tier
  env RUN_BENCHMARK=1 python bench/bench_methods.py --sizes 100000 1000000

Sequences come from OLGA if available (realistic TCR CDR3), else seeded random.
"""
import argparse
import os
import random
import time

import seqtree

AA = "ACDEFGHIKLMNPQRSTVWY"


def gen_pool(n, rng):
    """Realistic CDR3 aa from OLGA when available, else seeded random."""
    try:
        import shutil
        import subprocess
        import tempfile

        if shutil.which("olga-generate_sequences"):
            with tempfile.NamedTemporaryFile("r", suffix=".tsv") as fh:
                subprocess.run(
                    ["olga-generate_sequences", "--humanTRB", "-n", str(min(n, 50000)),
                     "-o", fh.name],
                    check=True, capture_output=True, timeout=600,
                )
                pool = [line.split("\t")[1].strip() for line in open(fh.name) if "\t" in line]
            pool = [s for s in pool if s and all(c in AA for c in s)]
            if pool:
                return pool
    except Exception:
        pass
    return ["".join(rng.choice(AA) for _ in range(rng.randint(12, 18))) for _ in range(min(n, 50000))]


def expand(pool, target, rng):
    refs = list(pool)
    while len(refs) < target:
        s = rng.choice(pool)
        j = rng.randrange(len(s))
        refs.append(s[:j] + rng.choice(AA) + s[j + 1:])
    return refs[:target]


def mutate(s, n_subs, n_indels, rng):
    s = list(s)
    for _ in range(n_subs):
        j = rng.randrange(len(s))
        s[j] = rng.choice(AA)
    for _ in range(n_indels):
        j = rng.randrange(len(s))
        if rng.random() < 0.5 and len(s) > 1:
            del s[j]
        else:
            s.insert(j, rng.choice(AA))
    return "".join(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="*")
    ap.add_argument("--queries", type=int, default=None)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    slow = os.environ.get("RUN_BENCHMARK")
    sizes = args.sizes or ([100_000, 1_000_000] if slow else [5_000, 50_000])
    nq = args.queries if args.queries is not None else (200_000 if slow else 5_000)

    pool = gen_pool(max(sizes), rng)
    print(f"# pool of {len(pool)} base sequences; threads={args.threads or 'all'}", flush=True)

    print("phase\tn_refs\tn_queries\tengine\tknob\tvalue\tqps\thits/q")
    for n in sizes:
        refs = expand(pool, n, rng)
        idx = seqtree.Index.build(refs, alphabet="aa")
        base = [refs[rng.randrange(n)] for _ in range(min(nq, n))]

        # ---- scope sweep: edit-count budget, both engines ----
        for scope in (1, 2, 3):
            queries = [mutate(s, scope, 0, rng) for s in base]
            for eng in ("seqtm", "seqtrie"):
                p = seqtree.SearchParams(max_subs=scope, max_total_edits=scope, engine=eng)
                t0 = time.perf_counter()
                res = idx.search_batch(queries, p, threads=args.threads)
                dt = time.perf_counter() - t0
                hpq = sum(len(r) for r in res) / len(queries)
                print(f"scope\t{n}\t{len(queries)}\t{eng}\tmax_edits\t{scope}\t"
                      f"{len(queries) / dt:,.0f}\t{hpq:.2f}", flush=True)

        # ---- indel scope (seqtm only; seqtrie via total budget) ----
        queries = [mutate(s, 1, 1, rng) for s in base]
        for eng, p in (("seqtm", seqtree.SearchParams(max_subs=1, max_ins=1, max_dels=1, engine="seqtm")),
                       ("seqtrie", seqtree.SearchParams(max_total_edits=2, engine="seqtrie"))):
            t0 = time.perf_counter()
            res = idx.search_batch(queries, p, threads=args.threads)
            dt = time.perf_counter() - t0
            print(f"indel\t{n}\t{len(queries)}\t{eng}\tsub+indel\t1+1\t"
                  f"{len(queries) / dt:,.0f}\t{sum(len(r) for r in res) / len(queries):.2f}", flush=True)

        # ---- BLOSUM62 score budget (seqtrie) ----
        queries = [mutate(s, 2, 0, rng) for s in base]
        for budget in (6, 12, 20):
            p = seqtree.SearchParams(matrix="BLOSUM62", max_penalty=budget, engine="seqtrie", gap_open=8)
            t0 = time.perf_counter()
            res = idx.search_batch(queries, p, threads=args.threads)
            dt = time.perf_counter() - t0
            print(f"blosum\t{n}\t{len(queries)}\tseqtrie\tmax_penalty\t{budget}\t"
                  f"{len(queries) / dt:,.0f}\t{sum(len(r) for r in res) / len(queries):.2f}", flush=True)

        # ---- alignment fetch cost ----
        queries = [mutate(s, 1, 0, rng) for s in base[: min(5000, len(base))]]
        p = seqtree.SearchParams(max_subs=1, engine="seqtm")
        res = idx.search_batch(queries, p, threads=args.threads)
        pairs = [(h.ref_id, q) for q, hits in zip(queries, res) for h in hits]
        t0 = time.perf_counter()
        for ref_id, q in pairs:
            idx.align(ref_id, q, p)
        dt = time.perf_counter() - t0
        if pairs:
            print(f"align\t{n}\t{len(pairs)}\t-\talign_calls\t{len(pairs)}\t"
                  f"{len(pairs) / dt:,.0f}\t-", flush=True)


if __name__ == "__main__":
    main()
