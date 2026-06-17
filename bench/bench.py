#!/usr/bin/env python3
"""Benchmark seqtree on real VDJdb data with ground-truth recall.

Downloads the AIRR VDJdb slim table, extracts CDR3 (aa) and epitope (aa) pools,
mutates known references to build queries with known answers, optionally
synthesizes up to 1M refs/queries, and reports recall, throughput, and peak RSS.

Fast tier (default):   python bench/bench.py
1M tier:               env RUN_BENCHMARK=1 python bench/bench.py --sizes 1000000 --queries 1000000

Needs: huggingface_hub, psutil (pip install -e ".[bench]").
"""
import argparse
import csv
import gzip
import os
import random
import time

import seqtree

csv.field_size_limit(10**7)  # slim table has wide metadata fields

REPO = "isalgo/airr_benchmark"
SLIM = "vdjdb/vdjdb-2025-12-29/vdjdb.slim.txt.gz"
CDR3_COLS = ("cdr3", "junction_aa", "cdr3_aa")
EPI_COLS = ("antigen.epitope", "epitope", "antigen_epitope")
AA = "ACDEFGHIKLMNPQRSTVWY"


def load_pools():
    from huggingface_hub import hf_hub_download

    path = hf_hub_download(repo_id=REPO, filename=SLIM, repo_type="dataset")
    cdr3, epi = set(), set()
    with gzip.open(path, "rt") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = reader.fieldnames or []
        cdr3_col = next((c for c in CDR3_COLS if c in cols), None)
        epi_col = next((c for c in EPI_COLS if c in cols), None)
        for row in reader:
            if cdr3_col:
                v = (row.get(cdr3_col) or "").strip().upper()
                if v and all(c in AA for c in v):
                    cdr3.add(v)
            if epi_col:
                v = (row.get(epi_col) or "").strip().upper()
                if v and all(c in AA for c in v):
                    epi.add(v)
    return sorted(cdr3), sorted(epi)


def mutate(seq, n_subs, rng):
    if not seq:
        return seq
    s = list(seq)
    for _ in range(n_subs):
        i = rng.randrange(len(s))
        s[i] = rng.choice(AA)
    return "".join(s)


def make_refs(pool, target, rng):
    refs = list(pool)
    while len(refs) < target:
        refs.append(mutate(rng.choice(pool), rng.randint(0, 3), rng))
    return refs[:target]


def make_queries(refs, n, n_subs, rng):
    """Return (queries, ground_truth_ref_ids)."""
    qs, gts = [], []
    for _ in range(n):
        rid = rng.randrange(len(refs))
        qs.append(mutate(refs[rid], n_subs, rng))
        gts.append(rid)
    return qs, gts


def peak_rss_mb():
    try:
        import psutil

        return psutil.Process().memory_info().rss / (1024 * 1024)
    except Exception:
        return 0.0


def run(name, pool, n_refs, n_queries, n_subs, threads, rng):
    refs = make_refs(pool, n_refs, rng)
    queries, gts = make_queries(refs, n_queries, n_subs, rng)

    t0 = time.perf_counter()
    idx = seqtree.Index.build(refs, alphabet="aa")
    build_s = time.perf_counter() - t0

    p = seqtree.SearchParams(max_subs=n_subs, engine="seqtm")
    t0 = time.perf_counter()
    results = idx.search_batch(queries, p, threads=threads)
    search_s = time.perf_counter() - t0

    in_set = top1 = 0
    for hits, gt in zip(results, gts):
        ids = {h.ref_id for h in hits}
        if gt in ids:
            in_set += 1
        if hits:
            best = min(hits, key=lambda h: (h.score, h.ref_id))
            if best.ref_id == gt:
                top1 += 1

    qps = n_queries / search_s if search_s else float("inf")
    print(
        f"{name}\t{n_refs}\t{n_queries}\t{n_subs}\t{build_s:.2f}\t"
        f"{qps:,.0f}\t{in_set / n_queries:.3f}\t{top1 / n_queries:.3f}\t{peak_rss_mb():.0f}"
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="*", help="reference-set sizes")
    ap.add_argument("--queries", type=int, default=None)
    ap.add_argument("--subs", type=int, default=1, help="substitutions per query")
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    cdr3, epi = load_pools()
    print(f"# loaded {len(cdr3)} CDR3 and {len(epi)} epitope sequences", flush=True)

    slow = os.environ.get("RUN_BENCHMARK")
    if args.sizes:
        sizes = args.sizes
    elif slow:
        sizes = [1_000_000]
    else:
        sizes = [min(len(cdr3), 20_000)]
    nq = args.queries if args.queries is not None else (1_000_000 if slow else 5_000)

    print("pool\tn_refs\tn_queries\tsubs\tbuild_s\tqps\trecall_set\trecall_top1\trss_mb")
    for sz in sizes:
        run("cdr3", cdr3, sz, min(nq, sz * 10 or nq), args.subs, args.threads, rng)
    if epi:
        run("epitope", epi, min(len(epi), sizes[0]), min(nq, 5_000), args.subs, args.threads, rng)


if __name__ == "__main__":
    main()
