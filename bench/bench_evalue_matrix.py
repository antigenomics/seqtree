#!/usr/bin/env python3
"""Comprehensive E-value benchmark: reference structure x control size x query x scope.

Reference sets (target D, all built to the same size N so E=(N/M)*n_control is comparable):
  vdjdb        -- VDJdb CDR3 (antigen-selected, clustered)
  olga         -- OLGA-generated TRB (generation null, weak structure)
  vdjdb+noise  -- half VDJdb, half uniform-random noise (clusters diluted 2x)
  olga+noise   -- half OLGA, half uniform-random noise
Background controls: OLGA-generated, sizes 1M / 2M / 10M (from bench/gen_olga.py cache).
Query sets: vdjdb or olga. Scope: 1..3 substitutions.

For each (ref, query, control, scope) it reports mean neighbours (exact excluded), median
E-value, and the fraction of neighbours called significant at a fixed E<1 cutoff and after
Benjamini-Hochberg FDR<0.05. Emits a full TSV table and a stacked SVG (one panel per query
set: fraction significant vs scope, one line per reference set, at the largest control).

    python bench/bench_evalue_matrix.py                       # control = 1M (fast)
    env RUN_BENCHMARK=1 python bench/bench_evalue_matrix.py   # controls 1M/2M/10M

Needs gnuplot, cached OLGA controls, and the VDJdb dataset.
"""
import argparse
import gzip
import os
import random
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import seqtree as st
from seqtree.evalue import _poisson_sf

from bench_evalue import bh_reject
from bench_gnuplot import mutate, render, style, vdjdb_cdr3

AA = "ACDEFGHIKLMNPQRSTVWY"
REF_SETS = ("vdjdb", "olga", "vdjdb+noise", "olga+noise")
COLOR_REF = {"vdjdb": "#d62728", "olga": "#1f77b4", "vdjdb+noise": "#ff9896", "olga+noise": "#aec7e8"}
CONTROLS = {"1M": "bench/cache/olga_1M.txt.gz", "2M": "bench/cache/olga_2M.txt.gz",
            "10M": "bench/cache/olga_10M.txt.gz"}


def load_cache(path):
    with gzip.open(path, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def rand_seqs(n, rng):
    return ["".join(rng.choice(AA) for _ in range(rng.randint(12, 17))) for _ in range(n)]


def make_ref(kind, size, vdjdb_pool, olga_pool, rng):
    half = size // 2
    if kind == "vdjdb":
        base = rng.sample(vdjdb_pool, min(size, len(vdjdb_pool)))
    elif kind == "olga":
        base = rng.sample(olga_pool, min(size, len(olga_pool)))
    elif kind == "vdjdb+noise":
        base = rng.sample(vdjdb_pool, half) + rand_seqs(size - half, rng)
    else:  # olga+noise
        base = rng.sample(olga_pool, half) + rand_seqs(size - half, rng)
    return list(dict.fromkeys(base))  # unique clonotypes


def neighbour_counts(idx, queries, scope, threads):
    """Per-query count of distinct hits with positive distance (exact/self excluded)."""
    res = idx.search_batch(queries, st.SearchParams(max_subs=scope, engine="seqtm"), threads)
    return [sum(1 for h in r if h.score > 0) for r in res]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", type=int, default=20000, help="reference set size N")
    ap.add_argument("--queries", type=int, default=1000)
    ap.add_argument("--scopes", type=int, nargs="*", default=[1, 2, 3])
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="bench/figures")
    args = ap.parse_args()

    slow = os.environ.get("RUN_BENCHMARK")
    ctrl_names = ["1M", "2M", "10M"] if slow else ["1M"]
    ctrl_names = [c for c in ctrl_names if os.path.exists(CONTROLS[c])]
    rng = random.Random(args.seed)

    vdjdb_pool = vdjdb_cdr3()
    olga_pool = load_cache(CONTROLS["1M"])
    if not vdjdb_pool:
        raise SystemExit("VDJdb pool empty")

    # query sets
    qsets = {"vdjdb": random.Random(1).sample(vdjdb_pool, min(args.queries, len(vdjdb_pool))),
             "olga": random.Random(2).sample(olga_pool, args.queries)}

    # build controls and precompute n_control[qset][ctrl][scope]
    controls = {c: st.Index.build(load_cache(CONTROLS[c]), alphabet="aa") for c in ctrl_names}
    print(f"# controls: " + ", ".join(f"{c}={len(controls[c]):,}" for c in ctrl_names), flush=True)
    n_ctrl = {qs: {c: {} for c in ctrl_names} for qs in qsets}
    for qs, qq in qsets.items():
        for c in ctrl_names:
            for sc in args.scopes:
                n_ctrl[qs][c][sc] = neighbour_counts(controls[c], qq, sc, args.threads)
        print(f"# control counts done for query={qs}", flush=True)

    # build references and precompute n_target[ref][qset][scope]
    refs = {r: st.Index.build(make_ref(r, args.size, vdjdb_pool, olga_pool, rng), alphabet="aa")
            for r in REF_SETS}
    n_tgt = {r: {qs: {} for qs in qsets} for r in REF_SETS}
    for r in REF_SETS:
        for qs, qq in qsets.items():
            for sc in args.scopes:
                n_tgt[r][qs][sc] = neighbour_counts(refs[r], qq, sc, args.threads)
    print(f"# reference sizes: " + ", ".join(f"{r}={len(refs[r])}" for r in REF_SETS), flush=True)

    # assemble table
    print("\nref\tquery\tcontrol\tscope\tN\tM\tmean_nbr\tmedian_E\tfracSig_E<1\tfracSig_BH<0.05")
    table = {}  # (ref,qset,ctrl,scope) -> fracBH
    for r in REF_SETS:
        N = len(refs[r])
        for qs in qsets:
            for c in ctrl_names:
                M = len(controls[c])
                for sc in args.scopes:
                    nt = n_tgt[r][qs][sc]
                    nc = n_ctrl[qs][c][sc]
                    Es, ps = [], []
                    for a, b in zip(nt, nc):
                        E = (3.0 if b == 0 else float(b)) * N / M
                        Es.append(E)
                        ps.append(_poisson_sf(a, E))
                    total = sum(nt) or 1
                    frac_e = sum(a for a, E in zip(nt, Es) if E < 1.0) / total
                    rej = bh_reject(ps, 0.05)
                    frac_bh = sum(a for a, rj in zip(nt, rej) if rj) / total
                    med_e = sorted(Es)[len(Es) // 2]
                    table[(r, qs, c, sc)] = frac_bh
                    print(f"{r}\t{qs}\t{c}\t{sc}\t{N}\t{M}\t{sum(nt)/len(nt):.2f}\t"
                          f"{med_e:.3g}\t{frac_e:.3f}\t{frac_bh:.3f}", flush=True)

    # figure: one panel per query set; x=scope, y=fracSig(BH); line per ref set; largest control
    cbig = ctrl_names[-1]
    panels = []
    for qs in qsets:
        panels.append({
            "title": f"query={qs}: significant fraction vs scope (control OLGA {cbig}, BH FDR<0.05)",
            "xlabel": "scope (substitutions)", "ylabel": "fraction significant",
            "xs": args.scopes, "xtics": [(str(s), s) for s in args.scopes],
            "series": [(r, [table[(r, qs, cbig, sc)] for sc in args.scopes], style(COLOR_REF[r], "seqtm"))
                       for r in REF_SETS]})
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    render(out, "evalue_matrix", panels)
    print(f"\nWrote {out}/evalue_matrix.svg")


if __name__ == "__main__":
    main()
