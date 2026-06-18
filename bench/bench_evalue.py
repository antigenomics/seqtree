#!/usr/bin/env python3
"""True E-value benchmark.

For a target repertoire (VDJdb, antigen-selected) scored against a background
control (airr_control), at each scope/budget we report:

  * raw_paths   -- total edit paths explored (unique hits + collisions)
  * collisions  -- references re-reached via a different edit path (seqtm w/ indels)
  * unique_hits -- distinct references in the ball, summed over queries
  * fracSig     -- fraction of unique hits whose query E-value is below a threshold

Antigen-selected (VDJdb) queries should yield a far higher significant fraction
than background (control) queries -- that contrast is the benchmark.

    python bench/bench_evalue.py
    env RUN_BENCHMARK=1 python bench/bench_evalue.py --target-size 200000 --control-size 2000000
"""
import argparse
import os
import random
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import seqtree
from seqtree.evalue import evalues

from bench import load_pools  # VDJdb CDR3/epitope loader (downloads on first use)

E_THRESHOLDS = (1.0, 0.1, 0.01)


def scopes(slow):
    sp = seqtree.SearchParams
    grid = [
        ("subs=1", sp(max_subs=1, engine="seqtm")),
        ("subs=2", sp(max_subs=2, engine="seqtm")),
        ("subs1+indel1", sp(max_subs=1, max_ins=1, max_dels=1, max_total_edits=1, engine="seqtm")),
        ("blosum<=30", sp(matrix="BLOSUM62", max_penalty=30, gap_open=20, engine="seqtrie")),
    ]
    return grid


def bh_reject(pvals, alpha):
    """Benjamini-Hochberg: boolean mask of hypotheses rejected at FDR <= alpha."""
    n = len(pvals)
    if n == 0:
        return []
    order = sorted(range(n), key=lambda i: pvals[i])
    kmax = 0
    for rank, idx in enumerate(order, start=1):
        if pvals[idx] <= rank / n * alpha:
            kmax = rank
    reject = [False] * n
    for rank, idx in enumerate(order, start=1):
        if rank <= kmax:
            reject[idx] = True
    return reject


def run_block(name, target, control, queries, threads):
    N = len(target)
    print(f"\n## {name}: target N={N}  control M={len(control)}  queries={len(queries)}")
    cols = ["scope", "neighbours", "exact", "collisions", "median_E"]
    cols += [f"fracSig(E<{t})" for t in E_THRESHOLDS] + ["fracSig(BH<0.05)", "qps"]
    print("\t".join(cols))
    for label, params in scopes(False):
        t0 = time.perf_counter()
        res = target.search_batch(queries, params, threads)
        coll = target.collisions_batch(queries, params, threads)
        dt = time.perf_counter() - t0
        # exclude exact (distance-0 / self) hits -- queries may be members of target/control
        ev = evalues(target, control, queries, params, threads, exclude_exact=True)

        nbr = [sum(1 for h in r if h.score > 0) for r in res]   # neighbours (excl. exact)
        n_exact = sum(len(r) for r in res) - sum(nbr)
        total_nbr = sum(nbr) or 1
        Es = sorted(e["E"] for e in ev)
        med_e = Es[len(Es) // 2] if Es else 0.0
        e_fracs = [sum(b for b, e in zip(nbr, ev) if e["E"] < t) / total_nbr for t in E_THRESHOLDS]
        rej = bh_reject([e["p_enrichment"] for e in ev], 0.05)
        bh_frac = sum(b for b, r in zip(nbr, rej) if r) / total_nbr
        qps = len(queries) / dt if dt else float("inf")
        row = [label, sum(nbr), n_exact, sum(coll), f"{med_e:.3g}"]
        row += [f"{f:.3f}" for f in e_fracs] + [f"{bh_frac:.3f}", f"{qps:,.0f}"]
        print("\t".join(str(x) for x in row), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--target-size", type=int, default=None)
    ap.add_argument("--control-size", type=int, default=None)
    ap.add_argument("--n-queries", type=int, default=2000)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    slow = bool(os.environ.get("RUN_BENCHMARK"))
    rng = random.Random(args.seed)

    cdr3, _ = load_pools()
    cdr3 = list(dict.fromkeys(cdr3))  # unique clonotypes
    tsize = args.target_size or (len(cdr3) if slow else min(len(cdr3), 20_000))
    target_seqs = cdr3[:tsize]
    target = seqtree.Index.build(target_seqs, alphabet="aa")

    csize = args.control_size or (2_000_000 if slow else None)  # None = full bundled subset
    control = seqtree.load_control("human_trb_aa", size=csize)

    vdjdb_q = [rng.choice(target_seqs) for _ in range(args.n_queries)]
    bg_q = [control.ref_seq(rng.randrange(len(control))) for _ in range(args.n_queries)]

    run_block("VDJdb queries (antigen-selected)", target, control, vdjdb_q, args.threads)
    run_block("background queries (control)", target, control, bg_q, args.threads)


if __name__ == "__main__":
    main()
