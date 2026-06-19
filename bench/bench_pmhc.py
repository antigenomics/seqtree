#!/usr/bin/env python3
"""pMHC epitope homology benchmark.

Validates the Dolton et al. (Cell 2023) HLA-A*02:01 cross-reactive trio as mutual
TCR-facing homologs, then (if the isalgo/pmhc_data table is available) builds a real
A*02:01 MHC-I store and reports homology-search throughput, plus a neoantigen
mimic-discovery example with presentation-aware E-values.

    python bench/bench_pmhc.py
    python bench/bench_pmhc.py --pmhc /Users/.../pmhc_full.tsv.gz --allele HLA-A*02:01

Needs `pip install -e ".[pmhc]"` (huggingface_hub) only for the auto-download path.
"""
import argparse
import gzip
import csv
import os
import time

from seqtree import pmhc

TRIPLE = ["EAAGIGILTV", "LLLGIGILVL", "NLSALGIFST"]  # Melan-A / BST2 / IMP2, A*02:01


def cell_triple():
    recs = [{"epitope": e, "mhc": "HLA-A*02:01", "mhc_class": "MHCI", "gene": g}
            for e, g in zip(TRIPLE, ["MLANA", "BST2", "IGF2BP2"])]
    store = pmhc.PMHCStore.from_records(recs, k=4)
    print("# Dolton A*02:01 trio — mutual TCR-facing homologs (k=4, max_subs=2):")
    for q in TRIPLE:
        hits = store.search_homologs(q, "mhc1", mhc="HLA-A*02:01", max_subs=2, min_shared=1)
        print(f"  {q}: " + ", ".join(f"{h.epitope}(n={h.shared_kmers},s={h.score})" for h in hits))


def load_pmhc_table(path):
    csv.field_size_limit(10**7)
    op = gzip.open if str(path).endswith(".gz") else open
    rows = []
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default="/Users/mikesh/hf/pmhc_data/pmhc_full.tsv.gz")
    ap.add_argument("--allele", default="HLA-A*02:01")
    ap.add_argument("--queries", type=int, default=2000)
    args = ap.parse_args()

    cell_triple()

    if not os.path.exists(args.pmhc):
        print(f"\n# {args.pmhc} not found — skipping the real-data benchmark.")
        return

    rows = [r for r in load_pmhc_table(args.pmhc) if str(r.get("mhc_class")) == "MHCI"]
    t0 = time.perf_counter()
    store = pmhc.PMHCStore.from_records(rows, k=4)
    build_s = time.perf_counter() - t0
    n1 = store.size("mhc1")
    print(f"\n# MHC-I store: {n1} unique epitopes, built in {build_s:.1f}s")

    import random
    rng = random.Random(0)
    allele_eps = [r["epitope"].strip().upper() for r in rows
                  if args.allele in str(r.get("mhc_a", ""))]
    allele_eps = list(dict.fromkeys(allele_eps))
    if not allele_eps:
        print(f"# no {args.allele} epitopes found"); return
    q = [rng.choice(allele_eps) for _ in range(min(args.queries, len(allele_eps) * 5))]
    t0 = time.perf_counter()
    total = 0
    for s in q:
        total += len(store.search_homologs(s, "mhc1", mhc=args.allele, max_subs=1, exclude_self=True))
    dt = time.perf_counter() - t0
    print(f"# {args.allele}: {len(allele_eps)} epitopes; {len(q)} queries in {dt:.2f}s "
          f"({len(q)/dt:,.0f} q/s), mean {total/len(q):.1f} homologs/query (max_subs=1)")

    # mimic-discovery example: a normal-complexity 9-mer A*02:01 epitope vs the presented set + a decoy
    neo = next((e for e in allele_eps if len(e) == 9 and len(set(e)) >= 7), allele_eps[0])
    bact = rng.sample(allele_eps, min(2000, len(allele_eps)))
    res = pmhc.find_mimics(neo, self_set=allele_eps[1:5000], bacterial_sets={"decoy": bact},
                           max_subs=1, min_shared=1)
    print(f"# find_mimics({neo}):")
    for src, r in res.items():
        print(f"    {src}: n_hits={len(r['hits'])} E={r['E']:.3g} p_enrich={r['p_enrichment']:.3g}")


if __name__ == "__main__":
    main()
