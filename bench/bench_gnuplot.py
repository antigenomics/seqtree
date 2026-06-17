#!/usr/bin/env python3
"""Simple max-edit-3 throughput benchmark for seqtree, rendered to SVG with gnuplot.

Queries are references mutated by up to 3 edits. Throughput is reported as
**queries per millisecond** against reference sets of increasing size, for both
engines under matched scope:

  * seqtm 3 subs            -- max_subs=3 (Hamming-style)
  * seqtm 1 sub + 1 indel   -- max_subs=1, max_ins=1, max_dels=1
  * seqtrie edits<=3        -- max_total_edits=3 (banded DP)
  * seqtrie BLOSUM62        -- max_penalty=15 (~3 conservative subs)

Plus two figures: peak RSS after index build, and the per-call cost of fetching
an alignment CIGAR. Each figure -> ``<out>/<key>.svg`` + the ``<key>.tsv`` behind it.

  python bench/bench_gnuplot.py                          # fast tier (small, seconds)
  env RUN_BENCHMARK=1 python bench/bench_gnuplot.py      # full tier: 100k/1M/10M, 1M queries

Needs gnuplot on PATH and ``pip install -e ".[bench]"`` (psutil; huggingface_hub
optional). ``--random`` skips the HuggingFace CDR3 fetch and uses seeded random seqs.
"""
import argparse
import os
import random
import shutil
import subprocess
import time
from pathlib import Path

import seqtree

AA = "ACDEFGHIKLMNPQRSTVWY"

# Fixed max-edit-3 configs: (label, engine, params-kwargs, (n_subs, n_indels) for queries).
CONFIGS = [
    ("seqtm 3 subs", dict(engine="seqtm", max_subs=3), (3, 0)),
    ("seqtm 1 sub+1 indel", dict(engine="seqtm", max_subs=1, max_ins=1, max_dels=1), (1, 1)),
    ("seqtrie edits<=3", dict(engine="seqtrie", max_total_edits=3), (2, 1)),
    ("seqtrie BLOSUM62 p<=15", dict(engine="seqtrie", matrix="BLOSUM62", max_penalty=15, gap_open=8), (3, 0)),
]


def gen_pool(n, rng, use_hf=True):
    """Real CDR3 (aa) from the cached VDJdb slim table, else seeded random seqs."""
    if not use_hf:
        return ["".join(rng.choice(AA) for _ in range(rng.randint(12, 18))) for _ in range(min(n, 50_000))]
    try:
        import csv
        import gzip

        from huggingface_hub import hf_hub_download

        csv.field_size_limit(10**7)
        path = hf_hub_download(
            repo_id="isalgo/airr_benchmark",
            filename="vdjdb/vdjdb-2025-12-29/vdjdb.slim.txt.gz",
            repo_type="dataset",
        )
        pool = set()
        with gzip.open(path, "rt") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            col = next((c for c in ("cdr3", "junction_aa", "cdr3_aa") if c in (reader.fieldnames or [])), None)
            for row in reader:
                v = (row.get(col) or "").strip().upper() if col else ""
                if v and all(c in AA for c in v):
                    pool.add(v)
        if pool:
            return sorted(pool)
    except Exception:
        pass
    return ["".join(rng.choice(AA) for _ in range(rng.randint(12, 18))) for _ in range(min(n, 50_000))]


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


def peak_rss_mb():
    try:
        import psutil

        return psutil.Process().memory_info().rss / (1024 * 1024)
    except Exception:
        return 0.0


def write_svg(out: Path, key: str, title: str, ylabel: str, sizes, labels, rows, logy=False):
    """rows[i] = [size_i, val_label_0, ...]; '' -> NaN (skipped by gnuplot)."""
    tsv = out / f"{key}.tsv"
    body = "\n".join("\t".join(("NaN" if v == "" else f"{v:g}") if isinstance(v, float) else str(v)
                                for v in r) for r in rows)
    tsv.write_text("size\t" + "\t".join(labels) + "\n" + body + "\n")

    plot = ", ".join(f"'{tsv.name}' using 1:{i + 2} with linespoints lw 2 pt 7 ps 0.8 "
                     f"title columnheader({i + 2})" for i in range(len(labels)))
    (out / f"{key}.gp").write_text(
        f"set terminal svg size 760,480 font 'Helvetica,13' background rgb 'white'\n"
        f"set output '{key}.svg'\n"
        f'set title "{title}"\n'
        f"set xlabel 'reference set size'\nset ylabel '{ylabel}'\n"
        f"set logscale x\n{'set logscale y' if logy else ''}\n"
        f"set grid\nset key outside right top\n"
        f'set datafile separator "\\t"\nset datafile missing "NaN"\n'
        f"set xtics ({', '.join(f'{chr(39)}{s:,}{chr(39)} {s}' for s in sizes)})\n"
        f"plot {plot}\n"
    )
    subprocess.run(["gnuplot", f"{key}.gp"], cwd=out, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="*")
    ap.add_argument("--queries", type=int, default=None)
    ap.add_argument("--threads", type=int, default=0)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="bench/figures")
    ap.add_argument("--random", action="store_true", help="seeded random pool (skip HuggingFace fetch)")
    args = ap.parse_args()

    if not shutil.which("gnuplot"):
        raise SystemExit("gnuplot not found on PATH (brew install gnuplot / conda install gnuplot)")

    rng = random.Random(args.seed)
    slow = os.environ.get("RUN_BENCHMARK")
    sizes = args.sizes or ([100_000, 1_000_000, 10_000_000] if slow else [10_000, 100_000])
    nq = args.queries if args.queries is not None else (1_000_000 if slow else 5_000)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    pool = gen_pool(max(sizes), rng, use_hf=not args.random)
    print(f"# pool {len(pool)} seqs; sizes={sizes}; queries={nq}; threads={args.threads or 'all'}", flush=True)

    tput = {lbl: [] for lbl, _, _ in CONFIGS}
    rss_vals, align_vals = [], []

    for n in sizes:
        refs = expand(pool, n, rng)
        t0 = time.perf_counter()
        idx = seqtree.Index.build(refs, alphabet="aa")
        rss_vals.append(round(peak_rss_mb(), 1))
        print(f"# built {n:,} refs in {time.perf_counter() - t0:.1f}s, RSS {rss_vals[-1]:.0f} MB", flush=True)
        base = [refs[rng.randrange(n)] for _ in range(min(nq, n))]

        qcache = {}  # one query list per distinct (subs, indels); generation dwarfs the C++ search
        for lbl, kw, (ns, ni) in CONFIGS:
            if (ns, ni) not in qcache:
                qcache[(ns, ni)] = [mutate(s, ns, ni, rng) for s in base]
            queries = qcache[(ns, ni)]
            t0 = time.perf_counter()
            idx.search_batch(queries, seqtree.SearchParams(**kw), threads=args.threads)
            dt_ms = (time.perf_counter() - t0) * 1000.0
            qpms = len(queries) / dt_ms if dt_ms else 0.0
            tput[lbl].append(round(qpms, 3))
            print(f"  {n:>10,}  {lbl:<22} {qpms:9.2f} q/ms", flush=True)

        # alignment CIGAR fetch cost (us per align() call returning the ops string)
        aq = [mutate(s, 1, 0, rng) for s in base[: min(5000, len(base))]]
        p = seqtree.SearchParams(max_subs=1, engine="seqtm")
        res = idx.search_batch(aq, p, threads=args.threads)
        pairs = [(h.ref_id, q) for q, hits in zip(aq, res) for h in hits]
        t0 = time.perf_counter()
        for ref_id, q in pairs:
            idx.align(ref_id, q, p).ops
        us = (time.perf_counter() - t0) * 1e6 / len(pairs) if pairs else 0.0
        align_vals.append(round(us, 3))
        print(f"  {n:>10,}  align CIGAR fetch       {us:9.3f} us/call ({len(pairs)} calls)", flush=True)

    labels = [lbl for lbl, _, _ in CONFIGS]
    rows = [[sizes[i]] + [tput[lbl][i] for lbl in labels] for i in range(len(sizes))]
    write_svg(out, "throughput", "Throughput at max edit distance 3", "queries / ms", sizes, labels, rows)
    write_svg(out, "ram", "Peak RSS after index build", "peak RSS (MB)", sizes, ["peak RSS"],
              [[sizes[i], rss_vals[i]] for i in range(len(sizes))], logy=True)
    write_svg(out, "align_fetch", "Alignment CIGAR fetch cost", "us / align() call", sizes, ["align()"],
              [[sizes[i], align_vals[i]] for i in range(len(sizes))])
    print(f"\nWrote 3 SVG figures + TSVs to {out}/", flush=True)


if __name__ == "__main__":
    main()
