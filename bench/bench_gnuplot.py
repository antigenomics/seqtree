#!/usr/bin/env python3
"""Multi-axis seqtree benchmark on two TCR-beta reference families, rendered to SVG.

Two reference families are measured **separately**:

  * ``olga``  -- OLGA-generated human TRB CDR3 (generative, *no antigen motif*)
  * ``vdjdb`` -- VDJdb CDR3 mutated (antigen-specific, *with motif* / shared structure)

Queries are 1000 fresh OLGA TRB (olga family) / 1000 held-out VDJdb (vdjdb family);
timings are over those 1000 queries. Figures are vertically stacked 2-panel SVGs
written to ``bench/figures/<key>.svg`` (+ ``<key>_<panel>.tsv``):

  scaling                 throughput vs reference size, per engine x threads (top olga / bottom vdjdb)
  scope                   throughput + matches per query vs edit budget 1..5
  matrix                  seqtm throughput: unit vs BLOSUM62 vs PAM50 (top olga / bottom vdjdb)
  selectivity_collisions  matches vs penalty budget (top); seqtm collisions vs edit budget (bottom)
  perop                   align() CIGAR fetch us/call (top); peak RSS (bottom)

Engine line styles: seqtm = long dash, seqtrie = dash-dot. Dataset = colour.
The "all reference sizes" figures (scaling/matrix/perop) span 10k..10M; the
edit-budget sweeps (scope/selectivity/collisions) run at one representative size.

  python bench/bench_gnuplot.py                       # fast tier (10k/100k)
  env RUN_BENCHMARK=1 python bench/bench_gnuplot.py   # full tier: 10k/100k/1M/10M

Needs gnuplot + olga-generate_sequences on PATH and ``pip install -e ".[bench]"``.
"""
import argparse
import csv
import gzip
import os
import random
import shutil
import subprocess
import tempfile
import time
from pathlib import Path

import seqtree as st

AA = "ACDEFGHIKLMNPQRSTVWY"
DATASETS = ("olga", "vdjdb")

COLOR_DS = {"olga": "#1f77b4", "vdjdb": "#d62728"}      # dataset -> colour
COLOR_ENG = {"seqtm": "#d62728", "seqtrie": "#1f77b4"}  # engine -> colour (scaling panels)
COLOR_M = {"unit": "#2ca02c", "BLOSUM62": "#1f77b4", "PAM50": "#d62728"}
DT = {"seqtm": "(25,12)", "seqtrie": "(16,8,3,8)"}      # engine -> dashtype (long-dash / dash-dot)
PT = {1: 6, 4: 4, 8: 8}                                 # threads -> point type


def style(color, engine=None, pt=7):
    dt = f" dt {DT[engine]}" if engine else ""
    return f"with linespoints lw 2 lc rgb '{color}'{dt} pt {pt} ps 0.8"


# --- data --------------------------------------------------------------------
def olga_trb(n, seed=None):
    if not shutil.which("olga-generate_sequences"):
        raise SystemExit("olga-generate_sequences not on PATH (pip install olga)")
    d = tempfile.mkdtemp()
    path = os.path.join(d, "olga.tsv")  # OLGA refuses to overwrite an existing file
    cmd = ["olga-generate_sequences", "--humanTRB", "-n", str(n), "-o", path]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    subprocess.run(cmd, check=True, capture_output=True, timeout=2400)
    pool = [ln.split("\t")[1].strip() for ln in open(path) if "\t" in ln]
    shutil.rmtree(d, ignore_errors=True)
    return [s for s in pool if s and all(c in AA for c in s)]


def vdjdb_cdr3():
    try:
        from huggingface_hub import hf_hub_download

        csv.field_size_limit(10**7)
        path = hf_hub_download(repo_id="isalgo/airr_benchmark",
                               filename="vdjdb/vdjdb-2025-12-29/vdjdb.slim.txt.gz",
                               repo_type="dataset")
        out = set()
        with gzip.open(path, "rt") as fh:
            r = csv.DictReader(fh, delimiter="\t")
            col = next((c for c in ("cdr3", "junction_aa", "cdr3_aa") if c in (r.fieldnames or [])), None)
            for row in r:
                v = (row.get(col) or "").strip().upper() if col else ""
                if v and all(c in AA for c in v):
                    out.add(v)
        return sorted(out)
    except Exception:
        return []


def mutate(s, n_subs, rng):
    s = list(s)
    for _ in range(n_subs):
        j = rng.randrange(len(s))
        s[j] = rng.choice(AA)
    return "".join(s)


def make_refs(kind, size, olga_pool, vdjdb_pool, rng):
    """OLGA refs are the generative pool; VDJdb refs are mutated VDJdb CDR3. Both are
    expanded to ``size`` by substitution-mutating random picks from the source pool."""
    src = olga_pool if kind == "olga" else vdjdb_pool
    refs = list(src) if kind == "olga" else [mutate(s, rng.randint(1, 2), rng) for s in src]
    while len(refs) < size:
        refs.append(mutate(rng.choice(src), rng.randint(0, 3), rng))
    rng.shuffle(refs)
    return refs[:size]


def peak_rss_mb():
    try:
        import psutil

        return psutil.Process().memory_info().rss / (1024 * 1024)
    except Exception:
        return 0.0


def qpms(idx, queries, params, threads):
    t0 = time.perf_counter()
    res = idx.search_batch(queries, params, threads=threads)
    dt_ms = (time.perf_counter() - t0) * 1000.0
    return (len(queries) / dt_ms if dt_ms else 0.0, sum(len(r) for r in res) / len(queries))


# --- gnuplot (vertically stacked panels) -------------------------------------
def write_tsv(path, xs, series):
    head = "x\t" + "\t".join(lbl for lbl, _, _ in series)
    rows = []
    for i, x in enumerate(xs):
        cells = [str(x)] + [("NaN" if s[1][i] == "" else f"{s[1][i]:g}") for s in series]
        rows.append("\t".join(cells))
    path.write_text(head + "\n" + "\n".join(rows) + "\n")


def render(out: Path, key, panels, width=720, panel_h=360):
    """panels: list of dicts(title,xlabel,ylabel,xs,series,logx?,logy?,xtics?)."""
    n = len(panels)
    lines = [f"set terminal svg size {width},{panel_h * n} font 'Helvetica,12' background rgb 'white'",
             f"set output '{key}.svg'", 'set datafile separator "\\t"',
             'set datafile missing "NaN"', "set grid", "set key outside right top"]
    if n > 1:
        lines.append(f"set multiplot layout {n},1")
    for pi, p in enumerate(panels):
        tsv = f"{key}_{pi}.tsv"
        write_tsv(out / tsv, p["xs"], p["series"])
        lines += [f'set title "{p["title"]}"', f"set xlabel '{p['xlabel']}'",
                  f"set ylabel '{p['ylabel']}'",
                  "set logscale x" if p.get("logx") else "unset logscale x",
                  "set logscale y" if p.get("logy") else "unset logscale y"]
        if p.get("xtics"):
            lines.append("set xtics (" + ", ".join(f"'{l}' {v}" for l, v in p["xtics"]) + ")")
        else:
            lines.append("set xtics auto")
        # An empty series label -> notitle (drops it from the legend, e.g. for the
        # dotted "predicted" companions that share a colour with their dashed "observed").
        lines.append("plot " + ", ".join(
            f"'{tsv}' using 1:{i + 2} {stl} "
            + (f"title columnheader({i + 2})" if lbl else "notitle")
            for i, (lbl, _, stl) in enumerate(p["series"])))
    if n > 1:
        lines.append("unset multiplot")
    (out / f"{key}.gp").write_text("\n".join(lines) + "\n")
    subprocess.run(["gnuplot", f"{key}.gp"], cwd=out, check=True)


# --- main --------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="*")
    ap.add_argument("--queries", type=int, default=1000)
    ap.add_argument("--threads", type=int, nargs="*", default=[1, 4, 8])
    ap.add_argument("--max-edits", type=int, default=5)
    ap.add_argument("--sweep-size", type=int, default=None, help="size for edit-budget sweeps")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="bench/figures")
    args = ap.parse_args()

    if not shutil.which("gnuplot"):
        raise SystemExit("gnuplot not found on PATH (brew install gnuplot / conda install gnuplot)")
    rng = random.Random(args.seed)
    slow = os.environ.get("RUN_BENCHMARK")
    sizes = args.sizes or ([10_000, 100_000, 1_000_000, 10_000_000] if slow else [10_000, 100_000])
    edits = list(range(1, args.max_edits + 1))
    penalties = [4, 8, 12, 16, 20, 24]
    sweep_size = args.sweep_size or min(max(sizes), 100_000)  # keep edit sweeps tractable
    tmax = max(args.threads)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    olga_n = min(max(sizes[0], 50_000), 200_000)
    print(f"# OLGA pool ({olga_n}) + query sets; sizes={sizes}; sweep@{sweep_size:,}; threads={args.threads}",
          flush=True)
    olga_pool = olga_trb(olga_n, seed=args.seed)
    vdjdb_pool = vdjdb_cdr3()
    if not vdjdb_pool:
        raise SystemExit("VDJdb pool empty (need cached HuggingFace dataset for the vdjdb family)")
    qset = {"olga": olga_trb(args.queries, seed=args.seed + 1)[: args.queries],
            "vdjdb": random.Random(args.seed + 2).sample(vdjdb_pool, min(args.queries, len(vdjdb_pool)))}
    print(f"# olga pool {len(olga_pool)}, vdjdb pool {len(vdjdb_pool)}", flush=True)

    # accumulators
    scaling = {ds: {(e, t): [] for e in ("seqtm", "seqtrie") for t in args.threads} for ds in DATASETS}
    mat = {ds: {m: [] for m in ("unit", "BLOSUM62", "PAM50")} for ds in DATASETS}
    rss = {ds: [] for ds in DATASETS}
    align = {ds: [] for ds in DATASETS}
    scope_q = {ds: {"seqtm": [], "seqtrie": []} for ds in DATASETS}    # throughput vs edits
    scope_h = {ds: {"seqtm": [], "seqtrie": []} for ds in DATASETS}    # matches/q vs edits
    sel_sc = {ds: {"BLOSUM62": [], "PAM50": []} for ds in DATASETS}    # matches/q vs penalty
    coll = {ds: [] for ds in DATASETS}                                 # collisions/q vs edits

    for ds in DATASETS:
        queries = qset[ds]
        for n in sizes:
            refs = make_refs(ds, n, olga_pool, vdjdb_pool, rng)
            t0 = time.perf_counter()
            idx = st.Index.build(refs, alphabet="aa")
            rss[ds].append(round(peak_rss_mb(), 1))
            print(f"# [{ds} {n:,}] built in {time.perf_counter()-t0:.1f}s, RSS {rss[ds][-1]:.0f} MB", flush=True)

            for eng in ("seqtm", "seqtrie"):
                for t in args.threads:
                    q, _ = qpms(idx, queries, st.SearchParams(max_subs=2, max_total_edits=2, engine=eng), t)
                    scaling[ds][(eng, t)].append(round(q, 3))
            for m in ("unit", "BLOSUM62", "PAM50"):
                kw = dict(max_subs=3, engine="seqtm")
                if m != "unit":
                    kw["matrix"] = m
                q, _ = qpms(idx, queries, st.SearchParams(**kw), tmax)
                mat[ds][m].append(round(q, 3))
            ap_ = st.SearchParams(matrix="BLOSUM62", max_penalty=20, gap_open=8, engine="seqtrie")
            res = idx.search_batch(queries, ap_, threads=tmax)
            pairs = [(h.ref_id, qy) for qy, hits in zip(queries, res) for h in hits][:2000]
            t0 = time.perf_counter()
            for rid, qy in pairs:
                idx.align(rid, qy, ap_).ops
            align[ds].append(round((time.perf_counter() - t0) * 1e6 / max(1, len(pairs)), 3))

            if n == sweep_size:
                for e in edits:
                    q1, h1 = qpms(idx, queries, st.SearchParams(max_subs=e, max_total_edits=e, engine="seqtm"), tmax)
                    q2, h2 = qpms(idx, queries, st.SearchParams(max_total_edits=e, engine="seqtrie"), tmax)
                    scope_q[ds]["seqtm"].append(round(q1, 3)); scope_q[ds]["seqtrie"].append(round(q2, 3))
                    scope_h[ds]["seqtm"].append(round(h1, 2)); scope_h[ds]["seqtrie"].append(round(h2, 2))
                    c = idx.collisions_batch(
                        queries, st.SearchParams(max_subs=e, max_ins=1, max_dels=1, max_total_edits=e, engine="seqtm"), tmax)
                    coll[ds].append(round(sum(c) / len(queries), 2))
                for p in penalties:
                    _, hb = qpms(idx, queries, st.SearchParams(matrix="BLOSUM62", max_penalty=p, gap_open=8, engine="seqtrie"), tmax)
                    _, hp = qpms(idx, queries, st.SearchParams(matrix="PAM50", max_penalty=p, gap_open=8, engine="seqtrie"), tmax)
                    sel_sc[ds]["BLOSUM62"].append(round(hb, 2)); sel_sc[ds]["PAM50"].append(round(hp, 2))
                print(f"# [{ds} {n:,}] sweep done: e5 seqtm {scope_q[ds]['seqtm'][-1]} q/ms, "
                      f"coll/q {coll[ds]}", flush=True)

    # --- render figures ---
    et = [(str(e), e) for e in edits]
    render(out, "scaling", [
        {"title": f"{ds}: throughput vs reference size (2 substitutions)", "xlabel": "reference set size",
         "ylabel": "queries / ms", "xs": sizes, "logx": True,
         "series": [(f"{eng} t={t}", scaling[ds][(eng, t)], style(COLOR_ENG[eng], eng, PT[t]))
                    for eng in ("seqtm", "seqtrie") for t in args.threads]}
        for ds in DATASETS])
    render(out, "matrix", [
        {"title": f"{ds}: seqtm matrix-scoring cost (3 substitutions)", "xlabel": "reference set size",
         "ylabel": "queries / ms", "xs": sizes, "logx": True,
         "series": [(f"seqtm {m}", mat[ds][m], style(COLOR_M[m], "seqtm", 6 + 2 * i))
                    for i, m in enumerate(("unit", "BLOSUM62", "PAM50"))]}
        for ds in DATASETS])
    render(out, "perop", [
        {"title": "seqtrie global-alignment CIGAR fetch (C++)", "xlabel": "reference set size",
         "ylabel": "us / align() call", "xs": sizes, "logx": True,
         "series": [(ds, align[ds], style(COLOR_DS[ds])) for ds in DATASETS]},
        {"title": "peak RSS after index build", "xlabel": "reference set size",
         "ylabel": "peak RSS (MB)", "xs": sizes, "logx": True, "logy": True,
         "series": [(ds, rss[ds], style(COLOR_DS[ds])) for ds in DATASETS]}])
    render(out, "scope", [
        {"title": f"throughput vs edit budget ({sweep_size:,} refs, t={tmax})", "xlabel": "edit budget",
         "ylabel": "queries / ms", "xs": edits, "logy": True, "xtics": et,
         "series": [(f"{ds} {eng}", scope_q[ds][eng], style(COLOR_DS[ds], eng))
                    for ds in DATASETS for eng in ("seqtm", "seqtrie")]},
        {"title": f"matches per query vs edit budget ({sweep_size:,} refs)", "xlabel": "edit budget",
         "ylabel": "matches / query", "xs": edits, "logy": True, "xtics": et,
         "series": [(f"{ds} {eng}", scope_h[ds][eng], style(COLOR_DS[ds], eng))
                    for ds in DATASETS for eng in ("seqtm", "seqtrie")]}])
    render(out, "selectivity_collisions", [
        {"title": f"matches per query vs score budget ({sweep_size:,} refs, gap_open=8)",
         "xlabel": "max penalty budget", "ylabel": "matches / query", "xs": penalties, "logy": True,
         "series": [(f"{ds} {m}", sel_sc[ds][m], style(COLOR_DS[ds], "seqtrie", 6 if m == "BLOSUM62" else 8))
                    for ds in DATASETS for m in ("BLOSUM62", "PAM50")]},
        {"title": f"seqtm collisions vs edit budget ({sweep_size:,} refs, +/-1 indel)",
         "xlabel": "edit budget", "ylabel": "collisions / query", "xs": edits, "logy": True, "xtics": et,
         "series": [(ds, coll[ds], style(COLOR_DS[ds], "seqtm")) for ds in DATASETS]}])
    print(f"\nWrote 5 stacked SVG figures + TSVs to {out}/", flush=True)


if __name__ == "__main__":
    main()
