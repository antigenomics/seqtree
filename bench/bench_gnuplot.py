#!/usr/bin/env python3
"""Multi-axis seqtree benchmark on TCR-beta data, rendered to SVG with gnuplot.

References are OLGA-generated human TRB CDR3 (amino acid) plus mutated VDJdb CDR3;
queries are 1000 fresh OLGA TRB sequences. Time is measured over those 1000
queries. Each figure is single-purpose and written to ``bench/figures/<key>.svg``
(+ the ``<key>.tsv`` it was drawn from):

  scaling          throughput (queries/ms) vs reference size, per engine x threads (1/4/8)
  scope            throughput vs edit budget 1..5 (seqtm Hamming ball, seqtrie edit ball)
  matrix_overhead  seqtm throughput: unit vs BLOSUM62 vs PAM50 (cost of matrix scoring)
  align_fetch      seqtrie C++ global-alignment CIGAR fetch cost (us / call)
  ram              peak RSS after index build
  selectivity_scope  matches per query vs edit budget (both engines)
  selectivity_score  matches per query vs BLOSUM62/PAM50 penalty budget (seqtrie)

Engine line styles: seqtm = long dash, seqtrie = dash-dot.

  python bench/bench_gnuplot.py                       # fast tier (10k/100k, seconds-minutes)
  env RUN_BENCHMARK=1 python bench/bench_gnuplot.py   # full tier: 10k/100k/1M/10M

Needs gnuplot on PATH, ``olga-generate_sequences`` on PATH, and
``pip install -e ".[bench]"`` (psutil; huggingface_hub for the VDJdb mix, optional).
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

# Engine line styles (gnuplot dashtype + colour). seqtm long-dash, seqtrie dash-dot.
ENGINE_DT = {"seqtm": "(25,12)", "seqtrie": "(16,8,3,8)"}
ENGINE_LC = {"seqtm": "#d62728", "seqtrie": "#1f77b4"}
THREAD_PT = {1: 6, 4: 4, 8: 8}  # gnuplot point types for the 1/4/8-thread variants


# --- data --------------------------------------------------------------------
def olga_trb(n, seed=None):
    """Generate ``n`` human TRB CDR3 (aa) with OLGA; returns the valid-aa list."""
    if not shutil.which("olga-generate_sequences"):
        raise SystemExit("olga-generate_sequences not on PATH (pip install olga)")
    d = tempfile.mkdtemp()
    path = os.path.join(d, "olga.tsv")  # OLGA refuses to overwrite an existing file
    cmd = ["olga-generate_sequences", "--humanTRB", "-n", str(n), "-o", path]
    if seed is not None:
        cmd += ["--seed", str(seed)]
    subprocess.run(cmd, check=True, capture_output=True, timeout=1200)
    pool = [ln.split("\t")[1].strip() for ln in open(path) if "\t" in ln]
    shutil.rmtree(d, ignore_errors=True)
    return [s for s in pool if s and all(c in AA for c in s)]


def vdjdb_cdr3():
    """VDJdb CDR3 (aa) from the cached HuggingFace slim table; [] if unavailable."""
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


def build_refs(target, olga_pool, vdjdb_pool, rng):
    """Reference set = OLGA TRB + mutated VDJdb CDR3, expanded to ``target`` by
    substitution-mutating random picks (keeps sequences valid amino acid)."""
    refs = list(olga_pool)
    refs += [mutate(s, rng.randint(1, 2), 0, rng) for s in vdjdb_pool[: max(1, target // 4)]]
    while len(refs) < target:
        refs.append(mutate(rng.choice(olga_pool), rng.randint(0, 3), 0, rng))
    rng.shuffle(refs)
    return refs[:target]


def peak_rss_mb():
    try:
        import psutil

        return psutil.Process().memory_info().rss / (1024 * 1024)
    except Exception:
        return 0.0


def qpms(idx, queries, params, threads):
    """Throughput in queries/ms and mean hits/query over ``queries``."""
    t0 = time.perf_counter()
    res = idx.search_batch(queries, params, threads=threads)
    dt_ms = (time.perf_counter() - t0) * 1000.0
    return (len(queries) / dt_ms if dt_ms else 0.0,
            sum(len(r) for r in res) / len(queries))


# --- gnuplot -----------------------------------------------------------------
def plot(out: Path, key, title, xlabel, ylabel, xs, series, logx=True, logy=False, xtics=None):
    """series: list of (label, y_values, style_suffix). '' in y -> NaN (skipped)."""
    tsv = out / f"{key}.tsv"
    header = "x\t" + "\t".join(lbl for lbl, _, _ in series)
    rows = []
    for i, x in enumerate(xs):
        cells = [str(x)] + [("NaN" if s[1][i] == "" else f"{s[1][i]:g}") for s in series]
        rows.append("\t".join(cells))
    tsv.write_text(header + "\n" + "\n".join(rows) + "\n")

    plot_cmd = ", ".join(
        f"'{tsv.name}' using 1:{i + 2} {style} title columnheader({i + 2})"
        for i, (_, _, style) in enumerate(series))
    xtic_line = ""
    if xtics:
        xtic_line = "set xtics (" + ", ".join(f"'{lab}' {v}" for lab, v in xtics) + ")\n"
    (out / f"{key}.gp").write_text(
        "set terminal svg size 820,500 font 'Helvetica,13' background rgb 'white'\n"
        f"set output '{key}.svg'\n"
        f'set title "{title}"\n'
        f"set xlabel '{xlabel}'\nset ylabel '{ylabel}'\n"
        f"{'set logscale x' if logx else ''}\n{'set logscale y' if logy else ''}\n"
        "set grid\nset key outside right top\n"
        'set datafile separator "\\t"\nset datafile missing "NaN"\n'
        f"{xtic_line}"
        f"plot {plot_cmd}\n")
    subprocess.run(["gnuplot", f"{key}.gp"], cwd=out, check=True)


def eng_style(engine, pt=7):
    return (f"with linespoints lw 2 dt {ENGINE_DT[engine]} "
            f"lc rgb '{ENGINE_LC[engine]}' pt {pt} ps 0.9")


# --- main --------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sizes", type=int, nargs="*")
    ap.add_argument("--queries", type=int, default=1000)
    ap.add_argument("--threads", type=int, nargs="*", default=[1, 4, 8])
    ap.add_argument("--max-edits", type=int, default=5)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="bench/figures")
    args = ap.parse_args()

    if not shutil.which("gnuplot"):
        raise SystemExit("gnuplot not found on PATH (brew install gnuplot / conda install gnuplot)")

    rng = random.Random(args.seed)
    slow = os.environ.get("RUN_BENCHMARK")
    sizes = args.sizes or ([10_000, 100_000, 1_000_000, 10_000_000] if slow else [10_000, 100_000])
    edits = list(range(1, args.max_edits + 1))
    penalties = [4, 8, 12, 16, 20, 24]  # BLOSUM/PAM penalty budgets (gap_open=8)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    olga_n = min(max(sizes[0], 50_000), 200_000)
    print(f"# generating OLGA TRB pool ({olga_n}) + 1000 query set ...", flush=True)
    olga_pool = olga_trb(olga_n, seed=args.seed)
    queries = olga_trb(args.queries, seed=args.seed + 1)[: args.queries]
    vdjdb = vdjdb_cdr3()
    print(f"# olga pool {len(olga_pool)}, vdjdb {len(vdjdb)}, queries {len(queries)}; "
          f"sizes={sizes}; threads={args.threads}", flush=True)

    # accumulators keyed by size index
    scaling = {(eng, t): [] for eng in ("seqtm", "seqtrie") for t in args.threads}  # scope=2
    mat_over = {m: [] for m in ("unit", "BLOSUM62", "PAM50")}  # seqtm scope=3, t=max
    rss_vals, align_vals = [], []
    tmax = max(args.threads)
    scope_data = sel_scope = sel_score = None

    for si, n in enumerate(sizes):
        refs = build_refs(n, olga_pool, vdjdb, rng)
        t0 = time.perf_counter()
        idx = st.Index.build(refs, alphabet="aa")
        rss_vals.append(round(peak_rss_mb(), 1))
        print(f"# [{n:,}] built in {time.perf_counter()-t0:.1f}s, RSS {rss_vals[-1]:.0f} MB", flush=True)

        # scaling: fixed scope=2 substitutions, every engine x thread count
        for eng in ("seqtm", "seqtrie"):
            for t in args.threads:
                q, _ = qpms(idx, queries, st.SearchParams(max_subs=2, max_total_edits=2, engine=eng), t)
                scaling[(eng, t)].append(round(q, 3))
        # matrix-scoring overhead: seqtm scope=3, unit vs BLOSUM62 vs PAM50
        for m in ("unit", "BLOSUM62", "PAM50"):
            kw = dict(max_subs=3, engine="seqtm")
            if m != "unit":
                kw["matrix"] = m
            q, _ = qpms(idx, queries, st.SearchParams(**kw), tmax)
            mat_over[m].append(round(q, 3))
        # seqtrie C++ global-alignment CIGAR fetch cost
        ap_ = st.SearchParams(matrix="BLOSUM62", max_penalty=20, gap_open=8, engine="seqtrie")
        res = idx.search_batch(queries, ap_, threads=tmax)
        pairs = [(h.ref_id, qy) for qy, hits in zip(queries, res) for h in hits][:2000]
        t0 = time.perf_counter()
        for rid, qy in pairs:
            idx.align(rid, qy, ap_).ops
        align_vals.append(round((time.perf_counter() - t0) * 1e6 / max(1, len(pairs)), 3))
        print(f"# [{n:,}] align fetch {align_vals[-1]:.2f} us/call ({len(pairs)} calls)", flush=True)

        # scope + selectivity sweeps only at the largest size (heaviest, most informative)
        if si == len(sizes) - 1:
            sm, sr = [], []
            sel_sm, sel_sr = [], []
            for e in edits:
                q1, h1 = qpms(idx, queries, st.SearchParams(max_subs=e, max_total_edits=e, engine="seqtm"), tmax)
                q2, h2 = qpms(idx, queries, st.SearchParams(max_total_edits=e, engine="seqtrie"), tmax)
                sm.append(round(q1, 3)); sr.append(round(q2, 3))
                sel_sm.append(round(h1, 2)); sel_sr.append(round(h2, 2))
                print(f"# [{n:,}] scope e={e}: seqtm {q1:.2f} q/ms ({h1:.1f} h/q), "
                      f"seqtrie {q2:.2f} q/ms ({h2:.1f} h/q)", flush=True)
            scope_data = (sm, sr)
            sel_scope = (sel_sm, sel_sr)
            bl, pm = [], []
            for p in penalties:
                _, hb = qpms(idx, queries, st.SearchParams(matrix="BLOSUM62", max_penalty=p, gap_open=8, engine="seqtrie"), tmax)
                _, hp = qpms(idx, queries, st.SearchParams(matrix="PAM50", max_penalty=p, gap_open=8, engine="seqtrie"), tmax)
                bl.append(round(hb, 2)); pm.append(round(hp, 2))
            sel_score = (bl, pm)

    # --- render ---
    plot(out, "scaling", "Throughput vs reference size (2 substitutions)",
         "reference set size", "queries / ms", sizes,
         [(f"{eng} t={t}", scaling[(eng, t)], eng_style(eng, THREAD_PT[t]))
          for eng in ("seqtm", "seqtrie") for t in args.threads])
    plot(out, "matrix_overhead", "seqtm matrix-scoring cost (3 substitutions)",
         "reference set size", "queries / ms", sizes,
         [("seqtm unit", mat_over["unit"], eng_style("seqtm", 6)),
          ("seqtm BLOSUM62", mat_over["BLOSUM62"], eng_style("seqtm", 4)),
          ("seqtm PAM50", mat_over["PAM50"], eng_style("seqtm", 8))])
    plot(out, "align_fetch", "seqtrie global-alignment CIGAR fetch (C++)",
         "reference set size", "us / align() call", sizes,
         [("align()", align_vals, eng_style("seqtrie", 7))])
    plot(out, "ram", "Peak RSS after index build", "reference set size", "peak RSS (MB)", sizes,
         [("peak RSS", rss_vals, "with linespoints lw 2 lc rgb '#2ca02c' pt 7 ps 0.9")], logy=True)

    if scope_data:
        biggest = sizes[-1]
        xt = [(str(e), e) for e in edits]
        plot(out, "scope", f"Throughput vs edit budget ({biggest:,} refs, t={tmax})",
             "edit budget", "queries / ms", edits,
             [("seqtm (subs)", scope_data[0], eng_style("seqtm")),
              ("seqtrie (edit dist)", scope_data[1], eng_style("seqtrie"))],
             logx=False, logy=True, xtics=xt)
        plot(out, "selectivity_scope", f"Matches per query vs edit budget ({biggest:,} refs)",
             "edit budget", "matches / query", edits,
             [("seqtm (Hamming ball)", sel_scope[0], eng_style("seqtm")),
              ("seqtrie (edit ball)", sel_scope[1], eng_style("seqtrie"))],
             logx=False, logy=True, xtics=xt)
        plot(out, "selectivity_score", f"Matches per query vs score budget ({biggest:,} refs, gap_open=8)",
             "max penalty budget", "matches / query", penalties,
             [("seqtrie BLOSUM62", sel_score[0], eng_style("seqtrie", 7)),
              ("seqtrie PAM50", sel_score[1], eng_style("seqtrie", 4))],
             logx=False, logy=True)
    print(f"\nWrote SVG figures + TSVs to {out}/", flush=True)


if __name__ == "__main__":
    main()
