#!/usr/bin/env python3
"""Epitope-specific detection complexity: A*02 NLV (diverse) vs GIL (one big cluster).

Tests the detectability formalism of the E-value appendix on two HLA-A*02 epitopes
with opposite repertoire structure:

  * NLVPMVATV (CMV pp65)   -- many small, diverse TCR clusters
  * GILGFVFTL (influenza M1) -- one dominant convergent (public) cluster

For each epitope we (1) characterize the within-epitope neighbour graph at a scope
(neighbour density rho, degree distribution, connected-component / cluster sizes),
and (2) subsample the repertoire to depth n, scoring each sampled node's enrichment
against an OLGA background control, and report the fraction of the sample called
significant (Benjamini-Hochberg FDR < 0.05). The predicted detection curve from the
cluster-size distribution is overlaid. Output: a summary table + a stacked SVG.

    python bench/bench_epitope.py                # scope 1 and 2, OLGA 1M control
    python bench/bench_epitope.py --control bench/cache/olga_2M.txt.gz --scopes 1 2 3

Needs gnuplot, a cached OLGA control (bench/gen_olga.py), and the VDJdb dataset.
"""
import argparse
import csv
import gzip
import os
import random
import subprocess
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import seqtree as st
from seqtree.evalue import evalues

AA = set("ACDEFGHIKLMNPQRSTVWY")
EPITOPES = {"NLV (NLVPMVATV)": "NLVPMVATV", "GIL (GILGFVFTL)": "GILGFVFTL"}
COLOR = {"NLV (NLVPMVATV)": "#1f77b4", "GIL (GILGFVFTL)": "#d62728"}


def load_epitope_trb(epi):
    from huggingface_hub import hf_hub_download

    csv.field_size_limit(10**7)
    path = hf_hub_download("isalgo/airr_benchmark",
                           "vdjdb/vdjdb-2025-12-29/vdjdb.slim.txt.gz", repo_type="dataset")
    seqs = set()
    with gzip.open(path, "rt") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["gene"] == "TRB" and row["antigen.epitope"] == epi and "A*02" in row.get("mhc.a", ""):
                s = row["cdr3"].strip().upper()
                if s and all(c in AA for c in s):
                    seqs.add(s)
    return sorted(seqs)


def load_cache(path):
    with gzip.open(path, "rt") as fh:
        return [ln.strip() for ln in fh if ln.strip()]


def neighbour_graph(members, scope):
    """Return (degree[], cluster_sizes_desc) for the within-set scope-neighbour graph."""
    idx = st.Index.build(members, alphabet="aa")
    res = idx.search_batch(members, st.SearchParams(max_subs=scope, engine="seqtm"), 0)
    parent = list(range(len(members)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    deg = [0] * len(members)
    for i, hits in enumerate(res):
        for h in hits:
            if h.ref_id != i:
                deg[i] += 1
                ri, rj = find(i), find(h.ref_id)
                if ri != rj:
                    parent[ri] = rj
    sizes = sorted(Counter(find(i) for i in range(len(members))).values(), reverse=True)
    return deg, sizes


def detect_fraction(members, control, scope, n, rng, alpha=0.05):
    """Subsample n members, return fraction called significant (BH FDR<alpha) vs control."""
    sample = rng.sample(members, min(n, len(members)))
    sub = st.Index.build(sample, alphabet="aa")
    p = st.SearchParams(max_subs=scope, engine="seqtm")
    ev = evalues(sub, control, sample, p, exclude_exact=True)
    pvals = [e["p_enrichment"] for e in ev]
    # Benjamini-Hochberg
    order = sorted(range(len(pvals)), key=lambda i: pvals[i])
    kmax = 0
    for rank, idx in enumerate(order, 1):
        if pvals[idx] <= rank / len(pvals) * alpha:
            kmax = rank
    return kmax / len(sample) if sample else 0.0


def detect_fraction_ci(members, control, scope, n, reps, seed):
    """Bootstrap the detection fraction over independent depth-n subsamples; return
    (mean, lo, hi) with [lo, hi] the 95% percentile interval across `reps` resamples."""
    fr = sorted(detect_fraction(members, control, scope, n, random.Random(seed + r)) for r in range(reps))
    mean = sum(fr) / len(fr)
    return mean, fr[int(0.025 * (reps - 1))], fr[int(round(0.975 * (reps - 1)))]


def render_ci(out, key, panels, width=720, panel_h=360):
    """Like bench_gnuplot.render but draws a 95% CI band (filledcurves) under each
    'obs' series. Each series: dict(label, mean[], lo[], hi[], color, kind in {obs,pred})."""
    n = len(panels)
    L = [f"set terminal svg size {width},{panel_h * n} font 'Helvetica,12' background rgb 'white'",
         f"set output '{key}.svg'", "set datafile separator '\\t'", "set grid",
         "set key outside right top"]
    if n > 1:
        L.append(f"set multiplot layout {n},1")
    for pi, p in enumerate(panels):
        tsv = f"{key}_{pi}.tsv"
        rows = []
        for i, x in enumerate(p["xs"]):
            cells = [f"{x:g}"]
            for s in p["series"]:
                cells += [f"{s['mean'][i]:g}", f"{s['lo'][i]:g}", f"{s['hi'][i]:g}"]
            rows.append("\t".join(cells))
        (out / tsv).write_text("x\t" + "\t".join(f"m{i}\tlo{i}\thi{i}" for i in range(len(p["series"]))) +
                               "\n" + "\n".join(rows) + "\n")
        L += [f'set title "{p["title"]}"', f"set xlabel '{p['xlabel']}'", f"set ylabel '{p['ylabel']}'",
              "set logscale x" if p.get("logx") else "unset logscale x",
              f"set yrange [{p['yrange'][0]}:{p['yrange'][1]}]" if p.get("yrange") else "set yrange [*:*]"]
        plots = []
        for si, s in enumerate(p["series"]):
            yc = 2 + 3 * si
            if s["kind"] == "obs":
                plots.append(f"'{tsv}' using 1:{yc + 1}:{yc + 2} with filledcurves "
                             f"lc rgb '{s['color']}' fs transparent solid 0.15 notitle")
                plots.append(f"'{tsv}' using 1:{yc} with linespoints lw 2 lc rgb '{s['color']}' "
                             f"dt 2 pt 7 ps 0.6 title '{s['label']}'")
            else:  # pred: deterministic dotted line, no band
                plots.append(f"'{tsv}' using 1:{yc} with lines lw 2 lc rgb '{s['color']}' dt 3 notitle")
        L.append("plot " + ", ".join(plots))
    if n > 1:
        L.append("unset multiplot")
    (out / f"{key}.gp").write_text("\n".join(L) + "\n")
    subprocess.run(["gnuplot", f"{key}.gp"], cwd=out, check=True)


def predicted_fraction(deg, K, n, d_min=1.0):
    """Fraction of a depth-n sample expected significant from the degree distribution:
    a node of full-set degree d retains ~ d*(n-1)/(K-1) of its neighbours when n are
    sampled, so it is detectable (>= d_min co-sampled neighbours over a ~0 background)
    when d >= d_min*(K-1)/(n-1). The epitope's degree law thus sets the detection curve."""
    if n <= 1:
        return 0.0
    thr = d_min * (K - 1) / (n - 1)
    return sum(1 for d in deg if d >= thr) / K


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--control", default="bench/cache/olga_1M.txt.gz")
    ap.add_argument("--scopes", type=int, nargs="*", default=[1, 2])
    ap.add_argument("--reps", type=int, default=30, help="bootstrap resamples for the 95% CI band")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="bench/figures")
    args = ap.parse_args()

    rng = random.Random(args.seed)
    if not os.path.exists(args.control):
        raise SystemExit(f"control cache {args.control} missing -- run bench/gen_olga.py first")
    control = st.Index.build(load_cache(args.control), alphabet="aa")
    M = len(control)
    sets = {name: load_epitope_trb(epi) for name, epi in EPITOPES.items()}
    print(f"# control (OLGA) M={M:,}; " + "; ".join(f"{n}: K={len(s)}" for n, s in sets.items()))

    # one shared depth grid (<= every epitope's K) so all series align on the x-axis
    kmin = min(len(s) for s in sets.values())
    depths = [d for d in (50, 100, 200, 500, 1000, 2000, 5000) if d <= kmin]

    print("\nepitope\tscope\tK\trho\tmean_deg\tmax_deg\tn_clusters\tlargest\tsingleton_frac")
    panels = []
    for scope in args.scopes:
        series = []
        for name, members in sets.items():
            K = len(members)
            deg, sizes = neighbour_graph(members, scope)
            edges = sum(deg) // 2
            rho = edges / (K * (K - 1) / 2) if K > 1 else 0.0
            sing_frac = sum(1 for s in sizes if s == 1) / len(sizes)
            print(f"{name}\t{scope}\t{K}\t{rho:.2e}\t{sum(deg)/K:.2f}\t{max(deg)}\t"
                  f"{len(sizes)}\t{sizes[0]}\t{sing_frac:.3f}")
            means, los, his = [], [], []
            for d in depths:
                m, lo, hi = detect_fraction_ci(members, control, scope, d, args.reps, args.seed + d)
                means.append(m); los.append(lo); his.append(hi)
            pred = [predicted_fraction(deg, K, d) for d in depths]
            # observed = dashed line + 95% CI band; predicted = dotted (same colour per epitope);
            # legend shows only the two epitopes (NLV / GIL), predicted lines are unlabelled.
            short, c = name.split()[0], COLOR[name]
            series.append({"label": short, "mean": means, "lo": los, "hi": his, "color": c, "kind": "obs"})
            series.append({"label": "", "mean": pred, "lo": pred, "hi": pred, "color": c, "kind": "pred"})
        panels.append({"title": f"epitope detection vs sampling depth (scope {scope}, BH FDR<0.05; "
                                f"band = 95% CI, {args.reps} resamples)",
                       "xlabel": "sampled TCRs (depth n)", "ylabel": "fraction significant",
                       "xs": depths, "logx": True, "yrange": (0, 1), "series": series})

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    render_ci(out, "epitope_detection", panels)
    print(f"\nWrote {out}/epitope_detection.svg")


if __name__ == "__main__":
    main()
