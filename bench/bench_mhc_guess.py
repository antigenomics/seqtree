#!/usr/bin/env python3
"""MHC-allele guessing benchmark with ROC + PR curves, split by species.

Reverse problem (peptide -> presenting allele), evaluated like vdjmatch evaluates
TCR-antigen specificity: a per-(peptide, allele) binary task. For each held-out peptide
we widen the scope on its *presentation* (anchor) signature until it has 10-100 non-exact
neighbours, vote the neighbours' alleles, and score each panel allele by the enrichment of
its vote vs the background frequency (-log10 binomial tail). Positives = the peptide's true
allele(s) (multi-label, since peptides can be promiscuous); negatives = the rest. We report
ROC-AUC and PR-AUC for MHC-I and MHC-II, **human and mouse separately**, plus top-1 accuracy
and a real-vs-random noise-rejection AUROC. Class-II uses the register trick (single best
9-mer core register, layout.presentation_features(register='anchored')).

    python bench/bench_mhc_guess.py --pmhc /Users/mikesh/hf/pmhc_data/pmhc_full.tsv.gz

Targets (aspirational; tuned later in the mhcmatch package, as vdjmatch tunes TCR specificity):
ROC-AUC ~0.90 (MHC-I), ~0.80 (MHC-II). Output: TSV table + bench/figures/mhc1_rocpr.svg,
mhc2_rocpr.svg. Needs gnuplot + the pmhc_data table.
"""
import argparse
import csv
import gzip
import math
import os
import random
from collections import Counter, defaultdict
from pathlib import Path

from seqtree import KmerIndex, SearchParams, layout

AA = "ACDEFGHIKLMNPQRSTVWY"
SCOPES = [0, 1, 2, 3]
SPECIES = {"HomoSapiens": "human", "MusMusculus": "mouse"}
# frontend-design: one hue per species, dashed grey references, white background.
COLOR = {"human": "#2563eb", "mouse": "#f59e0b"}
MIN_ALLELE = {"human": 100, "mouse": 40}


def binom_sf(k, n, p):
    if k <= 0:
        return 1.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    return min(1.0, sum(math.comb(n, i) * p**i * (1 - p) ** (n - i) for i in range(k, n + 1)))


def auc(xs, ys):
    return sum((xs[i] - xs[i - 1]) * (ys[i] + ys[i - 1]) / 2 for i in range(1, len(xs)))


def roc_pr(scored):
    """scored: list of (score, label in {0,1}); higher score = more confident.
    Returns (roc_xs, roc_ys, roc_auc, pr_xs, pr_ys, pr_auc)."""
    scored = sorted(scored, key=lambda x: -x[0])
    P = sum(l for _, l in scored)
    N = len(scored) - P
    if P == 0 or N == 0:
        return [0, 1], [0, 1], 0.5, [0, 1], [P / max(1, P + N)] * 2, P / max(1, P + N)
    tp = fp = 0
    roc_x, roc_y, pr_r, pr_p = [0.0], [0.0], [0.0], [1.0]
    i = 0
    while i < len(scored):
        j = i
        while j < len(scored) and scored[j][0] == scored[i][0]:
            j += 1
        for t in range(i, j):
            if scored[t][1]:
                tp += 1
            else:
                fp += 1
        roc_x.append(fp / N)
        roc_y.append(tp / P)
        pr_r.append(tp / P)
        pr_p.append(tp / (tp + fp) if tp + fp else 1.0)
        i = j
    return roc_x, roc_y, auc(roc_x, roc_y), pr_r, pr_p, auc(pr_r, pr_p)


def load_rows(path):
    csv.field_size_limit(10**7)
    op = gzip.open if str(path).endswith(".gz") else open
    out = []
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            cls = {"MHCI": "mhc1", "MHCII": "mhc2"}.get(str(r.get("mhc_class")))
            sp = SPECIES.get(str(r.get("mhc_species")))
            ep, a = str(r.get("epitope", "")).strip().upper(), str(r.get("mhc_a", "")).strip()
            if cls and sp and ep and a and all(c in AA for c in ep):
                out.append((ep, a, cls, sp))
    return out


def eval_panel(cls, species, rows, n_queries, rng):
    by_allele = Counter(a for _, a in rows)
    keep = {a for a, c in by_allele.items() if c >= MIN_ALLELE[species]}
    rows = [(e, a) for e, a in rows if a in keep]
    rows = list(dict.fromkeys(rows))
    if len(keep) < 2 or len(rows) < 200:
        return None
    panel = sorted(keep)
    total = len(rows)
    allele_freq = {a: by_allele[a] / total for a in keep}
    pep_alleles = defaultdict(set)
    for e, a in rows:
        pep_alleles[e].add(a)

    allele_to_id, allele_ids, pid_meta, feats = {}, [], [], []
    for e, a in rows:
        allele_to_id.setdefault(a, len(allele_to_id))
        allele_ids.append(allele_to_id[a])
        pid_meta.append((a, e))
        feats.append(layout.presentation_features(e, cls, register="anchored"))
    index = KmerIndex.build(feats, alphabet="aa", allele_ids=allele_ids)

    test_eps = list(dict.fromkeys(e for e, _ in rows))
    rng.shuffle(test_eps)
    test_eps = test_eps[:n_queries]
    lengths = [len(e) for e, _ in rows]

    scored, top1_ok, n_eval, real_conf = [], 0, 0, []
    for ep in test_eps:
        tally = guess_tally(index, pid_meta, cls, ep)
        if tally is None:
            continue
        n_eval += 1
        n = sum(tally.values())
        # ranking score = neighbour vote fraction (posterior P(allele | neighbours)); robust to
        # panel skew. confidence = enrichment vs background (-log10 binomial tail) -> noise rejection.
        vote = {a: k / n for a, k in tally.items()}
        enr = {a: -math.log10(max(binom_sf(k, n, allele_freq[a]), 1e-300)) for a, k in tally.items()}
        truth = pep_alleles[ep]
        for a in panel:
            scored.append((vote.get(a, 0.0), 1 if a in truth else 0))
        best = max(panel, key=lambda a: vote.get(a, 0.0))
        top1_ok += best in truth
        real_conf.append(max(enr.values()) if enr else 0.0)

    rand_conf = []
    for _ in range(len(test_eps)):
        rp = "".join(rng.choice(AA) for _ in range(rng.choice(lengths)))
        tally = guess_tally(index, pid_meta, cls, rp)
        if not tally:
            rand_conf.append(0.0)
            continue
        n = sum(tally.values())
        rand_conf.append(max(-math.log10(max(binom_sf(k, n, allele_freq[a]), 1e-300))
                             for a, k in tally.items()))

    rx, ry, r_auc, px, py, p_auc = roc_pr(scored)
    noise_auc = auc(*_roc_only(real_conf, rand_conf))
    return {"alleles": len(keep), "n_eval": n_eval, "top1": top1_ok / max(1, n_eval),
            "roc_auc": r_auc, "pr_auc": p_auc, "pr_base": sum(l for _, l in scored) / max(1, len(scored)),
            "noise_auc": noise_auc, "roc": (rx, ry), "pr": (px, py)}


def guess_tally(index, pid_meta, cls, query, lo=10, hi=100):
    feats = layout.presentation_features(query, cls, register="anchored")
    cands = []
    for sc in SCOPES:
        p = SearchParams(max_subs=sc, engine="seqtm")
        cands = [c for c in index.seed_and_gather([feats], p, 1, -1, 1)[0]
                 if pid_meta[c.peptide_id][1] != query]
        if len(cands) >= lo:
            break
    if not cands:
        return None
    return Counter(pid_meta[c.peptide_id][0] for c in cands[:hi])


def _roc_only(pos, neg):
    scored = sorted([(s, 1) for s in pos] + [(s, 0) for s in neg], key=lambda x: -x[0])
    P, N = len(pos), len(neg)
    if P == 0 or N == 0:
        return [0, 1], [0, 1]
    tp = fp = 0
    xs, ys = [0.0], [0.0]
    i = 0
    while i < len(scored):
        j = i
        while j < len(scored) and scored[j][0] == scored[i][0]:
            j += 1
        for t in range(i, j):
            tp += scored[t][1]
            fp += 1 - scored[t][1]
        xs.append(fp / N)
        ys.append(tp / P)
        i = j
    return xs, ys


def render_rocpr(out, key, title, curves):
    """curves: {species: result}. Two panels (ROC, PR), one line per species + references."""
    def tsv(name, xs, ys):
        (out / name).write_text("x\ty\n" + "\n".join(f"{x:g}\t{y:g}" for x, y in zip(xs, ys)) + "\n")

    lines = ["set terminal svg size 760,380 font 'Helvetica,12' background rgb 'white'",
             f"set output '{key}.svg'", "set datafile separator '\\t'", "set multiplot layout 1,2",
             "set grid lc rgb '#e5e7eb'", "set key bottom right box lc rgb '#d1d5db'",
             "set size square", "set xrange [0:1]", "set yrange [0:1.02]"]
    # ROC panel
    lines += [f"set title '{title} - ROC'", "set xlabel 'false positive rate'",
              "set ylabel 'true positive rate'"]
    plots = ["x with lines lc rgb '#9ca3af' dt 2 notitle"]
    for sp, r in curves.items():
        tsv(f"{key}_{sp}_roc.tsv", *r["roc"])
        plots.append(f"'{key}_{sp}_roc.tsv' u 1:2 w l lw 2.5 lc rgb '{COLOR[sp]}' "
                     f"title '{sp} (AUC {r['roc_auc']:.2f})'")
    lines.append("plot " + ", ".join(plots))
    # PR panel
    lines += [f"set title '{title} - precision-recall'", "set xlabel 'recall'",
              "set ylabel 'precision'", "set key top right box lc rgb '#d1d5db'"]
    plots = []
    for sp, r in curves.items():
        tsv(f"{key}_{sp}_pr.tsv", *r["pr"])
        plots.append(f"'{key}_{sp}_pr.tsv' u 1:2 w l lw 2.5 lc rgb '{COLOR[sp]}' "
                     f"title '{sp} (AUC {r['pr_auc']:.2f})'")
        plots.append(f"{r['pr_base']:g} w l lc rgb '{COLOR[sp]}' dt 3 notitle")
    lines.append("plot " + ", ".join(plots))
    lines.append("unset multiplot")
    (out / f"{key}.gp").write_text("\n".join(lines) + "\n")
    import subprocess
    subprocess.run(["gnuplot", f"{key}.gp"], cwd=out, check=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default="/Users/mikesh/hf/pmhc_data/pmhc_full.tsv.gz")
    ap.add_argument("--queries", type=int, default=800)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="bench/figures")
    args = ap.parse_args()
    if not os.path.exists(args.pmhc):
        raise SystemExit(f"pmhc table not found: {args.pmhc}")
    rng = random.Random(args.seed)
    rows_all = load_rows(args.pmhc)

    print("class\tspecies\talleles\tn_eval\ttop1\tROC_AUC\tPR_AUC\tPR_base\tnoise_AUROC")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    for cls, label in (("mhc1", "MHC-I"), ("mhc2", "MHC-II")):
        curves = {}
        for sp in ("human", "mouse"):
            rows = [(e, a) for e, a, c, s in rows_all if c == cls and s == sp]
            res = eval_panel(cls, sp, rows, args.queries, rng)
            if res is None:
                print(f"{cls}\t{sp}\t(insufficient data)")
                continue
            curves[sp] = res
            print(f"{cls}\t{sp}\t{res['alleles']}\t{res['n_eval']}\t{res['top1']:.3f}\t"
                  f"{res['roc_auc']:.3f}\t{res['pr_auc']:.3f}\t{res['pr_base']:.3g}\t{res['noise_auc']:.3f}")
        if curves:
            render_rocpr(out, f"{cls}_rocpr", f"{label} allele guessing", curves)
            print(f"# wrote {out}/{cls}_rocpr.svg")


if __name__ == "__main__":
    main()
