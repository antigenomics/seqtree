#!/usr/bin/env python3
"""MHC-allele guessing benchmark (class I and class II).

Reverse problem: can the alleles of a peptide's nearest TCR-/sequence-neighbours predict
its own restricting allele? For each held-out peptide we widen the scope until it has
10-100 non-exact homologs, tally their alleles, and test each allele's count against the
background allele frequency (binomial tail). The guessed allele is the most enriched; the
**aggregate E-value** = (#alleles tested) x p_top (Bonferroni over the allele panel), and
confidence = 1 - min(1, E).

A random-peptide arm (length-matched, uniform residues) is the noise control: random
peptides should get few/again-random neighbours -> high E -> filtered out. We report top-1
accuracy, the real-vs-random E-value separation (AUROC), and a fraction-confident-vs-E-value
curve (the figure). Output: TSV table + bench/figures/mhc_guess.svg.

    python bench/bench_mhc_guess.py --pmhc /Users/mikesh/hf/pmhc_data/pmhc_full.tsv.gz

Needs gnuplot + the pmhc_data table (local path or HuggingFace download).
"""
import argparse
import csv
import gzip
import math
import os
import random
from collections import Counter
from pathlib import Path

from seqtree import KmerIndex, SearchParams, layout
from bench_gnuplot import render, style  # shared gnuplot helpers

AA = "ACDEFGHIKLMNPQRSTVWY"
COLOR = {"real": "#d62728", "random": "#1f77b4"}
# scope ladder on the presentation (anchor) signature: max substitutions, strict -> loose
SCOPES = [0, 1, 2, 3]


def binom_sf(k, n, p):
    """P(Binomial(n, p) >= k), exact (n <= ~100 here)."""
    if k <= 0:
        return 1.0
    if p <= 0:
        return 0.0
    if p >= 1:
        return 1.0
    s = 0.0
    for i in range(k, n + 1):
        s += math.comb(n, i) * p**i * (1 - p) ** (n - i)
    return min(1.0, s)


def guess(index, peptide_allele, cls, query, allele_freq, n_alleles, lo=10, hi=100):
    """Guess the restricting allele from presentation-signature neighbours.

    Widen the scope (subs on the anchor signature) until the query has >= lo non-exact
    neighbours; tally their alleles; test the top allele's count vs background (binomial);
    return (guessed_allele, k, n, aggregate_E). `peptide_allele[pid]` maps neighbours to
    alleles. Excludes peptides identical to the query (self)."""
    feats = layout.presentation_features(query, cls)
    cands = []
    for sc in SCOPES:
        p = SearchParams(max_subs=sc, engine="seqtm")
        cands = index.seed_and_gather([feats], p, 1, -1, 1)[0]
        cands = [c for c in cands if peptide_allele[c.peptide_id][1] != query]  # drop self
        if len(cands) >= lo:
            break
    if not cands:
        return None
    cands = cands[:hi]
    cnt = Counter(peptide_allele[c.peptide_id][0] for c in cands)
    n = sum(cnt.values())
    best = None
    for allele, k in cnt.items():
        f = allele_freq.get(allele, 1.0 / max(1, n_alleles))
        pv = binom_sf(k, n, f)
        if best is None or pv < best[3]:
            best = (allele, k, n, pv)
    allele, k, n, pv = best
    E = min(n_alleles * pv, float(n_alleles))  # Bonferroni aggregate E-value over the panel
    return allele, k, n, E


def auroc(scores_pos, scores_neg):
    """AUROC with `scores` higher = more confident (here: -log10 E)."""
    labelled = [(s, 1) for s in scores_pos] + [(s, 0) for s in scores_neg]
    labelled.sort(key=lambda x: x[0])
    npos = len(scores_pos)
    nneg = len(scores_neg)
    if npos == 0 or nneg == 0:
        return float("nan")
    rank_sum = 0
    i = 0
    while i < len(labelled):
        j = i
        while j < len(labelled) and labelled[j][0] == labelled[i][0]:
            j += 1
        avg_rank = (i + 1 + j) / 2.0
        for t in range(i, j):
            if labelled[t][1] == 1:
                rank_sum += avg_rank
        i = j
    return (rank_sum - npos * (npos + 1) / 2.0) / (npos * nneg)


def load_rows(path, cls_key):
    csv.field_size_limit(10**7)
    op = gzip.open if str(path).endswith(".gz") else open
    want = "MHCI" if cls_key == "mhc1" else "MHCII"
    out = []
    with op(path, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if str(r.get("mhc_class")) == want and r.get("epitope") and r.get("mhc_a"):
                out.append((r["epitope"].strip().upper(), r["mhc_a"].strip()))
    return out


def run_class(cls, rows, n_queries, rng):
    # keep alleles with a solid background; dedup peptides to unique clonotypes
    by_allele = Counter(a for _, a in rows)
    keep = {a for a, c in by_allele.items() if c >= 200}
    rows = [(e, a) for e, a in rows if a in keep and all(c in AA for c in e)]
    rows = list(dict.fromkeys(rows))  # unique (epitope, allele)
    if not rows:
        return None
    total = len(rows)
    n_keep = len(keep)
    allele_freq = {a: by_allele[a] / total for a in keep}
    lengths = [len(e) for e, _ in rows]

    # presentation-feature index: anchors kept, TCR-facing dropped
    allele_to_id = {}
    allele_ids = []
    peptide_allele = []  # pid -> (allele, epitope)
    feats = []
    for e, a in rows:
        allele_to_id.setdefault(a, len(allele_to_id))
        allele_ids.append(allele_to_id[a])
        peptide_allele.append((a, e))
        feats.append(layout.presentation_features(e, cls))
    index = KmerIndex.build(feats, alphabet="aa", allele_ids=allele_ids)

    test = rng.sample(rows, min(n_queries, len(rows)))
    real_E, correct, n_eval = [], 0, 0
    conf_correct = conf_n = 0  # accuracy among confident calls (aggregate E < 1)
    for ep, true_allele in test:
        g = guess(index, peptide_allele, cls, ep, allele_freq, n_keep)
        if g is None:
            continue
        n_eval += 1
        real_E.append(g[3])
        ok = g[0] == true_allele
        correct += ok
        if g[3] < 1.0:
            conf_n += 1
            conf_correct += ok
    # random arm: length-matched uniform peptides (the noise control)
    rand_E = []
    for _ in range(len(test)):
        L = rng.choice(lengths)
        rp = "".join(rng.choice(AA) for _ in range(L))
        g = guess(index, peptide_allele, cls, rp, allele_freq, n_keep)
        if g is not None:
            rand_E.append(g[3])

    def nlog(es):
        return [-math.log10(max(e, 1e-12)) for e in es]

    acc = correct / n_eval if n_eval else 0.0
    acc_conf = conf_correct / conf_n if conf_n else float("nan")
    roc = auroc(nlog(real_E), nlog(rand_E))
    return {"store_n": index.num_peptides(), "alleles": n_keep, "n_eval": n_eval,
            "acc": acc, "acc_conf": acc_conf, "frac_conf": conf_n / n_eval if n_eval else 0.0,
            "auroc": roc, "real_E": real_E, "rand_E": rand_E,
            "mean_real_E": sum(real_E) / len(real_E) if real_E else float("nan"),
            "mean_rand_E": sum(rand_E) / len(rand_E) if rand_E else float("nan")}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pmhc", default="/Users/mikesh/hf/pmhc_data/pmhc_full.tsv.gz")
    ap.add_argument("--queries", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="bench/figures")
    args = ap.parse_args()
    if not os.path.exists(args.pmhc):
        raise SystemExit(f"pmhc table not found: {args.pmhc}")
    rng = random.Random(args.seed)

    thr = [10.0, 1.0, 0.1, 0.01, 0.001]  # E-value thresholds for the confident-fraction curve
    panels, summary = [], {}
    for cls, label in (("mhc1", "MHC class I"), ("mhc2", "MHC class II")):
        res = run_class(cls, load_rows(args.pmhc, cls), args.queries, rng)
        if res is None:
            print(f"# {cls}: no data"); continue
        summary[cls] = res
        # fraction with aggregate E <= threshold (confident), real vs random
        def frac_at(es):
            return [sum(1 for e in es if e <= t) / max(1, len(es)) for t in thr]
        panels.append({
            "title": f"{label}: allele-guess confidence (E<=thr), real vs random  "
                     f"[top-1 acc {res['acc']:.2f}, AUROC {res['auroc']:.2f}]",
            "xlabel": "aggregate E-value threshold", "ylabel": "fraction confident",
            "xs": thr, "logx": True,
            "series": [("real (presented)", frac_at(res["real_E"]), style(COLOR["real"], "seqtm")),
                       ("random peptides", frac_at(res["rand_E"]), style(COLOR["random"], "seqtrie"))]})

    print("class\tstore_n\talleles\tn_eval\ttop1_acc\tacc@E<1\tfrac@E<1\tAUROC\tmeanE_real\tmeanE_rand")
    for cls, r in summary.items():
        print(f"{cls}\t{r['store_n']}\t{r['alleles']}\t{r['n_eval']}\t{r['acc']:.3f}\t"
              f"{r['acc_conf']:.3f}\t{r['frac_conf']:.3f}\t{r['auroc']:.3f}\t"
              f"{r['mean_real_E']:.3g}\t{r['mean_rand_E']:.3g}")
    if panels:
        out = Path(args.out)
        out.mkdir(parents=True, exist_ok=True)
        render(out, "mhc_guess", panels)
        print(f"\nWrote {out}/mhc_guess.svg")


if __name__ == "__main__":
    main()
