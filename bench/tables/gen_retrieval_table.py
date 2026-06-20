#!/usr/bin/env python3
"""Produce the retrieval precision/recall table — deterministic, no external data.

A reference set of seeded-random amino-acid sequences is indexed; a held-out query
set is searched with **BLOSUM62-weighted** ``seqtrie`` over a penalty budget. Ground
truth is brute-force **Hamming distance <= D** (computed independently of seqtree).
Sweeping the penalty threshold yields a precision/recall/F1 table — a classic
retrieval curve, exactly reproducible because every count is integer-derived.

This is a *table producer*: it runs the compiled seqtree from the current repo and
writes a TSV to stdout (or ``--out``). Plotting is a separate step
(``bench/plots/plot_retrieval_pr.py``); the committed TSV doubles as the CI oracle.

  python bench/tables/gen_retrieval_table.py            # TSV to stdout
  python bench/tables/gen_retrieval_table.py --out bench/tables/retrieval_pr.tsv
"""
import argparse
import sys

import seqtree as st

AA = "ACDEFGHIKLMNPQRSTVWY"
# Fixed parameters define the oracle; changing them is a deliberate oracle update.
SEED, N_REFS, N_QUERIES, LENGTH, TRUTH_D, MAX_PENALTY = 0, 2000, 300, 14, 2, 40


def lcg_pool(n, length, seed):
    """Self-contained deterministic sequence pool (no dependence on the host RNG
    implementation, so the oracle is identical on every platform)."""
    state = (seed * 2862933555777941757 + 3037000493) & (2**64 - 1)
    seqs = []
    for _ in range(n):
        cs = []
        for _ in range(length):
            state = (state * 6364136223846793005 + 1442695040888963407) & (2**64 - 1)
            cs.append(AA[(state >> 33) % len(AA)])
        seqs.append("".join(cs))
    return seqs


def mutate(seq, n_subs, state):
    cs = list(seq)
    for _ in range(n_subs):
        state = (state * 6364136223846793005 + 1442695040888963407) & (2**64 - 1)
        i = (state >> 33) % len(cs)
        state = (state * 6364136223846793005 + 1442695040888963407) & (2**64 - 1)
        cs[i] = AA[(state >> 33) % len(AA)]
    return "".join(cs), state


def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="-", help="output TSV path ('-' = stdout)")
    args = ap.parse_args()

    refs = lcg_pool(N_REFS, LENGTH, SEED)
    # Queries: distinct refs mutated by 1..3 substitutions (deterministic walk).
    state = (SEED * 6364136223846793005 + 1) & (2**64 - 1)
    queries = []
    for k in range(N_QUERIES):
        q, state = mutate(refs[k % N_REFS], 1 + (k % 3), state)
        queries.append(q)

    # Ground truth: brute-force Hamming <= D over query x ref (independent of seqtree).
    truth = [set(j for j, r in enumerate(refs) if hamming(q, r) <= TRUTH_D) for q in queries]
    total_pos = sum(len(t) for t in truth)

    # seqtree candidates: BLOSUM62-weighted seqtrie over the penalty budget.
    idx = st.Index.build(refs, alphabet="aa")
    p = st.SearchParams(matrix="BLOSUM62", max_penalty=MAX_PENALTY, gap_open=8, engine="seqtrie")
    results = idx.search_batch(queries, p, threads=1)

    # Bucket candidates by integer penalty: pos[k] / neg[k] = #candidates at penalty k
    # that are / are not within the Hamming-truth ball.
    pos = [0] * (MAX_PENALTY + 1)
    neg = [0] * (MAX_PENALTY + 1)
    for qi, hits in enumerate(results):
        for h in hits:
            k = int(h.score)
            if 0 <= k <= MAX_PENALTY:
                (pos if h.ref_id in truth[qi] else neg)[k] += 1

    out = sys.stdout if args.out == "-" else open(args.out, "w")
    out.write(f"# retrieval PR: seed={SEED} n_refs={N_REFS} n_queries={N_QUERIES} "
              f"length={LENGTH} truth_hamming<={TRUTH_D} matrix=BLOSUM62 gap_open=8 "
              f"total_positives={total_pos}\n")
    out.write("penalty\ttp\tfp\tfn\tprecision\trecall\tf1\n")
    ctp = cfp = 0
    for k in range(MAX_PENALTY + 1):
        ctp += pos[k]
        cfp += neg[k]
        fn = total_pos - ctp
        prec = ctp / (ctp + cfp) if (ctp + cfp) else 1.0
        rec = ctp / total_pos if total_pos else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        out.write(f"{k}\t{ctp}\t{cfp}\t{fn}\t{prec:.6f}\t{rec:.6f}\t{f1:.6f}\n")
    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
