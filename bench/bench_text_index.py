#!/usr/bin/env python
"""TextIndex: one index, every query length.

The question this answers is the one that motivated the class. `Index` builds a trie over
reference *strings*, so a text question means enumerating every length-L window as its own
string -- and again for every distinct query length. On a real corpus (445,466 peptides across
**45 distinct lengths**) that ran > 2 h 10 m without finishing and filled 225 GB of index
cache. `TextIndex` keys on a k-mer seed table instead, so `k` belongs to the index: one build
serves every length.

Three tables:

  A  build cost and index size vs text size and k
  B  ms/query and hits/query across the query-length x max_subs grid, from ONE index
  C  thread scaling

and, on every configuration in A and B, **recall against a brute-force scan** -- the number
that decides whether any of the rest is worth reading. A missed hit does not degrade the
caller's answer, it changes it, so recall below 1.000 anywhere means the work has failed.

Texts, in size steps. The small ones are generated so the script runs with no downloads; the
real proteomes come from the HuggingFace dataset `isalgo/pmhc_data` (see SOURCES.md) and are
used only when `huggingface_hub` is installed.

Usage:
    python bench/bench_text_index.py                    # fast tier, synthetic text
    RUN_BENCHMARK=1 python bench/bench_text_index.py     # adds mouse + human proteomes

2026-09-06
"""
from __future__ import annotations

import gzip
import os
import random
import resource
import sys
import time

from seqtree import TextIndex

RUN_BENCHMARK = bool(os.getenv("RUN_BENCHMARK"))
AA = "ACDEFGHIKLMNPQRSTVWY"


def peak_rss_mb() -> float:
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024


def synthetic(n_records: int, mean_len: int, seed: int = 20260906) -> list[str]:
    rng = random.Random(seed)
    return ["".join(rng.choice(AA) for _ in range(max(20, int(rng.gauss(mean_len, mean_len / 3)))))
            for _ in range(n_records)]


def from_fasta_gz(path: str) -> list[str]:
    """Minimal FASTA reader -- seqtree has no runtime dependencies and this is a bench script."""
    seqs, cur = [], []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return seqs


def proteome(name: str) -> list[str] | None:
    """`isalgo/pmhc_data` proteomes, if huggingface_hub is available. See SOURCES.md."""
    try:
        from huggingface_hub import hf_hub_download
    except ImportError:
        return None
    try:
        path = hf_hub_download("isalgo/pmhc_data", f"proteome/{name}.fasta.gz", repo_type="dataset")
    except Exception as exc:  # network, auth, or a renamed file -- all equally "skip this tier"
        print(f"  (skipping {name}: {exc})", file=sys.stderr)
        return None
    return from_fasta_gz(path)


def brute_force(refs: list[str], q: str, max_subs: int) -> set[tuple[int, int, int]]:
    out = set()
    for r, ref in enumerate(refs):
        for off in range(len(ref) - len(q) + 1):
            n = sum(a != b for a, b in zip(q, ref[off:off + len(q)]))
            if n <= max_subs:
                out.add((r, off, n))
    return out


def found(res, i) -> set[tuple[int, int, int]]:
    a = res.to_numpy()
    lo, hi = a["query_begin"][i], a["query_begin"][i + 1]
    return {(int(a["ref_id"][h]), int(a["offset"][h]), int(a["n_subs"][h])) for h in range(lo, hi)}


def sample_queries(refs: list[str], length: int, n: int, seed: int) -> list[str]:
    """Real windows of the text, so the measurement is of hits and not of empty lookups."""
    rng = random.Random(seed)
    long_enough = [r for r in refs if len(r) >= length]
    out = []
    for _ in range(n):
        r = rng.choice(long_enough)
        off = rng.randrange(len(r) - length + 1)
        out.append(r[off:off + length])
    return out


def recall(refs, ix, length: int, max_subs: int, n: int = 12) -> float:
    """Fraction of brute-force hits the index returns. Anything below 1.0 is a failure."""
    queries = sample_queries(refs, length, n, seed=7)
    res = ix.search_batch(queries, max_subs=max_subs, threads=1)
    want = total = 0
    for i, q in enumerate(queries):
        truth = brute_force(refs, q, max_subs)
        total += len(truth)
        want += len(truth & found(res, i))
    return want / total if total else 1.0


def table_a(texts: dict[str, list[str]]) -> None:
    print("\n## A. Build cost and index size\n")
    print("| text | records | residues | k | build s | peak RSS MB |")
    print("|---|--:|--:|--:|--:|--:|")
    for name, refs in texts.items():
        residues = sum(len(r) for r in refs)
        for k in (3, 4, 5):
            t0 = time.perf_counter()
            ix = TextIndex.build(refs, alphabet="aa", k=k)
            dt = time.perf_counter() - t0
            print(f"| {name} | {len(refs):,} | {residues:,} | {k} | {dt:.2f} | {peak_rss_mb():.0f} |")
            del ix


def table_b(name: str, refs: list[str], check_recall: bool) -> None:
    print(f"\n## B. One k=4 index, every query length -- {name}\n")
    ix = TextIndex.build(refs, alphabet="aa", k=4)
    print("| L | max_subs | ms/query | hits/query | recall vs brute force |")
    print("|--:|--:|--:|--:|--:|")
    for length in (8, 9, 10, 11, 12, 15, 20, 25, 31, 50):
        if not any(len(r) >= length for r in refs):
            continue
        for max_subs in (0, 1, 2, 3):
            # The ball path is entered when L / (max_subs + 1) < k, and its cost is
            # sum_i<=m C(k,i)(A-1)^i -- 3,267 probes at m = 2 but 51,935 at m = 3, which on a
            # 69.6 M-residue text touches a sixth of every posting list. Fewer queries there
            # keeps the cell bounded; ms/query does not depend on how many were timed.
            heavy = length // (max_subs + 1) < 4 and len(refs) > 5_000
            n_q = 25 if heavy else (200 if RUN_BENCHMARK else 50)
            queries = sample_queries(refs, length, n_q, seed=length)
            t0 = time.perf_counter()
            res = ix.search_batch(queries, max_subs=max_subs, threads=1)
            ms = (time.perf_counter() - t0) * 1000 / len(queries)
            hits = res.num_hits / len(queries)
            rec = f"{recall(refs, ix, length, max_subs):.3f}" if check_recall else "--"
            print(f"| {length} | {max_subs} | {ms:.3f} | {hits:.1f} | {rec} |")


def table_c(name: str, refs: list[str]) -> None:
    print(f"\n## C. Thread scaling -- {name}, L=12, max_subs=2\n")
    ix = TextIndex.build(refs, alphabet="aa", k=4)
    queries = sample_queries(refs, 12, 5000 if RUN_BENCHMARK else 500, seed=3)
    print("| threads | s | ms/query | speedup |")
    print("|--:|--:|--:|--:|")
    base = None
    for threads in (1, 2, 4, 8, 0):
        t0 = time.perf_counter()
        ix.search_batch(queries, max_subs=2, threads=threads)
        dt = time.perf_counter() - t0
        base = base or dt
        label = "all" if threads == 0 else str(threads)
        print(f"| {label} | {dt:.2f} | {dt * 1000 / len(queries):.3f} | {base / dt:.1f}x |")


def main() -> None:
    small = synthetic(400, 300)          # ~120 k residues, brute-forceable
    texts = {"synthetic 120k": small}
    if RUN_BENCHMARK:
        for name in ("mouse", "human"):
            refs = proteome(name)
            if refs:
                texts[name] = refs

    print("# TextIndex benchmark")
    print(f"\nRUN_BENCHMARK={'1' if RUN_BENCHMARK else '0'}; "
          f"texts: {', '.join(f'{k} ({sum(len(r) for r in v):,} residues)' for k, v in texts.items())}")

    table_a(texts)
    # Recall is only checked where a brute-force scan is affordable; it is the same code path
    # on every text, and the C++ suite proves set equality exhaustively on small ones.
    for name, refs in texts.items():
        table_b(name, refs, check_recall=(name == "synthetic 120k"))
    for name, refs in texts.items():
        table_c(name, refs)


if __name__ == "__main__":
    main()
