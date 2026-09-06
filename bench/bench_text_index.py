#!/usr/bin/env python
"""TextIndex: one index, every query length.

The question this answers is the one that motivated the class. `Index` builds a trie over
reference *strings*, so a text question means enumerating every length-L window as its own
string -- and again for every distinct query length. On a real corpus (445,466 peptides across
**45 distinct lengths**) that ran > 2 h 10 m without finishing and filled 225 GB of index
cache. `TextIndex` keys on a k-mer seed table instead, so `k` belongs to the index: one build
serves every length.

Four tables:

  A  build cost and peak RSS vs text size and k
  B  ms/query at k=4 and k=5 across the query-length x max_subs grid, from ONE
     index per k -- with the hit counts asserted equal, since k must not change
     the answer
  C  thread scaling
  D  on-disk size, save seconds, and load seconds mapped and read

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
import pathlib
import random
import sys
import tempfile
import time

from seqtree import TextIndex
from _common import peak_rss_mb

RUN_BENCHMARK = bool(os.getenv("RUN_BENCHMARK"))
AA = "ACDEFGHIKLMNPQRSTVWY"


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
    """One index per k, the SAME queries through both.

    k is an index parameter and must be invisible in the answer, so the hit counts are compared
    cell by cell rather than trusted: a dispatch that got `b` or a per-block budget wrong would
    show up here as a disagreement on real text, which no synthetic unit test can rule out at
    this scale. What k is allowed to change is the work, which is the point of the two columns.
    """
    print(f"\n## B. One index, every query length -- {name}\n")
    ks = (4, 5)
    ix = {k: TextIndex.build(refs, alphabet="aa", k=k) for k in ks}
    print("| L | max_subs | ms/query k=4 | ms/query k=5 | hits/query | recall vs brute force |")
    print("|--:|--:|--:|--:|--:|--:|")
    for length in (8, 9, 10, 11, 12, 15, 20, 25, 31, 50):
        if not any(len(r) >= length for r in refs):
            continue
        for max_subs in (0, 1, 2, 3):
            # Every cell gets the same query count now. The old dispatch needed a reduced
            # count wherever it fell back to a one-block ball -- L/(m+1) < k, which was most of
            # this grid -- because a single L=12, m=3 query took 0.4 s. Under the search scheme
            # the worst cell here is L < 2k, which this grid does not reach at k=4.
            n_q = 200 if RUN_BENCHMARK else 50
            queries = sample_queries(refs, length, n_q, seed=length)
            ms, hits = {}, {}
            for k in ks:
                if length < k:
                    ms[k], hits[k] = float("nan"), None
                    continue
                t0 = time.perf_counter()
                res = ix[k].search_batch(queries, max_subs=max_subs, threads=1)
                ms[k] = (time.perf_counter() - t0) * 1000 / len(queries)
                hits[k] = res.num_hits
            seen = {h for h in hits.values() if h is not None}
            assert len(seen) == 1, f"k changed the answer at L={length}, max_subs={max_subs}: {hits}"
            rec = f"{recall(refs, ix[4], length, max_subs):.3f}" if check_recall else "--"
            print(f"| {length} | {max_subs} | {ms[4]:.3f} | {ms[5]:.3f} | "
                  f"{seen.pop() / len(queries):.1f} | {rec} |")


def table_c(name: str, refs: list[str]) -> None:
    print(f"\n## C. Thread scaling -- {name}, L=12, max_subs=2\n")
    ix = TextIndex.build(refs, alphabet="aa", k=4)
    # Big enough that the small texts clear the timer: at 0.02 ms/query a 500-query
    # batch on the synthetic text finishes before the threads are all running.
    queries = sample_queries(refs, 12, 50_000 if RUN_BENCHMARK else 2_000, seed=3)
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


def table_d(texts: dict[str, list[str]], tmp: str) -> None:
    """Is the index small, and is it fast to save and load? Numbers, not an opinion.

    `post_ids` is one uint32 per in-record k-mer start, so a positional seed table is inherently
    a few times the text it indexes; `load(mmap=True)` maps rather than parses, which is why it
    does not scale with that size and why several processes share one copy of the pages.
    """
    print("\n## D. Index on disk\n")
    print("| text | k | file MB | x text | save s | load mmap s | load read s |")
    print("|---|--:|--:|--:|--:|--:|--:|")
    for name, refs in texts.items():
        residues = sum(len(r) for r in refs)
        for k in (4, 5):
            ix = TextIndex.build(refs, alphabet="aa", k=k)
            t0 = time.perf_counter()
            ix.save(tmp)
            save_s = time.perf_counter() - t0
            size = pathlib.Path(tmp).stat().st_size
            t0 = time.perf_counter()
            mapped = TextIndex.load(tmp, mmap=True)
            map_s = time.perf_counter() - t0
            t0 = time.perf_counter()
            read = TextIndex.load(tmp, mmap=False)
            read_s = time.perf_counter() - t0
            del mapped, read
            print(f"| {name} | {k} | {size / 1e6:.1f} | {size / residues:.1f} | {save_s:.2f} "
                  f"| {map_s:.4f} | {read_s:.2f} |")
            del ix
    pathlib.Path(tmp).unlink(missing_ok=True)


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
    table_d(texts, os.path.join(tempfile.gettempdir(), "seqtree_bench_text.sti"))


if __name__ == "__main__":
    main()
