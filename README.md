# seqtree

[![CI](https://github.com/antigenomics/seqtree/actions/workflows/ci.yml/badge.svg)](https://github.com/antigenomics/seqtree/actions/workflows/ci.yml)
[![Docs](https://github.com/antigenomics/seqtree/actions/workflows/docs.yml/badge.svg)](https://antigenomics.github.io/seqtree/)

Fast fuzzy search over biological sequences (amino-acid or nucleotide), as a C++
core with a minimal Python binding. Build an immutable index once, then search
single queries or massive batches in parallel.

Two search engines over one trie:

- **`seqtm`** — branch-and-bound enumeration. Exact per-type edit caps
  (`max_subs` / `max_ins` / `max_dels`) and a fast Hamming-only path. Best for
  small edit distances (UMI collapse, error correction, CDR3/epitope matching).
- **`seqtrie`** — banded edit-distance DP. Matrix-weighted score budgets
  (BLOSUM62 + gap costs), cost independent of edit count. Best for
  similarity-scored searches.

`engine="auto"` picks one per query. Results are payload-agnostic:
`(ref_id, score, n_subs, n_ins, n_dels)`. Downstream libraries map `ref_id` back
to their own payloads (V gene, MHC, counts) and filter.

## Build

```fish
bash setup.sh            # repo-local .venv + editable install
bash setup.sh --tests    # + pytest
bash setup.sh --bench     # + benchmark deps (huggingface_hub, pandas, psutil)
```

## Quickstart

```python
import seqtree

idx = seqtree.Index.build(["CASSLAPGATNEKLFF", "CASSLELGATNEKLFF"], alphabet="aa")

p = seqtree.SearchParams(max_subs=2, engine="seqtm")
for hit in idx.search("CASSLAPGATNEKLFF", p):
    print(hit.ref_id, hit.score, hit.n_subs)

# parallel batch (releases the GIL)
results = idx.search_batch(queries, p, threads=0)   # 0 = all cores

# matrix-weighted budget
pm = seqtree.SearchParams(matrix="BLOSUM62", max_penalty=12, engine="seqtrie")
top = idx.search_top("CASSLAPGATNEKLFF", pm, k=5)

# alignment on demand
aln = idx.align(0, "CASSLELGATNEKLFF", p)
print(aln.aligned_query, aln.aligned_ref, aln.ops)

# batch-vs-batch (auto-indexes the larger set)
pairs = seqtree.pairwise_batch(query_set, db_set, p, alphabet="aa")
```

## Tests

```fish
cmake -S . -B build -G Ninja -DSEQTREE_TESTS=ON
cmake --build build
ctest --test-dir build           # C++ unit tests
pytest tests/python              # Python tests
```

## Benchmarks

```fish
python bench/bench.py                                   # fast tier (real VDJdb data)
env RUN_BENCHMARK=1 python bench/bench.py --sizes 1000000 --queries 1000000 --threads 16
```

## Development

This repo follows **git-flow**:

- `master` — stable, release-ready; CI + docs deploy run here.
- `dev` — integration branch for day-to-day work.
- feature branches branch off `dev` and merge back via PR; releases merge `dev` → `master`.

Roadmap (affine gaps, position-specific matrices, e-value / significance via
control-set and tf-idf, succinct memory packing) lives in [docs/roadmap.rst](docs/roadmap.rst).
