# seqtree

[![PyPI](https://img.shields.io/pypi/v/seqtree.svg)](https://pypi.org/project/seqtree/)
[![Python](https://img.shields.io/pypi/pyversions/seqtree.svg)](https://pypi.org/project/seqtree/)
[![License](https://img.shields.io/badge/license-GPLv3-green)](LICENSE)
[![CI](https://github.com/antigenomics/seqtree/actions/workflows/ci.yml/badge.svg)](https://github.com/antigenomics/seqtree/actions/workflows/ci.yml)
[![Docs](https://github.com/antigenomics/seqtree/actions/workflows/docs.yml/badge.svg?branch=dev)](https://antigenomics.github.io/seqtree/)

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

Beyond search, seqtree ships:

- **Substitution matrices** — built-in `identity`, `BLOSUM62`, `PAM250`, `PAM100`, and `structural`
  — a **Miyazawa–Jernigan interaction-strength** matrix: each residue's strength `q(a)=mean_b e(a,b)`
  is read off the MJ contact potential, so substitutions between residues of like interaction strength
  are cheap. It separates strong (hydrophobic `F W C L Y M I V`) from weak (polar/charged
  `S Q D E K`) interactors — the strong/weak-interactor axis of TCR-recognition models
  ([Košmrlj et al., *PNAS* 2008](https://doi.org/10.1073/pnas.0808081105); MJ contact energies from
  Miyazawa & Jernigan, *J Mol Biol* 1996) — letting dissimilar-but-chemically-equivalent loops align.
  Plus custom matrices via `SubstitutionMatrix.from_similarity` (Gram penalty `s(a,a)+s(b,b)−2·s(a,b)`).
- **E-values / significance** — calibrate hit counts against a background control repertoire
  (`load_control` + `evalues`), the TCRNET approach on a finite-sample footing. See the
  [E-value guide](https://antigenomics.github.io/seqtree/evalue.html).

## Install

```fish
pip install seqtree       # prebuilt wheels for CPython 3.10–3.13
```

Prebuilt wheels cover **Linux x86-64**, **macOS arm64 (Apple Silicon)**, and **Windows x86-64**.
There are **no Intel/x86-64 macOS wheels** — Intel Macs build from source (see below), which just
needs a C++17 compiler and CMake (pulled in automatically by the build).

## Build from source

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

# E-values against a background control repertoire (TCRNET-style significance)
control = seqtree.load_control("human_trb_aa", size=1_000_000)
target = seqtree.Index.build(vdjdb_cdr3s, alphabet="aa")
for q, r in zip(queries, seqtree.evalues(target, control, queries, p)):
    if r["p_enrichment"] < 1e-3:
        print(q, r["E"], r["n_target"], r["n_control"])
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
python bench/bench_gnuplot.py        # throughput / scaling / matrix / collisions → SVG (needs gnuplot)
python bench/bench.py                # recall vs ground truth (real VDJdb data)
python bench/bench_evalue.py         # true E-value benchmark (target vs background control)
python bench/bench_evalue_matrix.py  # significance across reference/control/query/scope grid
python bench/bench_epitope.py        # epitope detection-complexity (GIL vs NLV)
```

Figures (throughput, scaling, matrix-scoring overhead, collisions, E-value matrix, epitope
detection) and the full methodology are in the [benchmarks docs](https://antigenomics.github.io/seqtree/benchmarks.html).
Set `RUN_BENCHMARK=1` for the large tiers.

## Development

This repo follows **git-flow**:

- `master` — stable, release-ready; CI + docs deploy run here.
- `dev` — integration branch for day-to-day work.
- feature branches branch off `dev` and merge back via PR; releases merge `dev` → `master`.

Roadmap (affine gaps, position-specific matrices, succinct memory packing) lives in
[docs/roadmap.rst](docs/roadmap.rst). Control-set E-values already ship — see the
[E-value guide](https://antigenomics.github.io/seqtree/evalue.html).
