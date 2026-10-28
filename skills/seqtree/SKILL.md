---
name: seqtree
description: Use seqtree, a C++20 sequence-search and alignment library with Python bindings, for fuzzy indexed search, global/local alignment, sequence distances, substring search, gap-block scoring, and batched reductions.
---

# Using seqtree

seqtree is a payload-agnostic C++20 library with nanobind Python bindings and
zero Python runtime dependencies. Keep annotations in caller-owned tables and
map returned `ref_id` values back to input reference order. It supplies sequence
operations, not domain datasets, clustering policy, or biological interpretation.

## Choose the operation

| Task | Python API | C++ entry point/header |
|---|---|---|
| Bounded fuzzy search over whole sequences | `Index`, `SearchParams` | `Index`, `Searcher`, `seqtree/seqtree.hpp` |
| Global/local affine similarity alignment | `seqtree.pairwise.score`, `align`, `score_matrix` | `align_score`, `align_pair`, `align_score_matrix`, `seqtree/seqtree.hpp` |
| Hamming/Levenshtein distance | `seqtree.distance` | `hamming`, `levenshtein`, matrix variants, `seqtree/seqtree.hpp` |
| Substring occurrences in whole records | `TextIndex` | `seqtree/text_index.hpp` |
| One contiguous gap block | `seqtree.gapblock` | gap-block functions in `seqtree/seqtree.hpp` |
| K-mer seed lookup | `KmerIndex` | `seqtree/kmer_index.hpp` |
| Control-counted significance | `seqtree.evalue`, text E-value helpers | compose native searches/counts; statistical helpers are Python |

Read the installed signatures or the public headers for the version being used.
Python-only convenience layers do not imply matching C++ classes.
In a source checkout, consult [API contracts](../../docs/api.rst),
[gap-block operations](../../docs/gapblock.rst), and
[E-value methodology](../../docs/evalue.rst) for advanced reductions and layouts.

## C++ integration

Use the source tree's CMake target; it does not require building Python bindings:

```cmake
add_subdirectory(path/to/seqtree)
target_link_libraries(my_program PRIVATE seqtree_core)
```

The library requires C++20 and links the platform thread library through its target.
`SEQTREE_PYTHON=ON` enables nanobind bindings; tests and benchmarks are separate
`SEQTREE_TESTS` and `SEQTREE_BENCH` options.

```cpp
#include <seqtree/seqtree.hpp>

int main() {
    auto index = seqtree::Index::build({"CAT", "CAR", "WST"},
                                     seqtree::Alphabet::AminoAcid);
    seqtree::SearchParams params;
    params.engine = seqtree::Engine::SeqTm;
    params.max_substitutions = 1;
    seqtree::Searcher searcher(*index);
    auto hits = searcher.search("CAT", params);
    return hits.size() == 2 ? 0 : 1;
}
```

C++ field names differ from Python names: for example,
`max_substitutions` corresponds to Python `max_subs`.
`Index::align(query, ref_id, params)` uses query-first argument order;
Python `index.align(ref_id, query, params)` uses reference-first order.

## Python usage

Install the package with `pip install seqtree` in the project's environment.
For source development, use the repository build instructions; editable installs
have automatic native rebuilding disabled, so reinstall after C++ edits.

```python
import seqtree
from seqtree.pairwise import score, align

index = seqtree.Index.build(["CAT", "CAR", "WST"], alphabet="aa")
params = seqtree.SearchParams(max_subs=1, engine="seqtm")
rows = index.search_batch(["CAT", "CAR"], params, threads=1)
assert {hit.ref_id for hit in rows[0]} == {0, 1}

matrix = seqtree.SubstitutionMatrix.blosum62()
similarity = score("CAT", "CAR", matrix, mode="global", gap_open=11, gap_extend=1)
alignment = align("CAT", "CAR", matrix, mode="local", gap_open=11, gap_extend=1)
```

## Score and search contracts

- Search minimizes nonnegative penalties. Pairwise global/local alignment maximizes
  signed substitution similarity. `pairwise_batch` is bounded trie search, not an
  exhaustive aligner; use `pairwise.score_matrix` for every pair's similarity.
- `auto` selects `seqtm`, which enforces edit caps and penalty limits.
  `seqtrie` ignores edit caps and uses the penalty budget; specify its ceiling.
  Default edit caps are zero: supplying a matrix alone does not enable substitutions.
- Matrices retain both similarity and Gram-transformed penalty views.
  Gap costs must use the chosen view's units: a useful penalty starting point is
  `2 * matrix.scale()`; raw BLOSUM62 similarity commonly uses open 11, extend 1.
  A gap of length L costs `open + (L-1)*extend`; equal costs give linear gaps.
- Gap-block scoring permits one contiguous block, not arbitrary affine gaps.
  Positional weights scale mismatch penalties; zero weights do not supply match evidence.
- Hamming requires equal lengths. Plain distance functions compare case-sensitive
  bytes; encoded search uses its alphabet's normalization and rejects invalid queries.
  Amino-acid `B/Z/X/*` are literal symbols, not wildcard instructions.

## Batch operations and memory

Build an index once and submit query cohorts through native batch APIs. Do not
wrap an internally threaded search in a per-query worker pool.
Native batch calls release the Python GIL. Positive `threads=B` sets the total
CPU budget including the caller; nonpositive values select the detected allocation.
Divide that allocation between independent simultaneous calls.

Dense score matrices require space proportional to queries times references.
Use existing top-k, count, grouped, consumer, or target-max reductions when the
full matrix is not the requested output. Consult their documented tie handling,
byte limits, identity predicates and output ordering before composing reductions.

`TextIndex` indexes whole records once for multiple query lengths.
Ungapped search is bounded by substitutions; gapped search additionally requires
sufficient query length for its exact seeds. Check `truncated` before interpreting
capped results as complete. `matrix=` scores accepted ungapped text hits rather
than changing the Hamming acceptance predicate.

## Optional statistical calibration

`evalue_result(n_target, n_control, n_ref, m_control)` reports
`E = (n_ref / m_control) * n_control`, with Poisson any-hit and enrichment tails.
An empty control ball returns the approximate upper bound `3*n_ref/m_control`
and flags `rule_of_three`; an empty control population yields infinite E.
A meaningful control must match the intended null. These helpers do not fit
negative-binomial models or perform multiple-testing correction.

Keep exact hits for ordinary annotation. Exclusion filters hits without rebuilding
indices, refitting models, or changing population sizes. The legacy `evalues`
helper excludes zero-score hits, which can include nonidentical sequences under
masked penalties; other reducers use encoded identity. Check the specific API.
