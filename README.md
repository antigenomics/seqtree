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
- **`seqtrie`** — full-width edit-distance DP carried down the trie. Honours the
  `max_penalty` score budget only; it **ignores the per-type edit caps**. Use it
  when the budget is the whole specification.

`engine="auto"` always picks `seqtm`, because it is the only engine that enforces
the caps you asked for. Results are payload-agnostic:
`(ref_id, score, n_subs, n_ins, n_dels)`. Downstream libraries map `ref_id` back
to their own payloads (V gene, MHC, counts) and filter.

Beyond search, seqtree ships:

- **Substitution matrices** — built-in `identity`, `BLOSUM45`, `BLOSUM62`, `BLOSUM80`, `PAM250`, `PAM100`, and `structural`
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
- **Calibrated cutoffs** — `threshold_for_evalue` inverts the E-value into the score cutoff that
  achieves it, **per query**. A fixed cutoff is not a calibrated one: a control repertoire is
  dense near germline and sparse among rare junctions, so the same threshold buys a common query
  far more chance neighbours than a rare one.
- **Gap-block alignment** — `gapblock` restricts alignment to one contiguous indel, which is the
  right model for a V(D)J junction and, measured against unrestricted affine alignment, is
  exactly optimal on **98.8%** of genuinely related pairs at a calibrated `gap_open`. A gap
  prior (`central_prior`, `profile_prior`, `frame_prior`) chooses where the block goes — a
  sequence score alone cannot. `score_matrix` scores a whole query set against a whole reference
  set in one GIL-released C++ call (**532 M pairs/s** on 16 cores; `numpy.asarray` wraps the
  result with no copy), the shape a prototype-distance embedding needs.
- **Pairwise alignment without BioPython** — `seqtree.pairwise` is Needleman–Wunsch
  (`mode="global"`) and Smith–Waterman (`mode="local"`) with affine or linear gaps, on the raw
  log-odds scale. It is a **drop-in for `Bio.Align.PairwiseAligner`** — verified against it as an
  oracle across three matrices, ten gap/mode settings and sixty sequence shapes with **zero
  disagreements** — and **65–87× faster**, since there is no Python in the per-pair loop.
  `dist_matrix` gives `d = s(a,a) + s(b,b) − 2·s(a,b)` directly. BioPython is a *test-only*
  dependency; seqtree still needs nothing at runtime.
- **Island profiles** — `IslandProfile.fit` builds a position weight matrix over a set of
  frame-aligned junctions (an *island*) and scores a query column by column against the island
  consensus, as a non-negative penalty that flows through `threshold_for_evalue` unchanged. At a
  repertoire-scale cutoff it recovers **48.5%** of held-out members against **37.6%** for
  min-over-members; at a loose cutoff the two are indistinguishable, so it earns its keep only
  where the cutoff is strict.

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
bash setup.sh --bench     # + benchmark deps (huggingface_hub, psutil)
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

# ...and the cutoff that achieves a target E, per query (-1 = unreachable at this control size)
ceiling = seqtree.SearchParams(max_subs=14, max_penalty=50, matrix="BLOSUM62", engine="seqtm")
thetas = seqtree.threshold_for_evalue(target, control, queries, ceiling, e_target=0.05)

# one contiguous gap block, placed by a prior rather than by the score alone
from seqtree.gapblock import GapBlockIndex, central_prior, embed_in_frame

gbi = GapBlockIndex(cdr3s, "aa", d_max=2)
mat = seqtree.SubstitutionMatrix.blosum62()
for ref_id, score, block_len, block_pos in gbi.search(
        "CASSLGQAYEQYF", 40, mat, gap_open=2 * mat.scale(),
        gap_prior=central_prior(int(1.5 * mat.scale()))):
    ...

# a fixed frame column makes gap placement transitive -- and a column index, hence a PWM, possible
embed_in_frame("CASSGQAYEQYF", width=14, c=4)      # 'CASS--GQAYEQYF'

# a whole query set vs a whole reference set, in one GIL-released C++ call
from seqtree.gapblock import score_matrix, IslandProfile
sm = score_matrix(clonotypes, prototypes, mat, gap_open=2 * mat.scale(), threads=0)
import numpy as np
distances = np.asarray(sm)                          # (len(clonotypes), len(prototypes)) int32, zero-copy

# a position weight matrix over an island, still a non-negative penalty (feeds threshold_for_evalue)
profile = IslandProfile.fit(island_members)
profile.score("CASSLGQAYEQYF")                      # 0 on the consensus, > 0 for deviations

# ordinary pairwise alignment -- Needleman-Wunsch / Smith-Waterman, no BioPython
from seqtree.pairwise import align, score, dist_matrix
score("CASSLGQAYEQYF", "CASSPGQAYEQF", mat)                    # global, BLAST defaults (11/1)
score("WWWAAAWWW", "KKKAAAKKK", mat, mode="local")             # Smith-Waterman
score("AAA", "AAAAA", mat, gap_open=5, gap_extend=5)           # linear gaps: open == extend
aln = align("CASSLGQAYEQYF", "CASSPGQAYEQF", mat)              # + aligned strings and ops

d = np.asarray(dist_matrix(v_genes, v_genes, mat, threads=0))  # s(a,a)+s(b,b)-2s(a,b), zero diagonal
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
python bench/bench_gapblock.py       # the gap-freedom ladder: fixed centre → prior → flat → affine
python bench/bench_score_matrix.py   # dense batch gap-block throughput (µs/pair, M pairs/s, RSS)
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
