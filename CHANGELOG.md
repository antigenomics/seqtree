# Changelog

All notable changes to `seqtree`. Dates are release dates; the project is pre-1.0, so a **minor**
bump may carry breaking changes.

## [0.3.0] — unreleased

Gap-block alignment, calibrated cutoffs, seed significance — and the removal of several engine
paths that returned confident wrong answers.

### Breaking

- **`engine="auto"` now always resolves to `seqtm`.** It previously routed matrix-plus-indel
  searches to `seqtrie`, whose budget defaults to `INT_MAX/4`, so
  `SearchParams(max_subs=1, max_ins=1, matrix="BLOSUM62")` silently returned **every reference in
  the index**. `seqtrie` ignores per-type edit caps and can never honour them; `auto` will not pick
  it. Passing a matrix to `seqtrie` without an explicit `max_penalty` now raises.
- **`Mode::Local` / `mode="local"` deleted.** It was a no-op: the field was stored and never read,
  with zero call sites across `seqtree`, `vdjmatch` and `mhcmatch`.
- **`GapPrior` takes `(block_start, block_length, longer_length)`** instead of
  `(block_start, shorter_length)`. `central_prior`'s output is bit-identical
  (`|2i − (m−d)| == |2i + d − m|`), and a frozen table pins that.
- **`gap_extend` is now honoured.** `Index.align` is a real Gotoh affine alignment; previously
  `gap_extend=1` and `gap_extend=99` produced identical scores and ops.

### Fixed

- `pairwise_batch` silently inverted `n_ins`/`n_dels` when transposing, but only when
  `len(a) >= len(b)` — a size-dependent inversion.
- `Index.align` did not validate the query alphabet (unlike `search_into`), and accepted negative
  gap costs.
- `Limits::max_hits` was set and never read; in `mode="all"` it truncated an **unsorted** list.
- Stale docstring on `SubstitutionMatrix.from_similarity` (claimed `max(s_aa, s_bb) − s_ab`; the
  code is the Gram form `s_aa + s_bb − 2·s_ab`).

### Added

- **`SubstitutionMatrix.scale()`** — the median mismatch penalty (BLOSUM62 → 14). Gap costs must
  live on the matrix's scale; the old `gap_open=1` default made gaps ~14× cheaper than
  substitutions, so `align()` would gap an equal-length pair rather than substitute. Use
  `gap_open = 2 * matrix.scale()`.
- **`seqtree.gapblock`** — single-contiguous-gap-block alignment for anchored junction loops.
  `gapblock_score` is the exact `O(min(m,n))` optimum (0 mismatches in 55,727 pairs against
  brute-force layout enumeration). `GapBlockIndex` reuses the existing Hamming engine over
  deletion variants; no new C++.
- **Gap priors** — `central_prior(lam)`, `profile_prior(lam, w)`, `frame_prior(lam, c)` and
  `embed_in_frame(seq, width, c)`. A sequence score alone cannot place the block: a hard central
  pin agrees with the flat, score-only choice on only 10.6% of pairs.
- **`threshold_for_evalue` / `thetas_from_scores`** — invert `Ê = (N/M)·n_C` into the score cutoff
  that achieves a target E, **per query**. Exact rather than a root-find, because scores are
  integers. Returns `-1` where `e_target < 3N/M`, i.e. where the control is too small to certify
  the bar.
- **`seqtree.seeds`** — `core_kmers` and `SeedIndex` give control-calibrated E-values for shared
  core k-mers. A shared rare central k-mer is ~4× enriched among co-specific pairs, but covers only
  ~0.5% of them: seeds buy precision, not recall.
- **`bench/bench_gapblock.py`** — the gap-freedom ladder, from a hard central pin through priors to
  unrestricted affine.

### Measured

Numbers that constrain the API, all reproducible from `bench/` and the downstream repos:

- **One gap block is enough.** Against a model-independent structural oracle (iterative
  superposition + unrestricted affine DP) over 3,049 crystal junction pairs from 199 unique
  sequences, the true correspondence is a single contiguous block in **95.2–100%** of cases for
  every `d = 1..4`. Forcing one block costs no median CA-RMSD.
- **The restriction is free where it applies, and protective where it does not.** At
  `gap_open = 2*scale`, gap-block equals unrestricted affine on **98.8%** of related pairs (one
  indel + 0–2 substitutions). On *unrelated* pairs affine undercuts it by a median of **106**
  penalty units — affine inventing an alignment that does not exist. Extra gap freedom buys
  manufactured similarity.
- **A fixed score cutoff is not a calibrated cutoff.** Building islands on human TRB by union-find
  at `gapblock_score ≤ 60`, size-matched *random control* junctions land a larger fraction of nodes
  in components of ≥5 (**0.748**) than real same-epitope sequences do (**0.660**). Under per-query
  E-value edges at `E* = 0.05` the picture inverts: real 1.583 edges/node against control 0.022,
  and the control forms no component of size 3 at all.
- **Constraining the block is what buys precision.** Compared at a *matched* false-positive rate —
  each rung given the cutoff at which its own ball admits `E*` chance neighbours, since a freer rung
  finds lower scores and a fixed budget would reward it for that — retrieval precision on the
  length-different fraction of VDJdb same-epitope pairs is **0.65** for a hard central pin or a
  central prior, and **0.31** when the score chooses freely among `L+1` placements *or* among five
  plausible ones. Trying several positions and keeping the best score is worse than not trying.
- **Performance.** 91% of `GapBlockIndex.search` time is the query-deletion-variant branch, 9% the
  9.8M-entry auxiliary indices. Netting the prior out of each variant's budget cuts that branch
  from ~15 sub-searches to 2.5. Variant dedup (7–10% of variants before pruning, fewer after) and
  length-bucketing are **not** built. At budget 40 over 250k references, `d_max=2` gap-block search
  costs 2,562 µs/query — less than the plain Hamming ball at the same budget (3,051 µs/query).

### Docs

- `skills/seqtree/SKILL.md` — public API surface, invariants, and the gotchas that have bitten.
- `docs/gapblock.rst` — a worked guide: why one gap block, how to choose its position, why a
  placement rule is a column frame, and why a fixed score cutoff is not a calibrated one.
- README corrected: `seqtrie` is a full-width DP that ignores per-type caps, not a banded one, and
  `auto` does not choose between engines.
- `appendix/evalue.tex` gains §"The score model: one gap block, placed by a prior" (the appendix
  derived a theory of balls without ever saying what the score was) and a remark inverting the
  E-value into a per-query cutoff. The pMHC section is compacted from ~110 lines of prose to ~50,
  deferring to the `mhcmatch` appendix, which specialises this one rather than repeating it. Its
  empirical tables stay: they are this repo's own `bench/bench_mhc_guess.py` output.
- Test coverage: `gapblock.py`, `evalue.py` and `seeds.py` at **100%**; package total 88%. A new
  `tests/python/test_doc_coverage.py` fails the build if a public symbol is undocumented, missing
  from `__all__`, or unreachable from any docs page.

## [0.2.0]

- `structural` substitution matrix: Miyazawa–Jernigan interaction-strength similarity.
- Built-in matrix list: `identity`, `BLOSUM62`, `PAM250`, `PAM100`, `structural` (dropped PAM50).

## [0.1.0]

- `SubstitutionMatrix.penalty(a, b)` exposed to Python.

## [0.0.3]

- Reproducible table→plot benchmark pipeline with oracle + perf regression.
- pMHC non-binder E-value filter; class-II promiscuity notes.
