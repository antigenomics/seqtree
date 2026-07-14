# Changelog

All notable changes to `seqtree`. Dates are release dates; the project is pre-1.0, so a **minor**
bump may carry breaking changes.

## [0.4.0] — 2026-07-11

### Added

- **`seqtree.pairwise` — Needleman-Wunsch and Smith-Waterman, so ordinary protein alignment no
  longer needs BioPython.** Everything else in seqtree *minimises a non-negative penalty*, which
  is what a search ball and an E-value need. These **maximise a raw log-odds similarity**, the way
  BLAST and BioPython do, because that is what a pairwise alignment means.

  | | |
  |---|---|
  | `pairwise.score(q, r, matrix, mode=...)` | optimal score, `O(min(m,n))` memory |
  | `pairwise.align(...)` | plus the aligned strings and ops |
  | `pairwise.score_matrix(queries, refs, ...)` | dense `n × K`, GIL released, zero-copy numpy |
  | `pairwise.dist_matrix(...)` | `d = s(a,a) + s(b,b) − 2·s(a,b)`: non-negative, zero on the diagonal |

  `mode="global"` is Needleman-Wunsch, `mode="local"` Smith-Waterman, and **`gap_open == gap_extend`
  gives linear gaps** — no separate mode. A gap run of length `L` costs `gap_open + (L-1)·gap_extend`,
  and global charges end gaps (true NW, not semi-global).

  **It is a drop-in.** `tests/python/test_pairwise.py` runs it against `Bio.Align.PairwiseAligner`
  as an oracle over three matrices × ten gap/mode settings × sixty sequence shapes — **zero
  disagreements**, including on real germline V genes. BioPython is a *test-only* dependency;
  seqtree still has **zero required runtime dependencies** and never imports it.

  Measured on an M3 (all-against-all, BLOSUM62, global, 11/1):

  | sequence length | seqtree, 1 thread | seqtree, 16 threads | BioPython | speedup |
  |---|---|---|---|---|
  | 15 (a junction) | 1.7 M pairs/s | **20.1 M pairs/s** | 0.31 M pairs/s | **65×** |
  | 90 (a germline V gene) | 72 k pairs/s | **893 k pairs/s** | 10 k pairs/s | **87×** |

- **`SubstitutionMatrix.similarity(a, b)`** — the raw signed log-odds, alongside the existing
  non-negative `penalty(a, b)`. The Gram transform `pen = s(a,a) + s(b,b) − 2·s(a,b)` is **lossy**:
  it forces the diagonal to zero and destroys `s(a,a)`, so a similarity cannot be recovered from a
  penalty. Both views are now stored.

- **`SubstitutionMatrix.blosum45()` and `.blosum80()`**, and the names `"BLOSUM45"` / `"BLOSUM80"`
  wherever a matrix name is accepted. Shallower and deeper than BLOSUM62 — for remote and close
  homologs respectively.

## [0.3.1] — 2026-07-10 (never published; folded into 0.4.0)

> Tagged but not released to PyPI. Everything below ships in **0.4.0**, so upgrading from 0.3.0
> straight to 0.4.0 picks it all up. Kept as its own section because it is a distinct set of fixes.


### Fixed

- **A cold cache shared by concurrent processes could hand back a half-written index.**
  `Index::save` wrote straight into the destination, so for the whole duration of the write the
  file existed but was truncated. A second process that checked `os.path.exists(cache)` in that
  window loaded a stub and raised `RuntimeError: truncated or corrupt index`. On a 45 MB control
  index the window is ~55 ms, and a reader racing a writer hit it **10 times out of 10**.

  This is the first-use-only failure of any multi-process fan-out sharing `~/.cache`: pytest-xdist,
  a Snakemake or Nextflow pipeline calling `load_control` in parallel, a `multiprocessing` pool.
  Once the cache is warm it is read-only and was always safe. CI matrix jobs were never affected —
  separate runners, separate caches.

  `Index::save` and `KmerIndex::save` now serialize into a uniquely-named temporary beside the
  destination and `rename` it into place. Rename is atomic on the same filesystem, on POSIX and
  Windows alike, so a reader sees either the previous complete file or the new complete file and
  never a partial one. A failed save cleans up its temporary and leaves any pre-existing index
  intact.

- **A corrupt or stale cache now rebuilds instead of raising.** A file truncated by a full disk,
  left by a killed process, or written by an older seqtree sent `load_control` into an exception;
  it now falls back to rebuilding. The cache was always best-effort and now behaves that way.

- **The control cache is content-addressed, so a stale cache can no longer be served silently.**
  The key was `control_{name}_{size}.sqtree`, which named neither the **alphabet**, nor the
  **seed**, nor the **source data**. Three consequences, all live:

  - Two calls differing only in `seed` — which must draw *different* reservoir samples — shared one
    cache file, so the second silently received the first's sequences. Same for `alphabet`.
  - An upgrade that changed the bundled control kept the same filename, so a warm cache served the
    **previous release's** control. This is exactly how 0.3.0's corrected (uniform) control could be
    masked by a stale 0.2.0 (abundance-head) cache — and why 0.3.0's notes had to ask people to
    `rm ~/.cache/seqtree/control_*.sqtree` by hand.

  The key now carries a fingerprint of the bundled asset's own bytes (or, on the download path, the
  source and seed), so a control that changed simply misses the old cache. Superseded caches,
  including pre-fingerprint ones from earlier releases, are deleted on the next build.

  **You no longer need to clear `~/.cache/seqtree` when upgrading.** Doing so is harmless.

### Added

- **`load_control` takes an inter-process lock around build-and-save when `filelock` is available**
  (it arrives with `huggingface_hub`). This is an optimisation, not the fix: correctness comes from
  the atomic rename and holds with no lock at all. What the lock saves is work — without it, a cold
  fan-out of N workers has every worker build the same 250k-clonotype index and discard N−1 of them.
  seqtree still has **zero required runtime dependencies**; the import is guarded.

## [0.3.0] — 2026-07-10

Gap-block alignment, calibrated cutoffs, seed significance — the removal of several engine paths
that returned confident wrong answers, and a corrected background control that changes every
E-value.

### Breaking

- **The bundled control is a different set of sequences.** It was the *abundance head* of the
  upstream repertoire — the 250,000 most expanded clonotypes — because both `gen_control.py` and
  `_download` took the first `size` unique rows of a count-descending table. `appendix/evalue.tex`
  (`ass:indep`) assumes the control's unique clonotypes are i.i.d. from `P₀`. Measured against a
  uniform sample of the same size, the head is **25.8× more self-similar** (P(Hamming≤2 | equal
  length) 3.11×10⁻³ vs 1.20×10⁻⁴) and carries **3.1× the ball mass** at a BLOSUM62 budget of 40
  (mean n_C 110.1 vs 35.5). Both are now uniform reservoir samples over unique **productive**
  clonotypes, seeded and shuffled so any prefix is itself a valid sub-sample.

  **Every E-value moves.** Delete `~/.cache/seqtree/control_*.sqtree` after upgrading — a warm
  cache from 0.2.0 would otherwise be served in place of the corrected control. (Fixed in 0.3.1:
  the cache is now content-addressed and a stale one simply misses. Upgrading straight from 0.2.0
  to ≥0.3.1 needs no manual step.) Numbers derived from the control are corrected throughout this
  file, `seeds.py`, `SKILL.md` and the appendix.
- **Controls are filtered to productive clonotypes.** VDJtools marks out-of-frame rearrangements
  with `_` and in-frame stops with `*`; 13.7% of the mouse TRB table is out of frame. `_` cannot be
  repaired at the amino-acid level — VDJtools collapses a *run* of untranslatable positions into one
  character, so the residue count is already gone — and out-of-frame junctions escape thymic
  selection, making them an estimator of `P_gen`, which `lem:hierarchy` says is not `P₀`.
  `load_control("mouse_trb_aa")` previously raised on `_`; it now yields 694,241 clonotypes.
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
- **`gapblock.score_matrix` and `ScoreMatrix`** — the dense `n × K` counterpart of
  `GapBlockIndex.search`, for prototype-distance embeddings, where nothing can be pruned because
  the distance to every reference *is* the output. C++, GIL released, one thread per core. On an
  M3 against 3,000 prototypes: **51.3 M pairs/s** single-threaded, **532.7 M** on 16 cores, versus
  0.41 M for pure-Python `gapblock_score` — while evaluating all `L+1` block positions, not a
  fixed shortlist. The prior is flattened once into an `[m][d][i]` cube, so the kernel never
  re-enters Python. `ScoreMatrix` carries the CPython buffer protocol: `numpy.asarray` wraps it
  without copying, and seqtree keeps its zero runtime dependencies.
- **`gapblock.positions_prior(starts)`** — restrict the block to a fixed set of starts, negative
  values counting from the end, reproducing the `gap_positions=(3, 4, -4, -3)` convention that
  other junction aligners hardcode. Shipped for interoperability, not as a recommendation: at a
  matched false-positive rate on human TRB, candidate starts reach precision 0.156 against 0.414
  for a single hard-pinned centre.
- **`gapblock.IslandProfile`** — a per-island position weight matrix whose column penalty is
  measured against the column's own consensus, `pen(j, a) = round(lam·log(p_max_j / p_j(a)))`. A
  textbook log-odds score is signed and therefore not a ball; this one is `>= 0`, zero on the
  consensus, and flows through `thetas_from_scores` unchanged. The frame column defaults to the
  entropy-optimal one, which is modal at `c = 6` on real islands — where crystal structures put
  the block.

  Whether it beats scoring against every member **depends entirely on the cutoff**, which moves
  with `N`: the E-value's `k = floor(e_target·M/N)` is how many control neighbours the cutoff may
  admit, so the FPR is `k/M`. Over 108 calibrated VDJdb islands of ≥10 members (human TRB, three
  held-out splits, paired bootstrap over islands, 250k control negatives):

  | regime | FPR | min-over-members | `IslandProfile` | difference [95% CI] |
  |---|---|---|---|---|
  | loose reference | 1% | **99.5%** | 99.1% | −0.40 [−1.09, +0.14] |
  | per-epitope islands (`N`= group, median 88) | 0.0568% | 88.3% | **89.3%** | +0.93 [−0.80, +2.79] |
  | repertoire annotation (`N`≈20k) | 0.0012% | 37.6% | **48.5%** | +10.90 [+7.69, +14.21] |

  So: **no significant difference while building the islands**, a large one when using them to
  annotate a repertoire (on islands ≥50 members, 9.8% vs 22.6%). At `N≈20k` and `e_target=0.05`,
  `k=0` and `thetas_from_scores` returns `-1` — the rule of three certifies no `E` below
  `3N/M = 0.236`, and that is the cutoff the third row uses.

  It does **not** generalise: same-epitope junctions in a different island are recovered by
  neither representation. Nor is it a compression — 1,176 B against 182 B of member strings,
  break-even at 84 members.
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
  at `gapblock_score ≤ 60`, **31.7%** of size-matched *random control* junctions land in a component
  of ≥5 — structure invented by the threshold. Per-query E-value edges at `E* = 0.05` cut that to
  **0.000** while raising the real signal: 2.334 edges per node against 0.021 for the control, which
  forms 19,248 singletons, 223 pairs, three components of size 3–4, and nothing larger. The control
  arm's realised edge rate lands on `E*`, which is the check that the calibration is honest.
- **Mouse replicates it.** Against the mouse TRB control (694,241 productive clonotypes), 5 epitopes
  and 1,692 TCRs give 5.856 calibrated edges per node against 0.019 for the control, which again
  forms nothing larger than a pair. At a fixed θ=40 the control still lands 18.3% of its nodes in
  components of ≥5.
- **Constraining the block is what buys precision.** Compared at a *matched* false-positive rate —
  each rung given the cutoff at which its own ball admits `E*` chance neighbours, since a freer rung
  finds lower scores and a fixed budget would reward it for that — retrieval precision on the
  length-different fraction of VDJdb human TRB same-epitope pairs (2,000 queries, `E* = 0.1`):

  | rung | layouts | precision |
  |---|---|---|
  | fixed centre | 1 | **0.414** |
  | central prior λ=21 | ~1–2 effective | 0.336 |
  | flat (score alone) | L+1 | 0.176 |
  | candidates (after 3–4, before last 3–4, centre) | 5 | 0.156 |

  Trying several plausible positions and keeping the best score is worse than not trying: the score
  picks the structurally correct layout about a tenth of the time, so each extra candidate is mostly
  an opportunity to be wrong. Mouse replicates the ordering at `E* = 1.0` but its length-different
  stratum holds only 24–79 true positives, too few to separate the rungs.
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
