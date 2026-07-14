---
name: seqtree
description: Fast fuzzy search over biological sequences — C++20 arena trie + pybind11, with control-calibrated E-values, single-gap-block alignment, and seed significance.
---

# seqtree

Payload-agnostic fuzzy sequence search. Build an immutable index once; search single queries or
batches in parallel. Everything downstream (`vdjmatch`, `mhcmatch`) maps `ref_id` back to its own
payload — seqtree ships **zero runtime dependencies** and no domain data.

Repo: `/Users/mikesh/vcs/code/seqtree`. Venv: `.venv` (`bash setup.sh`). Docs build:
`env -C docs make html` (must be warning-free).

## Engines — pick with care

| engine | enforces `max_subs`/`max_ins`/`max_dels` | honours `max_penalty` | use when |
|---|---|---|---|
| `seqtm` | **yes** | yes | always, unless you specifically want budget-only |
| `seqtrie` | **no** | yes | the score budget is the entire specification |

`engine="auto"` resolves to `seqtm`. It does **not** choose per query; `seqtrie` silently ignores
edit caps, so `auto` can never select it. Passing a matrix to `seqtrie` without an explicit
`max_penalty` raises rather than scanning the whole table.

## Core API

```python
Index.build(refs: list[str], alphabet: str = "aa") -> Index      # 'aa' | 'nt' | 'iupac'
Index.search(query, params) -> list[Hit]                          # Hit(ref_id, score, n_subs, n_ins, n_dels)
Index.search_top(query, params, k=1) -> list[Hit]
Index.search_batch(queries, params, threads=0) -> list[list[Hit]] # releases the GIL
Index.align(ref_id, query, params) -> Alignment                   # true Gotoh affine, on demand
Index.save(path) / Index.load(path)
pairwise_batch(a, b, params, alphabet="aa", threads=0) -> list[list[Hit]]   # a-major

SearchParams(max_subs=0, max_ins=0, max_dels=0, max_total_edits=0, max_penalty=0,
             matrix="", gap_open=1, gap_extend=1, engine="auto", mode="all")
SearchParams.pos_matrix = PositionalMatrix | None
```

**`pairwise_batch` is a bounded trie search, not an aligner.** It returns nothing outside
`max_penalty`. The only affine aligner is `Index.align`.

## Substitution matrices

```python
SubstitutionMatrix.blosum62() / .blosum45() / .blosum80() / .pam250() / .pam100()
                  / .structural() / .unit(n)
SubstitutionMatrix.from_similarity(grid: list[list[int]])   # Gram: pen = s_aa + s_bb - 2*s_ab, clamped >= 0
SubstitutionMatrix.penalty(a, b) -> int                      # pen(a, a) == 0 always
SubstitutionMatrix.similarity(a, b) -> int                   # raw signed log-odds
SubstitutionMatrix.scale() -> int                            # median mismatch; BLOSUM62 == 14
PositionalMatrix.from_weights(base, weights) / .from_tables(size, width, data, masked=[])
```

`scale()` exists because gap costs must live on the matrix's scale. `gap_open=1` against a
BLOSUM62 mismatch of 14 makes gaps ~14× cheaper than substitutions and every alignment
degenerates to gaps. **Use `gap_open = 2 * matrix.scale()`.**

Because `pen(a, a) == 0`, a positional weight scales *mismatch* cost only. It carries **zero
match evidence** — evidence lives in the control-counted E-value, never in a weight.

## E-values and calibrated cutoffs — `seqtree.evalue`

```python
evalues(target, control, queries, params, threads=0, exclude_exact=False) -> list[dict]
    # keys: n_target, n_control, E, p_any, p_enrichment, rule_of_three
    # E = (N/M) * n_control ; rule of three (3N/M) when the control ball is empty

threshold_for_evalue(target, control, queries, params, e_target, threads=0,
                     exclude_exact=False) -> list[int]     # per query; -1 == unreachable
thetas_from_scores(control_scores, n_target, m_control, e_target, theta_max,
                   *, exclude_exact=False) -> list[int]    # the pure-data core
load_control(name="human_trb_aa", size=None) -> Index      # 250k bundled; larger via HuggingFace
    # Cold-cache safe under a multi-process fan-out (xdist / Snakemake / Nextflow): Index.save
    # writes a temp file and renames it into place, so a reader never sees a partial index.
    # A corrupt/stale cache is rebuilt, not raised. filelock (optional) additionally dedupes
    # the build work; correctness does not depend on it.
```

**A fixed score cutoff is not a calibrated cutoff.** The control is dense near germline and
sparse among rare junctions. Measured on human TRB: at `gapblock_score <= 60`, **31.7%** of random
control junctions land in a connected component of ≥5 — structure the threshold invents. Per-query
cutoffs cut that to 0.000 while raising the real signal to 2.334 edges/node. Always invert an
E-value; never pick θ by hand. `-1` means `e_target < 3N/M` — enlarge the control rather than
lowering the bar.

`threshold_for_evalue` needs `params.max_penalty > 0` as the ceiling to search the control at.
One control scan supplies every query's cutoff.

## Pairwise alignment — `seqtree.pairwise` (replaces BioPython)

Everything else here MINIMISES a non-negative penalty. This module MAXIMISES a raw log-odds
similarity, like BLAST/BioPython, because that is what pairwise alignment means.

```python
from seqtree.pairwise import score, align, score_matrix, dist_matrix
score(q, r, matrix, mode="global", gap_open=11, gap_extend=1) -> int   # NW; "local" = SW
align(...) -> Alignment          # + aligned strings/ops; .score is a SIMILARITY, not a penalty
score_matrix(queries, refs, ...) -> ScoreMatrix        # dense n x K, GIL released
dist_matrix(queries, refs, ...)  -> ScoreMatrix        # d = s(a,a)+s(b,b)-2s(a,b), 0 on diagonal
```

- Gap run of length L costs `gap_open + (L-1)*gap_extend`. **`gap_open == gap_extend` = linear
  gaps** (no separate mode). Global charges end gaps (true NW, not semi-global). Local floors at 0.
- **Drop-in for `Bio.Align.PairwiseAligner`**: 0 disagreements over 3 matrices x 10 gap/mode
  settings x 60 shapes, and on real germline V genes. **65-87x faster.** BioPython is a TEST-only
  dep; seqtree still has zero runtime deps.
- `SubstitutionMatrix` now stores BOTH views: `penalty()` (>=0, Gram, for search/E-values) and
  `similarity()` (signed log-odds, for these aligners). The Gram transform is lossy -- it zeroes
  the diagonal -- so a similarity cannot be recovered from a penalty.

## Gap-block alignment — `seqtree.gapblock`

A V(D)J junction's length variation is **one** contiguous indel event, so restrict alignment to
one gap block of length `d = |m - n|`.

```python
gapblock_score(q, r, matrix=None, gap_open=None, gap_extend=1, gap_prior=None) -> (score, block_pos)
gap_cost(d, gap_open, gap_extend) -> int                    # 0 when d == 0 -- MANDATORY guard
deletion_variants(q, d) -> list[(block_pos, variant)]
GapBlockIndex(refs, alphabet="aa", d_max=1).search(query, max_penalty, matrix, gap_open,
                                                   gap_extend=1, gap_prior=None)
    # -> list[(ref_id, score, block_len, block_pos)]

GapPrior = Callable[[block_start_i, block_length_d, longer_length_m], int]
central_prior(lam)          # lam * |block_midpoint - m/2|; lam ~ 1.5 * scale()
profile_prior(lam, w)       # lam * sum(w(j, m) for j in block); w in [0,1], callable
frame_prior(lam, c)         # lam * |i - c|; the block start does not depend on d
positions_prior(starts)     # 0 at each start (negatives from the end), UNREACHABLE elsewhere
embed_in_frame(seq, width, c, gap="-") -> str

# Dense n x K, for prototype-distance embeddings. C++, GIL released, buffer protocol.
score_matrix(queries, refs, matrix=None, gap_open=None, gap_extend=1, gap_prior=None,
             alphabet="aa", threads=0) -> ScoreMatrix     # numpy.asarray(sm) is zero-copy

# Per-island PWM. pen(j,a) = lam*log(p_max_j / p_j(a)): >= 0 and 0 on the consensus, so it is
# still a ball and still feeds thetas_from_scores. c defaults to the entropy-optimal column.
IslandProfile.fit(members, c=None, lam=1000, pseudocount=0.5) -> IslandProfile
    .score(seq) -> int          # UNREACHABLE if it does not fit the frame
    .score_batch(seqs) -> list[int]
    .consensus() -> str         # score(consensus()) == 0
```

**Invariants a prior must satisfy** (nothing else): `>= 0`, and `== 0` when `d == 0`. Otherwise
`s(q, q) != 0` and the score stops defining a ball, breaking the E-value theory and admissible
trie pruning. Monotonicity in `d` is **not** required and is false for `central_prior`.

### What the measurements say

- **One block is enough.** Against a model-independent structural oracle (iterative superposition
  + unrestricted affine DP) over 3,049 crystal junction pairs, the true correspondence is a single
  contiguous block in **95.2–100%** of cases for every `d = 1..4`. Forcing one block costs no
  median CA-RMSD.
- **The restriction is free on related sequences, and protective on unrelated ones.** At
  `gap_open = 2*scale`, gap-block equals unrestricted affine on **98.8%** of related pairs; on
  unrelated pairs affine undercuts it by a median of 106 penalty units, i.e. affine *invents* an
  alignment. Free gapping buys manufactured similarity.
- **A sequence score cannot place the block.** A hard central pin agrees with the flat,
  score-only choice on only 10.6% of pairs. The block sits at the **loop apex** (Cys-offset 6 for
  both TRA and TRB) and does not drift with `d`. `central_prior` hits it 42.4% (TRA) / 30.1%
  (TRB); a germline-untemplated-span rule hits 0.4% / 19.8% and was rejected.
- **A frame is transitive iff the block start is constant in `d`.** That is `frame_prior`.
  `central_prior`'s start drifts, so embedding two shorter members into a common frame relates
  them by *two* blocks — no consistent column index, no PWM. Use `embed_in_frame(seq, W, c)`:
  left-anchor the first `c` residues, right-anchor the rest.
- **Scoring several candidate placements and keeping the best loses.** At a matched FPR on human
  TRB, candidate starts `(3, 4, mid)` reach precision 0.156; a hard-pinned centre 0.414. That is
  what `positions_prior` implements, and it exists to reproduce other aligners' conventions —
  `mir.distances.aligner.JunctionAligner` hardcodes `(3, 4, -4, -3)` for all seven loci — not
  because it is the right rule.

- **A per-island PWM beats min-over-members only at a strict cutoff, and "strict" depends on `N`.**
  The E-value's `k = floor(E*·M/N)` is how many control neighbours the cutoff admits, so FPR = `k/M`.
  108 calibrated VDJdb islands >= 10, held-out members, 250k control negatives, paired bootstrap
  over islands:

  | regime | FPR | min | PWM | diff [95% CI] |
  |---|---|---|---|---|
  | loose reference | 1% | **99.5%** | 99.1% | −0.40 [−1.09, +0.14] |
  | per-epitope islands (`N` = group size, median 88) | 0.0568% | 88.3% | **89.3%** | +0.93 [−0.80, +2.79] |
  | repertoire annotation (`N` ≈ 20k) | 0.0012% | 37.6% | **48.5%** | +10.90 [+7.69, +14.21] |

  No significant difference *while building* islands; a large one when annotating a repertoire with
  them (islands >= 50: 9.8% vs 22.6%). Note `E*=0.05` at `N≈20k` gives `k=0` → `thetas_from_scores`
  returns `-1`; the rule of three forbids any `E < 3N/M = 0.236`, and that is the third row's cutoff.
  Neither representation generalises to same-epitope junctions in a *different* island (3.7% vs 3.5%
  at 1% FPR, ~0% at either operating point).

### Batch scoring — `score_matrix`

A search prunes; an embedding cannot, because the distance to every prototype *is* the output.
Measured on an M3, 3,000 refs, human TRB: pure-Python `gapblock_score` **0.41 M pairs/s**,
`score_matrix` **51.3 M** (1 thread) and **532.7 M** (16 threads). The prior costs nothing — it
is flattened once into an `[m][d][i]` cube. Reproduce with `bench/bench_score_matrix.py`.

The result *is* the distance: the Gram transform is applied per residue when the matrix is built,
so there is no `d = s(a,a) + s(b,b) - 2·s(a,b)` step, and non-negativity, symmetry and a zero
diagonal hold by construction. Budget `4 * n * K` bytes and chunk the queries.

### Performance

Profiled over the bundled 250k control: 91% of query time is the query-deletion-variant branch,
9% the auxiliary indices (9.8M entries at `d_max=3`). Netting the prior out of each variant's
budget cuts that branch from ~15 sub-searches to 2.5. **Variant dedup and length-bucketing are
not worth building** — measured, not guessed. At budget 40 over 250k refs, `d_max=2` gap-block
search costs 2,562 µs/query, *less* than the plain Hamming ball at the same budget.

## Seed E-values — `seqtree.seeds`

```python
core_kmers(seq, k, flank=4) -> set[str]                      # drops the germline-framed ends
SeedIndex(seqs, k=5, flank=4) / SeedIndex.from_index(index, k, flank)
    .count(seed) / .evalue(seed, n_target) / .seed_evalues(query, n_target)
    .significant(query, n_target, alpha=1.0) / .union_evalue(...) / .gather(query, seeds=None)
```

Seeds buy **precision, not recall**: a shared rare central k-mer is real evidence (≈4× enriched),
but only ~0.5% of co-specific cross-island pairs share one. Never assume a shared k-mer is
significant — a D-gene run like `LAGG` sits in 3,471/250,000 control sequences (`E_seed` = 1,388).
Count it in the control.

## Layout / pMHC — `seqtree.layout`, `seqtree.pmhc`

`AnchorSpec`, `mask_anchors`, `kmers`, `presentation_features`, `weight_profile` (a
*mismatch-tolerance* profile, not information weighting), `PMHCStore`, `find_mimics`.

## Gotchas

1. `gap_open` defaults to **1**, which is wrong for any real matrix. Pass `2 * matrix.scale()`.
2. `Index.build` rejects non-alphabet characters — filter `B/Z/X/*` out of CDR3s first.
3. The bundled control is a **uniform reservoir sample** of unique productive clonotypes, shuffled,
   so `bundled[:k]` is a valid sub-sample. It used to be the abundance head (25.8x more
   self-similar, 3.1x the ball mass); if you have a cached `~/.cache/seqtree/control_*.sqtree`
   from before that fix, delete it.
4. OLGA generates from `P_gen`, which is **not** the post-selection background `P_0`. Do not use
   an OLGA sample as the null for an operational cutoff.
5. Rebuild after pulling: a stale `_core.so` has silently shadowed new Python-visible C++ methods
   more than once (`pip install -e .`).

## Layering

seqtree stays generic and dependency-free. TCR germline data, clustering, structures and
benchmarks against VDJdb live downstream in `vdjmatch` / `tcren-ms`. Do not import domain data
into `bench/`.
