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
SubstitutionMatrix.blosum62() / .pam250() / .pam100() / .structural() / .unit(n)
SubstitutionMatrix.from_similarity(grid: list[list[int]])   # Gram: pen = s_aa + s_bb - 2*s_ab, clamped >= 0
SubstitutionMatrix.penalty(a, b) -> int                      # pen(a, a) == 0 always
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
```

**A fixed score cutoff is not a calibrated cutoff.** The control is dense near germline and
sparse among rare junctions. Measured on human TRB: at `gapblock_score <= 60`, random control
junctions cluster *harder* than real same-epitope ones. Always invert an E-value; never pick θ
by hand. `-1` means `e_target < 3N/M` — enlarge the control rather than lowering the bar.

`threshold_for_evalue` needs `params.max_penalty > 0` as the ceiling to search the control at.
One control scan supplies every query's cutoff.

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
embed_in_frame(seq, width, c, gap="-") -> str
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
3. The bundled control has detectable row-order structure. Shuffle before subsampling; using all
   250k rows is order-invariant and safe.
4. OLGA generates from `P_gen`, which is **not** the post-selection background `P_0`. Do not use
   an OLGA sample as the null for an operational cutoff.
5. Rebuild after pulling: a stale `_core.so` has silently shadowed new Python-visible C++ methods
   more than once (`pip install -e .`).

## Layering

seqtree stays generic and dependency-free. TCR germline data, clustering, structures and
benchmarks against VDJdb live downstream in `vdjmatch` / `tcren-ms`. Do not import domain data
into `bench/`.
