# Downstream wrappers roadmap: `vdjmatch` and `mhcmatch`

**Status:** living draft. Owner: @mikessh. Extend freely — sections marked _(TBD: owner)_ await
more detail.

**Purpose.** `seqtree` is the shared substrate — a fast, payload-agnostic fuzzy-search core (C++ +
Python) with a control-calibrated E-value theory (see `appendix/evalue.tex`). It deliberately stops
at the reference-implementation level: the *applied* tools — `vdjmatch` (TCR antigen specificity) and
`mhcmatch` (peptide–MHC) — are separate packages built on top of it. This file is the contract
between them, written so an agent developing either package can pick it up cold: it states what
`seqtree` already provides, where each wrapper plugs in, and which design decisions are still open.

**How to use this doc (for agents).** Before implementing a feature in `vdjmatch`/`mhcmatch`:
1. Check **§1 Substrate** for the `seqtree` primitive that already covers it — do not reimplement
   search, E-values, anchor masking, or k-mer indexing.
2. Read the matching subsection in **§2** (vdjmatch) or **§3** (mhcmatch) for the intended design and
   the relevant `appendix/evalue.tex` section.
3. Treat the math in the appendix as the spec; treat formulas here marked _(design)_ as direction,
   not final — confirm with the owner before committing to one.

---

## 1. Substrate: what `seqtree` provides

| Capability | Module / symbol | Notes |
|---|---|---|
| Fuzzy fixed-length search | `seqtree.Index`, engines `seqtm` / `seqtrie` | seqtm: per-type edit caps + Hamming fast path. seqtrie: budget-only full-width DP over the trie (not banded), ignores per-type caps |
| Single-gap-block loop alignment | `seqtree.gapblock` (Python) | anchored CDR3/junction alignment: one contiguous indel + positional gap prior |
| Substitution scoring | `SubstitutionMatrix` (Gram → squared-distance penalty `s_aa+s_bb−2·s_ab`) | payload-agnostic |
| Per-position scoring | `PositionalMatrix` (`penalty(pos,a,b)` = base × per-position weight; weight 0 = free/anchor) | the hook for PSSMs and anchor masking |
| k-mer seed index | `KmerIndex` (C++): unique-k-mer trie + CSR postings + per-peptide allele tag; `seed_and_gather` (GIL-released, parallel) | million-scale candidate generation |
| Control-calibrated E-values | `seqtree.evalues`, `seqtree.load_control` | `Ê = (N/M)·n_C`, Poisson tail, `exclude_exact` |
| Calibrated score cutoffs | `seqtree.threshold_for_evalue`, `thetas_from_scores` | inverts `Ê` to a **per-query** θ. A fixed θ is not calibrated: at `gapblock_score ≤ 60`, random control junctions cluster *harder* than real same-epitope ones |
| Gap-placement priors | `central_prior`, `profile_prior`, `frame_prior`, `embed_in_frame` | `prior(i, d, m) ≥ 0` and `= 0` at `d = 0`. Only a **constant-`i`** rule (`frame_prior`) makes a frame transitive — and hence a column index, and hence a PWM, possible |
| Seed significance | `seqtree.seeds`: `core_kmers`, `SeedIndex` | precision, not recall: ~0.5% cross-island coverage |
| Anchor / layout model | `seqtree.layout`: `AnchorSpec`, `DEFAULTS`, `mask_anchors`, `kmers`, `presentation_features`, `weight_profile` | parametrized anchors; class-II register trick |
| pMHC homology + reverse | `seqtree.pmhc`: `PMHCStore`, `search_homologs`, `assign_allele`, `find_mimics`, `build_kmer_index` | reference impl; mhcmatch productionizes |
| Presentation-aware E-values | `seqtree.pmhc_evalue.homolog_evalue` | per-allele null |

**E-value theory** lives in `appendix/evalue.tex`: the empirical-control null (§Setup, §Null), the
Poisson/Chen–Stein bound (§Poisson), multiple testing (§E-value), the closest-hit Gumbel law
(§Gumbel), Karlin–Altschul as the i.i.d. special case (§KA), the pMHC presentation-aware extension
(§Epitopes, §Reverse problem), and elementary applications — UMI/barcode birthday-collision and
CDR3-nt error clustering (§Related applications).

**Datasets.** `isalgo/pmhc_data` ships two tiers (see `appendix` Table "Scope of the two pmhc_data
tiers"): **full** (every IEDB-positive epitope–allele assay) and **shortlist** (epitope–allele pairs
with ≥2 supporting publications). `mhcmatch` should expose tier choice; the benchmark uses full.
TCR side: VDJdb (specificity-labelled) + `isalgo/airr_control` (matched control repertoire).

---

## 2. `vdjmatch` — TCR antigen-specificity

**Goal.** Given a query TCR (CDR3, optionally V/J, optionally paired α/β), score its similarity to
known antigen-specific TCRs and return a control-calibrated E-value per candidate epitope —
generalizing TCRNET-style neighbour counting (appendix §Intro) to a usable annotation tool.

### 2.1 Train / test / validation splits
- **Split by epitope, not by sequence.** Held-out *epitopes* test generalization to unseen
  specificities; held-out *sequences within an epitope* test within-specificity recall. Never let the
  same clonotype leak across splits (dedup to unique clonotypes first — appendix Prop. on collapse).
- **Three-way:** train (fit substitution/positional weights, §2.2), validation (tune scope/budget and
  E-value thresholds), test (report ROC/PR per epitope, as the MHC-guess benchmark does per allele).
- **Null/control** is independent of the split: `airr_control` supplies the background `P₀`; size it
  per appendix §"How large must the control be?".
- _(TBD: owner)_ exact epitope-level CV scheme, minimum cluster size, negative sampling.

### 2.2 Custom substitution matrices (epitope–paratope interaction-aware)
- **Hook:** `seqtree.PositionalMatrix` (per-position `penalty(pos,a,b)`) and/or a bespoke
  `SubstitutionMatrix`. seqtree already wires `PositionalMatrix` into the seqtm Hamming path.
- **Intent:** replace BLOSUM with a CDR3-position-weighted matrix reflecting which residues actually
  contact the peptide–MHC (paratope), learned from structural contacts or from specificity data
  (residues whose substitution most changes specificity get the most weight; framework positions →
  weight 0). Mirrors how `layout.weight_profile` builds anchor/TCR-facing weights for pMHC.
- _(TBD: owner)_ source of contact statistics (structures vs. learned), per-V-gene vs. global matrix.

### 2.3 Single-chain E-value (have today)
- Per chain (β, or α), this is exactly `seqtree.evalues`: count neighbours of the query among the
  epitope-specific set vs. the control, report `Ê` and `p_enrich`. No new theory needed.

### 2.4 Paired α/β E-value (both chains known)
- _(design)_ Under chain independence in the null, the **joint** null ball-mass factorizes:
  `π₀^{αβ} ≈ π₀^α · π₀^β`, so the paired intensity is `λ_αβ = N · π₀^α π₀^β` and the paired enrichment
  is a Poisson tail on the count of references matching **both** chains within budget. Equivalent to
  Fisher-combining the per-chain `p_enrich` when independence holds; the appendix `b₂` co-occupancy
  term measures the dependence to correct for. Implement as a joint `seed_and_gather` keyed on a
  paired peptide id.
- This is strictly more specific than either chain alone (joint null mass is tiny) → far smaller
  E-values for true pairs.

### 2.5 Paired estimate for a single chain (rarity of the unknown chain)
- _(design)_ When only one chain (say β) is observed, marginalize the paired E-value over the unknown
  α: weight the single-chain evidence by the **rarity** of the partner that would complete a known
  pair — a generation-probability / abundance prior on α (`Pgen`, appendix §thymic/selection). A β
  matching a rare-α-paired reference is more informative than one matching a common-α pairing. Yields
  a "paired-equivalent" E-value from one chain without observing the other.
- _(TBD: owner)_ the exact prior (OLGA/Murugan `Pgen` vs. empirical α-abundance), and calibration.

### 2.6 API sketch _(non-binding)_
```python
vm = vdjmatch.Annotator.from_vdjdb(chains=("beta",), matrix=paratope_matrix, control=airr_control)
hits = vm.annotate(cdr3b="CASS...", v="TRBV...", scope=2)          # -> [(epitope, Ê, p_enrich), ...]
hits = vm.annotate_paired(cdr3a=..., cdr3b=...)                    # §2.4
hits = vm.annotate(cdr3b=..., partner_prior="pgen")               # §2.5
```

---

## 3. `mhcmatch` — peptide–MHC

**Goal.** Productionize the `seqtree.pmhc` reference layer: epitope homology / cross-reactivity,
molecular mimicry, allele guessing, and non-binder filtering — with tuned thresholds and the
additions below.

### 3.1 Reuse from seqtree (have today)
- Homology / mimicry: `search_homologs`, `find_mimics` (anchor-masked TCR-facing k-mers; per-allele
  presentation-aware E-values via `pmhc_evalue.homolog_evalue`). Positive control: the Dolton et al.
  A\*02:01 cross-reactive trio.
- Allele guessing (reverse problem): `assign_allele` + the vote-fraction ranking / register trick
  validated in `bench/bench_mhc_guess.py` (appendix §Reverse problem). ROC-AUC 0.90–0.98.
- **Non-binder filter** (appendix §Reverse problem, "Filtering non-binders"): _binds no MHC_ →
  best-over-panel E-value high / no allele scores above background; _doesn't bind allele a_ →
  per-allele `E_a > α`. mhcmatch exposes both thresholds.
- Tier choice (full vs shortlist, §1 Datasets).

### 3.2 Motif logos
- For each allele (and each class-II core register), render a sequence logo of the anchor / pocket
  residues from its presented set — the visual counterpart of `presentation_features`. Use a standard
  information-content logo (bits per position). Class II: build the logo over the register-aligned
  9-mer cores (the register trick already picks the core), so anchors P1/P4/P6/P9 are columns.
- _(TBD: owner)_ logo library choice; whether to weight by `n_references` (shortlist confidence).

### 3.3 MHC pseudosequence → clustering & promiscuity
- **Pseudosequence:** represent each allele by the polymorphic residues that line the peptide-binding
  groove (the pseudosequence of peptide-contacting positions, as used for per-allele MHC binding
  prediction by Glynn, Ghersi & Singh, *PNAS* 2025, [doi:10.1073/pnas.2405106122](https://doi.org/10.1073/pnas.2405106122)).
  This turns "allele" from an opaque label into a short sequence the **same seqtree engine** can compare.
- **Cluster MHCs** by pseudosequence distance (seqtree fuzzy search on the pseudosequence "alphabet"
  of groove residues) → groups of functionally similar alleles. This is the principled way to express
  **cross-allele similarity**, which `seqtree.pmhc` deliberately does *not* model (appendix §Impl.
  limitation: "distinct alleles are distinct nulls"). mhcmatch is where that limitation is lifted.
- **Promiscuity** (esp. class II — appendix §"Class-II promiscuity"): a peptide presented across a
  pseudosequence-cluster of alleles is promiscuous; quantify as the spread of its positive alleles
  over the MHC clustering, and use the cluster to pool nulls for related alleles when data are thin.
- _(TBD: owner)_ pseudosequence position set per locus (HLA-A/B/C, DR/DQ/DP, mouse H-2), distance
  metric, cluster cut.

### 3.4 API sketch _(non-binding)_
```python
mm = mhcmatch.Store.from_pmhc(tier="shortlist")
mm.logo("HLA-A*02:01")                       # §3.2
clusters = mm.cluster_alleles(class_="mhc1")  # §3.3, via pseudosequence
mm.is_binder("SIINFEKL", allele="H-2Kb", alpha=0.05)   # §3.1 non-binder filter
mm.promiscuity("PKYVKQNTLKLAT")              # §3.3
```

---

## 4. Shared conventions
- **seqtree is upstream and stays generic.** New general-purpose primitives (a learned matrix loader,
  a logo helper, a pseudosequence comparator) may land in seqtree if reusable; tuned thresholds,
  predictors, and domain glue stay in the wrappers.
- **Anchors** are parametrized in `seqtree.layout` (presets per class, overridable) — wrappers pass
  `AnchorSpec`, they don't hardcode positions.
- **Citations:** never fabricate. Verify every DOI via a tool before adding it (PubMed/arXiv).
- **Versioning / gitflow:** feature branch → `dev` → `master`; end commit messages with the
  `Co-Authored-By` trailer; don't publish to PyPI without an explicit release.

## 5. Pointers
- E-value theory & all derivations: `appendix/evalue.tex` (compiled `appendix/evalue.pdf`).
- pMHC usage & limitations: `docs/pmhc.rst`. Benchmarks & figures: `docs/benchmarks.rst`,
  `bench/bench_mhc_guess.py`.
- seqtree internal roadmap (PSSM-graded d_TCR, native local align, Flashback build, predictor proteome
  scans): `docs/roadmap.rst`.
