# `TextIndex` — exact k-mismatch search over a concatenated reference text

**Status: implemented in `seqtree` 1.0.0, 2026-09-06.** Steps 1-6 of §8 shipped; §3.7 gapped
q-grams and §9's `mhcmatch` adapter did not. Six things below were corrected during
implementation and the code, not this document, is authoritative on them:

1. **The amino-acid codec is 24 symbols** (`ARNDCQEGHILKMFPSTWYVBZX*`), not 20. So `A^k` at
   `k = 4` is 331,776 buckets (§5.2's 160,001 assumed 20), and the ball's branching factor is
   23, not 19.
2. **The ball is enumerated over `q[0:k]`, not over the whole query** (§5.3). The two give the
   *same* candidate set -- the distinct leading k-mers of the full-query ball are exactly the
   ball of `q[0:k]` -- at 3,267 probes instead of 19,252 for `k = 4, m = 2`, and independently
   of `L`.
3. **§5.3's "reuses 0.7.0's neighbourhood work" does not hold.** That work is pure Python with
   a BFS and a shared `seen` set, whose dedup exists only because it unions over many centres.
   The C++ ball is new code; what transferred was the duplicate-free guarantee.
4. **Queries shorter than `k` are refused.** Seeds are indexed only where a whole k-mer fits
   inside a record, so a shorter query could match within `k - L` of a record end and be missed.
5. **An out-of-alphabet residue in the TEXT is a counted hole, not an error.** The human
   proteome has 36 `U`; refusing the build over them would make the class unusable on its own
   reference data. A *query* containing one is still refused, as §7.2 requires.

6. **§5.3's two-path dispatch is one rule, and that is where the speed was.** As written, a
   query took `m + 1` exact seeds when `L / (m+1) >= k` and otherwise a single radius-`m` ball
   over `q[0:k]` -- with nothing in between, so `L = 8, m = 2, k = 4`, where two disjoint 4-mers
   plainly fit, still paid for a 3,267-variant one-block ball. The shipped rule generalises both:

   > Split into `b` disjoint blocks of width `bw = L/b >= k`, probe block *j* at radius `c_j`.
   > **Lossless iff `sum_j c_j >= m - b + 1`** -- the cheapest violating error vector is
   > `e_j = c_j + 1`, of weight `sum(c_j) + b`, so none exists below that. Probe cost
   > `N(c) = sum_{i<=c} C(k,i)(A-1)^i` has steeply rising, block-independent increments, so the
   > cheapest legal scheme is `b = min(m+1, L/k)` with `r = max(0, m-b+1)` spread evenly.

   `b = m+1` recovers the seed path and `b = 1` the ball path. Measured on the human proteome,
   single-threaded, `k = 4` (`bench/bench_text_index.py`):

   | `L` | `m` | before, ms/query | after, ms/query | |
   |--:|--:|--:|--:|--:|
   | 8 | 2 | 34.505 | 1.379 | 25x |
   | 9 | 2 | 35.453 | 1.278 | 28x |
   | 10 | 2 | 40.940 | 1.311 | 31x |
   | 11 | 2 | 30.196 | 1.242 | 24x |
   | 12 | 3 | 415.578 | 1.804 | 230x |
   | 15 | 3 | 400.484 | 1.531 | 262x |

   Two things fell out of it that the document did not anticipate. **Spare budget units belong on
   the LAST blocks** -- worth a further 1.4x, because a candidate matches its own block's k-mer,
   so verification scanning left-to-right meets the unconstrained prefix first and early-exits,
   and it is the high-budget block that contributes nearly all the candidates. And **§5.4.5's
   candidate sort is not needed at all**: a start two blocks both reach is emitted by the
   lowest-indexed block that could have produced it, a test read off the mismatch positions
   verification has already computed.

   The framing is the "search scheme" of Kianfar, Pockrandt, Torkamandi, Luo & Reinert,
   *Optimum Search Schemes for Approximate String Matching Using Bidirectional FM-Index*,
   arXiv:1711.02035. Only the partition-and-budget half transfers: a direct-addressed table has
   no bidirectional extension, so there is no search *order* to optimise over.

§5.4's bit-packed verification was measured as the wrong lever and not built: the short path is
bound by one cache miss per candidate, not by the comparison loop, and §5.4.5's locality sort
made it **1.6x slower** when tried. Against §7.1's acceptance targets, `L >= 12, m = 2` is met
with room (0.127 ms/query against <= 0.5) and `L = 8-11, m = 2` lands at **1.24-1.38 ms against
<= 1.0** -- 1.3x over, from 30x over. At `k = 5` every `L >= 10` is under 0.2 ms; see
`docs/text-index.rst` for the measured `k` guidance.

Driven by a measured failure in `mhcmatch`, but the primitive is generic and stays generic
(`ROADMAP.md` §4, "seqtree is upstream and stays generic"). The reference implementation to beat is
**PEPMatch** (§3); this document sets out how, on seven axes, with the numbers for each.

---

## 1. Problem scope

### 1.1 The question

Given a **reference text** — a proteome, 10⁷–10⁸ residues in ~10⁵ records — and a **large batch of
short queries of many different lengths**, find **every** position where a query matches within `m`
substitutions. No gaps. Full length. **Exact answer, not a heuristic one.**

This is peptide-to-proteome mapping. It is not homology search: no scoring matrix in the predicate,
no E-value, no local alignment. It is Hamming-ball search over a text.

### 1.2 Why `Index` is the wrong structure for it

`Index::build(std::vector<std::string> refs, Alphabet)` builds a trie over reference **strings**. To
ask a text question with it the caller must enumerate every length-`L` window of the text as a
separate string — and do it again for every distinct query length. `mhcmatch.Proteome._index(L)`
does exactly that, and its own docstring records the cost:

> The human proteome has **68,389,335** 9-mer windows… 12.6 GB peak for the first length and ~3.6 GB
> for each further one. **Ask for the lengths you need** — `find_sources` builds one index per
> distinct query length.

Right for a fixed-length class-I register. Collapses the moment the query set is not fixed-length.

### 1.3 The measurement that prompted this

`bench/neoag/resolve_genes.py` in `2026-mhcmatch-benchmark`, over `bench/epic/neoantigens.parquet`
(886,047 rows), 2026-09-06:

| | |
|---|--:|
| peptides needing a source lookup | 445,466 human + 2,310 mouse |
| **distinct human query lengths** | **45** |
| index size, one human length | 3.7 – 9.7 GB |
| on-disk index cache after one *partial* run | **225 GB** |
| wall clock, human side, **never finished** | **> 2 h 10 m** |

Killed at 2 h 10 m with the disk at 99 % full. The length histogram is the whole story:

| lengths | peptides | share | what they are |
|---|--:|--:|---|
| L = 8 – 11 | 347,813 | 78.1 % | class-I registers |
| L = 12 | 83,517 | 18.7 % | class-II core + 3 |
| L = 13 – 50, **40 distinct lengths** | 14,136 | **3.2 %** | class-II, minigenes, constructs |

L = 43 has **one** peptide. L = 39 has three. L = 33 has six. Each triggered a multi-gigabyte
whole-proteome index build. **The cost is not the search — the index is keyed on query length when
the reference is a text that does not have one.**

### 1.4 The workarounds already in the tree, each a symptom

- `bench/mhc2_human/mmseqs_sources.py`: *"Class II admits 11-25, so that is fifteen indexes at
  ~5.5 GB each — **~82 GB measured**."* Routes class II through mmseqs2 instead. **See §3.4 — that
  substitution is probably losing most of its matches.**
- `bench/cassette/tesla_genes.py`: *"**12.6 GB peak** `find_sources` needs for the first one."*
- `mhcmatch.proteome.index_cache_dir`: an entire disk-cache subsystem with cross-process build
  hand-off (`O_EXCL` claim, bounded wait) exists *because* the builds are so expensive. A cheap index
  removes the reason for it.

### 1.5 Callers

`mhcmatch`: `Proteome.find_sources` / `find_source` / `find_exact_sources`, `vector.self_origin_risk`
(`vector.py:1107`, `max_subs=1`), `cli.py:483` (`mhcmatch sources`).
`2026-mhcmatch-benchmark`: `bench/neoag/resolve_genes.py`, `bench/mhc2_human/gene_recovery.py`,
`bench/cassette/tesla_genes.py`, and it retires `bench/mhc2_human/mmseqs_sources.py`.

---

## 2. Measured prototype

Python/NumPy, **single-threaded**, author's laptop (Apple M-series), human proteome
(`~/hf/pmhc_data/proteome/human.fasta.gz`, 147,506 records, **69,578,135 residues**), `max_subs = 2`.
Measurements, not estimates.

| approach | build | per query | verdict |
|---|--:|--:|---|
| naive vectorised Hamming scan, L = 9 | — | **1161 ms** | ruled out |
| naive vectorised Hamming scan, L = 31 | — | 3199 ms | ruled out |
| pigeonhole, blocks of ⌊L/3⌋ = 3, L = 9 | 25 s | 6.83 ms | 3-mers not specific enough |
| **shared k = 4 seed table, L = 12** | **42 s once** | **0.47 ms** | |
| **same table, L = 15 / 20 / 25 / 31 / 50** | *(reused)* | **0.46 – 0.53 ms** | **flat in L** |
| **exact L-mer ball, L = 9** | 49 s | **1.96 ms** | 13,168 codes/query |
| **exact L-mer ball, L = 11** | 52 s | **2.74 ms** | 20,065 codes/query |

Mean candidates surviving the k = 4 gather: **~3,200–3,600 per query**, near-constant L = 12…50.

Projected on the corpus of §1.3: L ≥ 12 (97,653 peptides) in **~90 s** from one table; L = 8–11
(347,813) in **~16 min**; **human total ≈ 17 min at ~1.1 GB peak**, against **> 2 h unfinished at
225 GB**.

### 2.1 Why the crossover is structural, not tuning

Pigeonhole at `m = 2` needs three disjoint blocks, so the seed is ⌊L/3⌋ residues:

| L | seed width | buckets (20 letters) | positions/bucket | candidates/query |
|--:|--:|--:|--:|--:|
| 9 | 3 | 8,000 | ~8,700 | ~13,000 |
| 12 | 4 | 160,000 | ~435 | ~3,500 |
| 31 | 10 | 20¹⁰ | ≪ 1 | ~3,400 |

**Short queries are the hard case.** They get their own path, and it is where most of the corpus is.

### 2.2 Where the short path actually spends its time — the profile that sets the design

| phase | ms/query | share |
|---|--:|--:|
| **ball enumeration** (incl. `np.unique`) | **1.039** | **65.4 %** |
| 2× `searchsorted` (`lo` and `hi`) | 0.548 | 34.6 % |
| 1× `searchsorted` alone | 0.330 | |

**Two-thirds of the cost is generating the ball, and nearly all of that is a sort that does not need
to happen.** Under *proper* substitutions (target residue ≠ original), the ≤ m ball is
**duplicate-free by construction**:

```
|B_m(q)| = Σ_{i=0..m} C(L, i) · 19^i
L = 9, m = 2  ->  1 + 9·19 + 36·361 = 1 + 171 + 12,996 = 13,168
```

which is **exactly** the 13,168 measured. The `np.unique` finds nothing to remove. In C++ the ball
is a direct write into a preallocated `13,168 × uint64` buffer: no sort, no hash set, no allocation
per query. And only `lo` is needed if bucket counts are stored, halving the lookup too.

---

## 3. The bar: PEPMatch, and the seven axes to beat it on

> Marrama D, Chronister WD, Westernberg L, Vita R, Koşaloğlu-Yalçın Z, Sette A, Nielsen M,
> Greenbaum JA, Peters B. **PEPMatch: a tool to identify short peptide sequence matches in large
> sets of proteins.** *BMC Bioinformatics* 2023;24(1):485.
> PMID [38110863](https://pubmed.ncbi.nlm.nih.gov/38110863/) ·
> doi [10.1186/s12859-023-05606-4](https://doi.org/10.1186/s12859-023-05606-4) ·
> <https://github.com/IEDB/PEPMatch>

PEPMatch is the right tool to measure against: same problem, same community, an explicit benchmark
framework, and 100 % recall. Its published numbers, from the paper (16-core i9-9900K, 32 GiB, NVMe):

| benchmark | queries | `m` | preprocess | **search** | per query | recall |
|---|--:|--:|--:|--:|--:|--:|
| exact, human proteome | 2,000 × 9-mer | 0 | 39.6 s | **0.08 s** | 0.04 ms | 100 % |
| neoepitope, human proteome | 620 × 15-mer | ≤ 3 | — | **18.4 s** | **29.7 ms** | 100 % |
| SARS-CoV-2, betacoronavirus | 628 × 8–15-mer | ≤ 2 | — | (0.77–115.9 s range) | | 100 % |
| milk allergen, best match | 677 × 15-mer | best | 203.7 s (NmerMatch) | | | 100 % |

Its exact-match search is excellent and hard to beat. **Its mismatch search is the target: 29.7 ms
per 15-mer at m = 3.** Our k = 4 seed path measures **0.47–0.53 ms** per query at m = 2 in
single-threaded NumPy. Different `m`, so not a like-for-like claim — §7 specifies the like-for-like
run — but the gap is two orders of magnitude and worth chasing.

### 3.1 One index for every length, not one per k

PEPMatch, from the paper: *"This preprocessing step is performed only once per proteome **and per
given k value**"*, with `k = ⌊L/(m+1)⌋`. So a query set spanning 45 lengths spans many `k`, and each
needs its own preprocessed database. The comparator NmerMatch is worse and the paper says so —
*"this step increased significantly for NmerMatch because it performs this step for each length of
the query peptides"*, **214 s** of preprocessing on one benchmark. That is precisely our §1.3
failure, in their table.

**Ours:** `k` is a property of the *index*, not the query. One `k = 4` table serves every `L ≥ 12`
(measured flat, §2), and `L < 12` dispatches to the ball path against the same table. **One build,
all lengths, all `m`.**

### 3.2 Best match without log₂(L) more databases

PEPMatch's best-match mode *"perform[s] the preprocessing step multiple times on a proteome for
different k values"* — `k = L, L/2, L/4, … 2` — then searches ascending. That is ~log₂(L) whole
proteome databases.

**Ours:** one index, ascending `m` with early exit per query. A query resolved at `m = 0` never
enters the `m = 1` pass. No extra structure at all.

### 3.3 Threads, which PEPMatch does not have

The paper's own future work: *"The utilization of GPU programming and **parallelization** may
significantly improve performance."* PEPMatch is single-threaded Python over SQLite/pickle.

**Ours:** C++, immutable index, one worker per thread, GIL released — the shape `Index::search_batch`
already has. On the 8–16 cores of a typical workstation this is a free order of magnitude *on top of*
the algorithmic difference, and it is the difference between 16 min and ~2 min on §1.3's corpus.

### 3.4 Recall — and a warning for `mhcmatch` today

PEPMatch's benchmark measured the alignment tools' recall on short peptides:

| method | exact 9-mers | 8–15-mers, m ≤ 2 | 15-mers, m ≤ 3 |
|---|--:|--:|--:|
| PEPMatch / NmerMatch | 100 % | 100 % | 100 % |
| BLAST | 98.3 % | 73.3 % | 58.1 % |
| DIAMOND | < 2 % | 6.4 % | 34.0 % |
| **MMseqs2** | **< 2 %** | **7.4 %** | **24.6 %** |

**`mhcmatch`'s `bench/mhc2_human/mmseqs_sources.py` routes class-II source lookup through MMseqs2.**
Its `-s 7.5` / `--min-ungapped-score 0` / permissive `-e` and its post-hoc predicate filter will do
better than the paper's `-s 7` defaults, but every one of those knobs sits *downstream* of a k-mer
prefilter that was measured at 7.4 % recall on exactly this length range. **This is worth measuring
on our own data before it is trusted, independent of whether `TextIndex` gets built** — and
`TextIndex` retires the script either way.

### 3.5 Storage: flat and mmap-able, not SQLite/pickle

PEPMatch stores the k-mer map *"in a SQLite database or a serialized pickle format"*. Ours is two
flat `uint32` arrays plus the text (§5.2) — `mmap`-able, zero-copy to load, shareable read-only
across processes. That is what makes 39.6 s of preprocessing a one-off rather than a per-process
cost, and it removes `mhcmatch`'s reason for its `O_EXCL` build hand-off.

### 3.6 Output that answers the next question too

PEPMatch reports *"the original sequence, the matched sequence, the UniProt protein ID, the protein
name, the gene symbol, the number of mismatches, and the residue positions where those mismatches
occur"*. Good, and the paper names what is missing: *"Another potential extension for the tool might
be **ranking matches based on an amino acid substitution score**. Certain amino acid substitutions
are more frequent than others"* — and *"peptides will often map to multiple proteins, especially
when accounting for residue substitutions, which is important in immunology."*

Both are cheap for us and neither exists there. See §6.

### 3.7 Better filtration than contiguous blocks

PEPMatch's blocks are contiguous and non-overlapping. **Gapped q-grams** give strictly better
lossless filtration for the k-mismatch problem at the same seed weight:

> Burkhardt S, Kärkkäinen J. **Better filtering with gapped q-grams.** *Fundamenta Informaticae*
> 2003;56(1–2):51–70. Preliminary version: *Combinatorial Pattern Matching (CPM 2001)*, LNCS 2089,
> doi [10.1007/3-540-48194-X_6](https://doi.org/10.1007/3-540-48194-X_6).

A gapped shape of weight `w` spanning `s > w` positions can be lossless for `m` mismatches where a
contiguous `w`-mer is not, because the gaps let one shape straddle mismatch positions. The paper's
own caveat is the reason this is rung 3 and not rung 1: *"the best shapes are rare and often possess
no apparent regularity"* — they are found by search, per `(L, m, w)`. **Specify it, ship without it,
measure it as an arm** (§8, step 7).

---

## 4. Reference data

All public, all fetchable. HuggingFace dataset **`isalgo/pmhc_data`**
(<https://huggingface.co/datasets/isalgo/pmhc_data>), mirrored at `~/hf/pmhc_data`.

| file | contents | size |
|---|---|--:|
| `proteome/fluA_PR8_UP000009255.fasta.gz` | ~5 kB of sequence — **brute-force ground truth** | 3.4 KB |
| `proteome/ecoli_K12_UP000000625.fasta.gz` | 1.4 M residues — size step | 988 KB |
| `proteome/mouse.fasta.gz` | UP000000589, 54,857 records / 23,131,234 residues | 15.0 MB |
| `proteome/human.fasta.gz` | UP000005640, 147,506 records / 69,578,135 residues | 38.5 MB |
| `neoantigens/neoag_tested.tsv.gz` | real queries with the §1.3 length profile | 7.1 MB |
| `immunogenicity/pathogen_tcell_mhc1_human.tsv.gz` | 16,841 class-I 8–11-mers | 417 KB |
| `immunogenicity/pathogen_tcell_mhc2_human.tsv.gz` | 7,948 class-II 12–20-mers | 202 KB |

```python
from huggingface_hub import hf_hub_download
p = hf_hub_download("isalgo/pmhc_data", "proteome/human.fasta.gz", repo_type="dataset")
```

PEPMatch's own benchmark sets (2,000 9-mers; 628 SARS-CoV-2; 620 neoepitope 15-mers; 677 milk
allergen 15-mers) are in `github.com/IEDB/PEPMatch` and are the like-for-like comparison (§7).

**Non-standard residues are present and are the accuracy trap (§6.2):**

| proteome | `X` | `U` | `B` | `Z` | total |
|---|--:|--:|--:|--:|--:|
| human | 8,417 | 36 | 0 | 0 | 8,453 |
| mouse | 4,263 | 33 | 17 | 23 | 4,336 |

---

## 5. Design

### 5.1 C++ surface

```cpp
// include/seqtree/text_index.hpp
#pragma once
#include "seqtree/seqtree.hpp"

namespace seqtree {

struct Mismatch { uint16_t pos; char query_aa; char text_aa; };

struct TextHit {
    uint32_t ref_id = 0;      // which input record
    uint32_t offset = 0;      // start within that record
    uint16_t n_subs = 0;      // Hamming distance; always <= max_subs
    int32_t  score  = 0;      // substitution-matrix score, 0 if none supplied
    // Mismatch detail lives in a parallel CSR (see TextResult) -- never one vector per hit.
};

struct TextQueryOpts {
    uint16_t max_subs      = 0;
    bool     exclude_exact = false;   // drop 0-mismatch hits (a candidate that IS the text)
    bool     best_only     = false;   // keep only the minimum-n_subs shell (still ALL of it)
    uint32_t max_hits      = 0;       // 0 = unlimited; a cap is reported, never silent
    const SubstitutionMatrix* matrix = nullptr;   // optional: score each hit
    int32_t  group_by      = -1;      // >=0: aggregate to group ids (see build)
};

class TextIndex {
public:
    // `refs` are whole records, NOT windows. `k` is the seed width (default 4, max 6).
    // `group_ids` is optional and arbitrary -- gene ids, species, cluster labels: seqtree
    // does not know what they mean, only that hits can be folded to them.
    static std::unique_ptr<TextIndex> build(const std::vector<std::string>& refs,
                                            Alphabet, uint8_t k = 4,
                                            const std::vector<uint32_t>& group_ids = {});
    ~TextIndex();

    uint32_t num_refs() const;  uint64_t num_residues() const;
    uint8_t  k() const;         Alphabet alphabet() const;

    // Flat CSR results -- see §5.5. threads <= 0 => hardware_concurrency.
    TextResult search_batch(const std::vector<std::string>& queries,
                            const TextQueryOpts&, int threads = 0) const;

    void save(const std::string& path) const;              // flat, mmap-able
    static std::unique_ptr<TextIndex> load(const std::string& path, bool mmap = true);
};
}  // namespace seqtree
```

### 5.2 Layout

```
text_       vector<uint8_t>    codec symbols; Codec::kInvalid between records
starts_     vector<uint64_t>   num_refs + 1; absolute offset of each record
group_      vector<uint32_t>   num_refs; caller labels (0 if unused)
post_begin_ vector<uint32_t>   A^k + 1   -- CSR bucket offsets, DIRECT-ADDRESSED
post_ids_   vector<uint32_t>   num_residues -- absolute positions grouped by k-mer
packed_     vector<uint64_t>   optional: 5-bit-packed text for SIMD verification (§5.4)
```

At `A = 20, k = 4`: **160,001 buckets, ~69.7 M postings ≈ 280 MB**, one allocation each, built by
**counting sort in two passes** — no comparison sort, no hashing, no `argsort`. (The NumPy prototype's
42 s build is dominated by an `argsort` this replaces.)

`Codec::kInvalid` between records is load-bearing twice: a k-mer containing it is never inserted, so
a **seed cannot span a record boundary**, and verification against it never matches, so a **hit
cannot span one either**. Structural, not a check.

Absolute position → `(ref_id, offset)` is one `std::upper_bound` on `starts_`.

### 5.3 Query dispatch

Per query, on `s = L / (max_subs + 1)`:

**`s >= k` → seed path.** Split into `max_subs + 1` disjoint blocks; probe each block's first `k`
residues; translate postings to candidate starts (`position − block_offset`); dedupe; verify.
Lossless by pigeonhole: a match within `max_subs` leaves ≥ 1 block untouched, so ≥ 1 seed is exact.

**`s < k` → ball path.** Enumerate the query's ≤ `max_subs` neighbourhood **without dedup or sort**
(§2.2 — it is duplicate-free by construction) into a preallocated buffer, and probe each variant's
first `k` residues against the *same* table, verifying as above. Reuses 0.7.0's neighbourhood work
(`d7a6e2b`) and adds no second structure.

**`best_only`** runs `m = 0, 1, … max_subs` and stops a query at the first `m` that yields a hit —
returning **all** hits at that `m`. One index, no extra preprocessing (contrast §3.2).

### 5.4 Constant factors worth specifying

1. **Bit-packed verification.** 20 letters fit in 5 bits; 12 residues per `uint64`. Verify by XOR of
   packed words, mask to 5-bit lanes, `popcount` of non-zero lanes. Early exit past `max_subs`.
   Replaces a byte loop with 1–4 word operations for L ≤ 48.
2. **Duplicate-free ball, no sort** (§2.2) — removes ~65 % of the measured short-path cost.
3. **Counting-sort build**, not `argsort` (§5.2).
4. **One-sided lookup:** store bucket counts so only `post_begin_[key]` is read, not `[key+1]`
   — the profile shows `hi` costs as much as `lo`.
5. **Sort queries by first-block bucket** before dispatch, so consecutive queries touch adjacent
   postings. Restore input order at the end. Pure locality win, no semantic change.
6. **`mmap` load** — a 280 MB index is paged in on demand and shared across processes.

### 5.5 Results: flat CSR, zero-copy to NumPy

Never a `vector<vector<Hit>>` with a `vector<Mismatch>` inside each hit. One allocation per array:

```cpp
struct TextResult {
    std::vector<uint32_t> query_begin;  // n_queries + 1, CSR into the hit arrays
    std::vector<uint32_t> ref_id, offset;
    std::vector<uint16_t> n_subs;
    std::vector<int32_t>  score;
    std::vector<uint32_t> mm_begin;     // n_hits + 1, CSR into the mismatch arrays
    std::vector<uint16_t> mm_pos;
    std::vector<char>     mm_query_aa, mm_text_aa;
    std::vector<uint8_t>  truncated;    // n_queries; 1 if max_hits capped this query
};
```

Every array becomes a zero-copy `numpy` view in the binding. A 445 k-query run returns a handful of
arrays, not millions of Python objects — which is most of the difference between "fast search" and
"fast tool".

---

## 6. Output — what to report beyond a position

### 6.1 Mismatch detail, and the substitution score PEPMatch names as future work

Per hit: `n_subs`, and per mismatch `(pos, query_aa, text_aa)`. The **pair**, not just the position:
a caller ranking by chemistry needs to know `L→I` from `L→D` without re-fetching the window.

With `TextQueryOpts::matrix` set, each hit also carries a `score` from an existing
`seqtree::SubstitutionMatrix`. `best_only` then has two orderings available — fewest substitutions,
or best-scoring — and the caller picks. This is exactly the extension the PEPMatch paper names and
does not have, and it costs one table lookup per mismatch on a hit we have already verified.

### 6.2 Group aggregation, with ties made explicit

`build(..., group_ids)` labels each record; `TextQueryOpts::group_by` folds hits to
`(group_id, min n_subs, n_hits)`. seqtree does not know a group is a gene.

This is worth having in the library rather than the caller because it is where the real work is:
`resolve_genes` currently pulls every hit into Python, takes `min(n_subs)`, maps protein → gene, and
**refuses the peptide if the nearest shell names more than one gene**. Folding in C++ turns
~13,000 candidate positions into a handful of group rows and makes the tie a first-class output
rather than something each caller re-derives. It also directly answers the PEPMatch paper's
*"peptides will often map to multiple proteins… which is important in immunology."*

### 6.3 Truncation is reported, never silent

`max_hits` exists because a degenerate query can match a repeat region a million times. When it
fires, `truncated[q] = 1`. A cap that is invisible is a recall bug wearing a performance costume.

---

## 7. Accuracy — what must not regress

The current path is **exact** and the caller's contract depends on it. `resolve_genes` refuses a
peptide whose nearest parents name more than one gene, on the stated grounds that *"`expr_norm` on
the wrong gene is worse than `expr_norm` missing"*. A missed hit does not degrade an answer there —
it **changes** it, from ambiguous to confidently wrong. Four obligations:

1. **Completeness is a proof obligation, not a benchmark.** Test against a **brute-force scan**, not
   against the old implementation: text = fluA (~5 kB, §4); every `L` from 6 to 30; every `max_subs`
   from 0 to 3; `k` ∈ {3, 4, 5}; queries = every text window, mutated copies, random strings, and
   strings containing `X`. Assert **set equality** of `(ref_id, offset, n_subs)` and of the mismatch
   detail — not merely that hits were found. Repeat on E. coli as a size step.
2. **The alphabet is where the silent wrong answer lives.** The NumPy prototype for §2 encoded
   residues as `ord(c) % 25`, under which **`A` (65) and `Z` (90) collide** — a false match no
   "hits were found" test would catch. The proteomes carry 8,453 human and 4,336 mouse non-standard
   residues (§4), so this is not hypothetical. Use `Codec::kInvalid`; assert a window containing it
   never matches, and that a *query* containing it is refused with a named error (0.6.1 set that
   precedent: *"gapblock alphabet error now names index/string"*).
3. **Return the whole nearest shell.** No truncation without `truncated` (§6.3), no early stop at
   the first hit. Tie detection is the caller's and needs every hit at the minimum distance.
4. **Deterministic order** — `(n_subs, ref_id, offset)`. Stable across runs and thread counts, so a
   downstream digest is stable.

### 7.1 The like-for-like benchmark

Run PEPMatch's own four benchmark datasets (§4) through both tools on one machine and report the
same metrics its framework does — preprocess time, search time, recall. Its framework accepts a
Python wrapper around any tool, so this is a wrapper, not a reimplementation. **Recall must be
100 % on all four**, or the work has failed regardless of speed.

Then our own axes:

| axis | values |
|---|---|
| text | fluA (5 kB), E. coli (1.4 M), mouse (23 M), human (69.6 M) |
| query length | 8, 9, 10, 11, 12, 15, 20, 25, 31, 50 |
| `max_subs` | 0, 1, 2, 3 |
| `k` | 3, 4, 5 |
| threads | 1, 8, all |
| report | build s, peak RSS, ms/query, mean candidates, hits, recall vs brute force |

**Acceptance targets** (single-threaded, human proteome, `m = 2`), from §2: ≤ 0.5 ms/query at
L ≥ 12 from one shared index; ≤ 1.0 ms/query at L = 8–11 (the ball path with §5.4.2 applied, against
1.96 ms measured in NumPy); ≤ 1.5 GB peak; build ≤ 20 s. And the end-to-end number that motivated
this: `resolve_genes` over 445,466 human peptides across 45 lengths, currently **> 2 h / 225 GB,
unfinished** — target **< 5 min / < 2 GB** on 8 threads.

---

## 8. Order of work

1. `TextIndex` skeleton — flat text, sentinel, `starts_`, CSR seed table by counting sort,
   `save`/`load` (mmap), accessors.
2. Seed path + byte-loop verification, single-threaded. **Brute-force equivalence test (§7.1) before
   anything else is added.**
3. Ball path, duplicate-free and unsorted (§2.2); same equivalence test at `L < k·(m+1)`.
4. `TextResult` CSR + mismatch detail + `best_only` + `group_by` + `truncated`.
5. Threading; bit-packed verification (§5.4.1); query-order locality (§5.4.5).
6. Python binding, `.pyi`, doctest, `docs/api.rst`; then the benchmark tables of §7.1.
7. **Arm, not default:** gapped q-grams (§3.7). Ship only if measured better on §7's grid.
8. **Only then** the `mhcmatch` adapter (§9), behind a keyword so the trie path stays reachable for
   one release and the two can be diffed on real data.

Gitflow per `ROADMAP.md` §4: `feature/text-index` → `dev` → `master`, `CHANGELOG.md` entry, 0.8.0.

---

## 9. Consuming it from `mhcmatch`

The point of the packaging is that `mhcmatch` gains a dependency and loses a subsystem.

### 9.1 Install

`seqtree` already builds wheels via `scikit-build-core` + `pybind11` and declares **no runtime
dependencies**. So:

```toml
# mhcmatch/pyproject.toml
dependencies = ["seqtree>=0.8.0", ...]
```

```bash
uv add "seqtree>=0.8.0"     # binary wheel; no CMake, no compiler, no C++ toolchain
```

Nothing else changes. `TextIndex` must therefore be exported from the top-level `seqtree` namespace
(`python/seqtree/__init__.py`) and covered by the same `wheel.packages` entry, and it must not add a
dependency — NumPy is used only for zero-copy *views* at the boundary, so it stays an optional
accessor (`.to_numpy()`), not an import at module scope.

### 9.2 The whole API a caller needs

```python
from seqtree import TextIndex

ix = TextIndex.build(list(seqs.values()), alphabet="protein", k=4,
                     group_ids=gene_id_per_record)          # group_ids optional
res = ix.search_batch(peptides, max_subs=2, exclude_exact=True,
                      best_only=True, group_by=True, threads=0)

for q, hits in zip(peptides, res):                          # res is iterable, hit-by-query
    for h in hits:
        h.ref_id, h.offset, h.n_subs, h.score
        h.mismatches                                        # [(pos, 'L', 'I'), ...]

res.groups[i]        # [(group_id, min_subs, n_hits), ...] when group_by
res.truncated[i]     # 1 if max_hits capped this query
res.to_numpy()       # dict of zero-copy arrays, for the 445k-query case
```

Convenience, so a caller needs no FASTA parser:

```python
ix = TextIndex.from_fasta("human.fasta.gz", k=4)   # names kept, gz handled
ix.save("human_k4.sti"); ix = TextIndex.load("human_k4.sti")   # mmap by default
```

### 9.3 The `mhcmatch` adapter

`Proteome.find_sources` becomes ~20 lines: build **one** `TextIndex` instead of one `Index` per
length, and map `(ref_id, offset)` → `(protein, position, window)` as it already does. Its return
shape (`{peptide: [SourceHit(name, pos, ref_peptide, n_subs, mutations)]}`) is unchanged — `mutations`
is now read straight from `mm_pos` / `mm_query_aa` / `mm_text_aa` instead of being recomputed in
Python — so `vector.py`, `cli.py` and the existing tests keep working.

What `mhcmatch` then deletes: `Proteome._index`, `_index_key`, `_index_from_disk`, `_index_to_disk`,
`_await_index`, `_index_building_path`, `index_cache_dir`, and the `O_EXCL` cross-process build
hand-off — an entire subsystem that exists only because the current index is expensive enough to be
worth racing over.

---

## 10. Out of scope

- **Gaps and scored alignment.** `TextIndex` answers a Hamming question over a text; `Index` and the
  seqtm/seqtrie engines keep the substitution-matrix and gap work. The optional matrix here *scores*
  a hit the Hamming predicate already accepted; it never changes which hits are returned.
- **`KmerIndex`.** Peptide-vs-peptide homology, keyed by peptide id, not position in a text. Stays.
- **Suffix arrays / FM-index.** They also serve any query length, but k-mismatch search on an
  FM-index is backtracking exponential in `m`, and the measured pigeonhole already meets the target.
  Revisit only if `max_subs` must exceed 3.
- **Edit distance.** Indels change the length and break the pigeonhole block argument entirely.
- **GPU.** PEPMatch names it as future work; §5.4 buys the same order of magnitude on CPU with none
  of the deployment cost.
