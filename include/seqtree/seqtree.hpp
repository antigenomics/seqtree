#pragma once
#include "seqtree/types.hpp"

#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

namespace seqtree {

struct Trie;     // internal, src/trie.hpp
struct Scratch;  // internal, src/engines.hpp
struct CountBuffer;

// Maps characters to small integer codes (<=32) and back. Case-insensitive on
// input. Invalid characters encode to kInvalid.
class Codec {
public:
    static constexpr uint8_t kInvalid = 0xFF;
    explicit Codec(Alphabet);
    uint8_t  size() const { return size_; }
    Alphabet alphabet() const { return alphabet_; }
    uint8_t  encode(char c) const { return enc_[static_cast<uint8_t>(c)]; }
    char     decode(uint8_t code) const { return dec_[code]; }

private:
    Alphabet alphabet_;
    uint8_t  size_ = 0;
    uint8_t  enc_[256];
    char     dec_[32];
};

// A substitution matrix in two views.
//
//   penalty(a,b)     non-negative, penalty(a,a)==0. What search and the E-value ball need:
//                    a distance-like cost that is minimised.
//   similarity(a,b)  the raw log-odds, signed. What Needleman-Wunsch and Smith-Waterman need:
//                    a score that is maximised.
//
// The penalty is the Gram / squared-distance transform of the similarity, which is LOSSY --
// it forces the diagonal to zero, destroying s(a,a) -- so the raw grid is retained rather
// than reconstructed.
class SubstitutionMatrix {
public:
    // Unit cost: 0 for a match, 1 for a mismatch (this is the "identity" matrix).
    static SubstitutionMatrix unit(uint8_t size);
    // Built-in amino-acid matrices (valid only for the AminoAcid codec order).
    static SubstitutionMatrix blosum62();
    static SubstitutionMatrix blosum45();  // shallower: remote homologs
    static SubstitutionMatrix blosum80();  // deeper: close homologs
    static SubstitutionMatrix pam250();   // EMBOSS EPAM250 (NCBI log-odds)
    static SubstitutionMatrix pam100();   // EMBOSS EPAM100 (NCBI log-odds)
    static SubstitutionMatrix structural();  // Miyazawa-Jernigan interaction-strength similarity
    // Convert a similarity matrix (row-major, size*size) to penalties via the Gram /
    // squared-distance transform pen[a][b] = sim[a][a] + sim[b][b] - 2*sim[a][b]. The raw
    // similarity is kept too.
    static SubstitutionMatrix from_similarity(uint8_t size, const int32_t* sim);

    uint8_t size() const { return size_; }
    int32_t penalty(uint8_t a, uint8_t b) const { return pen_[a * size_ + b]; }
    int32_t similarity(uint8_t a, uint8_t b) const { return sim_[a * size_ + b]; }
    // Median penalty over all mismatched symbol pairs -- the natural unit of this matrix.
    // Gap costs must be on this scale: the Gram transform makes a typical BLOSUM62 mismatch
    // cost ~15, so the default gap_open of 1 would make gaps ~15x cheaper than substitutions
    // and the aligner would gap rather than substitute. Use gap_open ~ 1-2 * scale().
    //
    // This is the PENALTY scale. A similarity-scoring aligner lives on the raw log-odds scale
    // instead, where the usual BLAST defaults (gap_open 11, gap_extend 1) apply.
    int32_t scale() const;

private:
    uint8_t size_ = 0;
    std::vector<int32_t> pen_;
    std::vector<int32_t> sim_;
};

// Per-position substitution penalties pen(pos, a, b) over a fixed frame width W.
// Two construction modes: a base SubstitutionMatrix scaled by per-position integer
// weights (weight 0 == masked/free, e.g. anchors; >1 == up-weighted, e.g. a TCR
// hotspot), or a full per-position PSSM table. masked(pos) reports weight-0 columns
// so they neither count as substitutions nor add penalty.
class PositionalMatrix {
public:
    // pen[pos][a][b] = weights[pos] * base.penalty(a, b); weights has length width.
    static PositionalMatrix from_weights(const SubstitutionMatrix& base,
                                         const std::vector<int32_t>& weights);
    // Full table: data is row-major [width][size][size]; masked[] optional (len width).
    static PositionalMatrix from_tables(uint8_t size, uint16_t width,
                                        const std::vector<int32_t>& data,
                                        const std::vector<uint8_t>& masked = {});

    uint8_t  size()  const { return size_; }
    uint16_t width() const { return width_; }
    bool     masked(uint16_t pos) const { return masked_[pos] != 0; }
    int32_t  penalty(uint16_t pos, uint8_t a, uint8_t b) const {
        return pen_[(size_t(pos) * size_ + a) * size_ + b];
    }

private:
    uint8_t  size_  = 0;
    uint16_t width_ = 0;
    std::vector<int32_t> pen_;     // [width * size * size]
    std::vector<uint8_t> masked_;  // [width]
};

struct AlignmentMaxResult {
    std::vector<double> scores, raw_similarity, normalized_similarity, pfit, informativeness;
    std::vector<uint64_t> counts, modeled_counts, aa_counts, exact_counts, candidate_counts, nw_pairs;
    std::vector<int64_t> best_ref_ids, best_lane_ids;
    std::vector<uint32_t> gap_counts;
    std::vector<uint8_t> query_supported;
};

// Immutable after build. Lock-free for concurrent reads; share one Index across
// threads, give each thread its own Searcher.
class Index {
public:
    static std::unique_ptr<Index> build(std::vector<std::string> refs, Alphabet);
    ~Index();

    uint32_t size() const;
    Alphabet alphabet() const;
    const Codec& codec() const;
    std::string_view ref_seq(uint32_t ref_id) const;

    // On-demand global alignment between query and a specific reference. Never
    // computed during search.
    Alignment align(std::string_view query, uint32_t ref_id, const SearchParams&) const;

    // Parallel batch: pure C++ (no Python), one Searcher per worker thread.
    // threads <= 0 => available CPU budget. Results aligned with input order.
    std::vector<std::vector<Hit>> search_batch(const std::vector<std::string>& queries,
                                               const SearchParams&, int threads = 0) const;

    // One pool, serialized owned primitive blocks; each query is emitted completely.
    // Columns: query index, ref_id, score, n_subs, n_ins, n_dels (all nonnegative).
    void search_batch_consume(const std::vector<std::string>& queries, const SearchParams&,
        const std::function<void(size_t, size_t, CountBuffer&&)>& consume,
        int threads = 0, size_t chunk_size = 32, size_t max_hit_bytes = 64*1024*1024) const;

    // Coordinate two indexes in one pool; both owned buffers cover the same queries.
    void search_paired_batch_consume(const Index& other,
        const std::vector<std::string>& queries, const std::vector<std::string>& other_queries,
        const SearchParams&,
        const std::function<void(size_t, size_t, CountBuffer&&, CountBuffer&&)>& consume,
        int threads = 0, size_t chunk_size = 32, size_t max_hit_bytes = 64*1024*1024, bool skip_empty = false) const;

    // Exact edit-count bins 0..max_total_edits; one reusable hit buffer per worker.
    std::vector<std::vector<uint64_t>> edit_histogram_batch(
        const std::vector<std::string>& queries, const SearchParams&,
        int threads = 0, bool exclude_exact = false, int max_indels = -1) const;

    // Numeric requested-label maxima; float NW selection is independent of candidate scope.
    AlignmentMaxResult alignment_max_batch(const std::vector<std::string>& queries,
        const SearchParams&, const std::vector<std::vector<double>>& similarity, const std::string& symbols, int max_indels,
        const std::vector<uint32_t>& query_target_offsets, const std::vector<uint32_t>& query_target_ids,
        const std::vector<uint32_t>& reference_target_offsets, const std::vector<uint32_t>& reference_target_ids,
        const std::vector<int32_t>& target_class_ids, const std::vector<uint32_t>& reference_class_offsets,
        const std::vector<int32_t>& reference_class_ids, const std::vector<std::vector<int32_t>>& query_group_ids,
        const std::vector<std::vector<int32_t>>& reference_group_ids, const std::vector<std::vector<std::vector<double>>>& group_offsets,
        const std::vector<uint64_t>& reference_other_counts, const std::vector<uint64_t>& reference_multiplicities,
        const std::vector<double>& coefficients, int threads=0, bool exclude_exact=false) const;

    // Per-query count of seqtm collisions: how often a reference was re-reached via a
    // different edit path during branch-and-bound (0 for the seqtrie engine / Hamming).
    std::vector<uint64_t> collisions_batch(const std::vector<std::string>& queries,
                                           const SearchParams&, int threads = 0) const;

    const Trie& trie() const { return *trie_; }

    // Serialize the frozen index to a flat binary file (little-endian) and load it
    // back. load() throws std::runtime_error on a missing/corrupt/version-mismatched file.
    void save(const std::string& path) const;
    static std::unique_ptr<Index> load(const std::string& path);

private:
    Index();
    std::unique_ptr<Trie> trie_;
};

// Per-thread search context holding reusable scratch buffers.
class Searcher {
public:
    explicit Searcher(const Index& idx);
    ~Searcher();

    std::vector<Hit> search(std::string_view query, const SearchParams&);
    void search_into(std::string_view query, const SearchParams&, std::vector<Hit>& out, int max_indels = -1);
    bool search_top(std::string_view query, const SearchParams&, Hit& out);

    // seqtm collisions from the most recent search_into on this Searcher.
    uint64_t last_collisions() const;

private:
    const Index& idx_;
    std::unique_ptr<Scratch> scratch_;
    std::vector<uint8_t> qcodes_;
};

// Batch-vs-batch. Indexes the larger set and streams the smaller, then returns
// results in a-major order (each inner list is for a[i]; Hit.ref_id indexes b).
std::vector<std::vector<Hit>> pairwise_batch(const std::vector<std::string>& a,
                                             const std::vector<std::string>& b,
                                             Alphabet, const SearchParams&, int threads = 0);

// ---------------------------------------------------------------------------------------
// Pairwise similarity alignment: Needleman-Wunsch (global) and Smith-Waterman (local).
//
// These MAXIMISE the raw log-odds similarity, unlike everything else in this header, which
// minimises a non-negative penalty. They exist so a caller does not need BioPython for an
// ordinary protein alignment; the conventions match Bio.Align.PairwiseAligner exactly:
//
//   * a gap of length L costs  gap_open + (L-1)*gap_extend   (so gap_open is the cost of the
//     first gap column, not an extra charge on top of it);
//   * gap_open == gap_extend gives LINEAR gaps -- no separate mode is needed;
//   * in global mode END gaps are charged like any other (true Needleman-Wunsch, not
//     semi-global / overlap);
//   * in local mode the score never drops below zero and the best local cell wins
//     (Smith-Waterman).
//
// Gap costs are given as POSITIVE magnitudes and subtracted. BLAST's protein defaults are
// gap_open = 11, gap_extend = 1.
enum class AlignMode { Global, Local };

// Optimal score only: O(min(m,n)) memory, no traceback. This is the hot path.
int32_t align_score(std::string_view query, std::string_view ref, const SubstitutionMatrix&,
                    Alphabet, AlignMode, int32_t gap_open, int32_t gap_extend);

// Optimal score plus the aligned strings and ops. O(m*n) memory for the traceback.
Alignment align_pair(std::string_view query, std::string_view ref, const SubstitutionMatrix&,
                     Alphabet, AlignMode, int32_t gap_open, int32_t gap_extend);

// Dense N*K similarity matrix, row-major (row i is queries[i] vs every ref). Parallel;
// threads <= 0 => available CPU budget.
std::vector<int32_t> align_score_matrix(const std::vector<std::string>& queries,
                                        const std::vector<std::string>& refs,
                                        const SubstitutionMatrix&, Alphabet, AlignMode,
                                        int32_t gap_open, int32_t gap_extend, int threads = 0);

// Dense N*K distance matrix d(a,b) = s(a,a) + s(b,b) - 2*s(a,b) >= 0, the Gram transform
// applied at the SEQUENCE level to the alignment scores above. This is the distance a
// prototype-embedding wants, and the one BioPython users hand-roll.
std::vector<int32_t> align_dist_matrix(const std::vector<std::string>& queries,
                                       const std::vector<std::string>& refs,
                                       const SubstitutionMatrix&, Alphabet, AlignMode,
                                       int32_t gap_open, int32_t gap_extend, int threads = 0);

// Dense N*K single-gap-block penalty matrix, row-major (row i is queries[i] vs every ref).
// Unlike pairwise_batch this is exhaustive: no budget, no trie, every cell scored.
//
// The block of length d = |m-n| is placed at whichever of the min(m,n)+1 positions minimises
// substitution cost plus `prior`. `prior` is the gap prior flattened to [m][d][i] with stride
// (prior_width+1); pass it empty for no prior. It is read only when d > 0, which is what keeps
// the diagonal zero. `matrix` may be null for unit cost. threads <= 0 => available CPU budget.
// Optional integer weights are indexed [longer length][longer position], symmetric in q/r.
// All input frames must be covered; empty tables retain the existing unweighted scorer.
// Alternatively, per-input rows index each sequence's full coordinates. Aligned
// substitution costs use max(query weight, reference weight); gaps are unchanged.
// Both axes are required; these rows cannot combine with length profiles.
// Optional group IDs index a nonnegative rectangular penalty table added per pair.
// Exact exclusion remains encoded junction identity, independent of these groups.
// Factored nonnegative offsets: paired feature codes, one shared code cost table.
// Round weight * sum(max(head cost, tail cost)) / depth to nearest, ties even,
// then multiply by scale. Input descriptors are copied into immutable batch state.
struct PairedFeatureOffsets {
    std::vector<std::vector<int32_t>> query_codes, reference_codes, costs;
    double weight;
    int32_t scale;
    PairedFeatureOffsets(std::vector<std::vector<int32_t>> q,
                         std::vector<std::vector<int32_t>> r,
                         std::vector<std::vector<int32_t>> c, double w = 1.0, int32_t s = 1)
        : query_codes(std::move(q)), reference_codes(std::move(r)), costs(std::move(c)),
          weight(w), scale(s) {}
};

std::vector<int32_t> gapblock_matrix(const std::vector<std::string>& queries,
                                     const std::vector<std::string>& refs, Alphabet,
                                     const SubstitutionMatrix* matrix, int32_t gap_open,
                                     int32_t gap_extend, const std::vector<int32_t>& prior,
                                     uint32_t prior_width, int threads = 0,
                                     const std::vector<std::vector<int32_t>>& weights = {},
                                     const std::vector<int32_t>& query_groups = {},
                                     const std::vector<int32_t>& reference_groups = {},
                                     const std::vector<std::vector<int32_t>>& group_distances = {},
        const std::vector<std::vector<int32_t>>& query_weights = {},
        const std::vector<std::vector<int32_t>>& reference_weights = {},
        const std::optional<PairedFeatureOffsets>& feature_offsets = std::nullopt);

// Exhaustive ball rows in reference order; hit edit fields are unset (zero).
std::vector<std::vector<Hit>> gapblock_ball_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, const std::vector<int32_t>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights = {},
        const std::vector<int32_t>& query_groups = {},
        const std::vector<int32_t>& reference_groups = {},
        const std::vector<std::vector<int32_t>>& group_distances = {},
        const std::vector<std::vector<int32_t>>& query_weights = {},
        const std::vector<std::vector<int32_t>>& reference_weights = {},
        const std::optional<PairedFeatureOffsets>& feature_offsets = std::nullopt);

// Exhaustive bounded reductions using the same gapblock cell. Top-k ties use
// (penalty, ref_id); paired penalty is max(lane scores). Exact exclusion uses
// encoded sequence identity (both lanes for pairs), never zero penalty alone.
// Count the component ball before/after group offsets; minimum full score within it.
// One threshold per query; nearest=-1 denotes an empty component ball. O(Q) output.
std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<int32_t>>
gapblock_group_audit_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, const std::vector<int32_t>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights = {},
        const std::vector<int32_t>& query_groups = {},
        const std::vector<int32_t>& reference_groups = {},
        const std::vector<std::vector<int32_t>>& group_distances = {},
        const std::vector<std::vector<int32_t>>& query_weights = {},
        const std::vector<std::vector<int32_t>>& reference_weights = {},
        const std::optional<PairedFeatureOffsets>& feature_offsets = std::nullopt,
        const std::vector<int32_t>& query_partitions = {},
        const std::vector<int32_t>& reference_partitions = {});

std::vector<std::vector<Hit>> gapblock_topk_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, uint32_t k, Alphabet alphabet,
        const SubstitutionMatrix* matrix, int32_t go, int32_t ge, const std::vector<int32_t>& prior,
        uint32_t width, int threads, bool exclude_exact);
std::vector<std::vector<std::tuple<uint32_t, int32_t, int32_t>>> gapblock_paired_topk_batch(
        const std::vector<std::string>& qa, const std::vector<std::string>& qb,
        const std::vector<std::string>& ra, const std::vector<std::string>& rb,
        uint32_t k, Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact);
std::vector<std::vector<uint64_t>> gapblock_count_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, const std::vector<std::vector<int32_t>>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights = {},
        const std::vector<int32_t>& query_groups = {},
        const std::vector<int32_t>& reference_groups = {},
        const std::vector<std::vector<int32_t>>& group_distances = {},
        const std::vector<std::vector<int32_t>>& query_weights = {},
        const std::vector<std::vector<int32_t>>& reference_weights = {},
        const std::optional<PairedFeatureOffsets>& feature_offsets = std::nullopt, bool shared_thresholds = false);
std::pair<std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>> gapblock_count_mass_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, const std::vector<std::vector<int32_t>>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights = {},
        const std::vector<int32_t>& query_groups = {},
        const std::vector<int32_t>& reference_groups = {},
        const std::vector<std::vector<int32_t>>& group_distances = {},
        const std::vector<std::vector<int32_t>>& query_weights = {},
        const std::vector<std::vector<int32_t>>& reference_weights = {},
        const std::optional<PairedFeatureOffsets>& feature_offsets = std::nullopt, bool shared_thresholds = false);
std::vector<std::vector<uint64_t>> gapblock_paired_count_batch(
        const std::vector<std::string>& qa, const std::vector<std::string>& qb,
        const std::vector<std::string>& ra, const std::vector<std::string>& rb,
        const std::vector<std::vector<int32_t>>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact);
std::pair<std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>> gapblock_paired_sum_count_batch(
        const std::vector<std::string>& qa, const std::vector<std::string>& qb,
        const std::vector<std::string>& ra, const std::vector<std::string>& rb,
        const std::vector<std::vector<int32_t>>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights_alpha,
        const std::vector<std::vector<int32_t>>& weights_beta,
        const std::vector<int32_t>& query_groups_alpha, const std::vector<int32_t>& query_groups_beta,
        const std::vector<int32_t>& reference_groups_alpha, const std::vector<int32_t>& reference_groups_beta,
        const std::vector<std::vector<int32_t>>& group_distances_alpha,
        const std::vector<std::vector<int32_t>>& group_distances_beta,
        const std::vector<std::vector<int32_t>>& query_weights_alpha = {},
        const std::vector<std::vector<int32_t>>& query_weights_beta = {},
        const std::vector<std::vector<int32_t>>& reference_weights_alpha = {},
        const std::vector<std::vector<int32_t>>& reference_weights_beta = {},
        const std::optional<PairedFeatureOffsets>& feature_offsets_alpha = std::nullopt,
        const std::optional<PairedFeatureOffsets>& feature_offsets_beta = std::nullopt, bool shared_thresholds = false);

// ---------------------------------------------------------------------------------------
// Plain string edit distances on raw characters, unit costs. Deliberately alphabet-agnostic:
// no Codec, no substitution matrix, no gap model -- the distance you want when you just need a
// number, not a scored alignment. Comparison is byte-for-byte and CASE-SENSITIVE; unlike the
// search Codec, these do not fold case.
//
//   hamming     -- number of differing positions; defined only for EQUAL-length sequences and
//                  throws std::invalid_argument otherwise;
//   levenshtein -- classic insertion/deletion/substitution distance, each edit cost 1,
//                  O(min(m,n)) memory.
int hamming(std::string_view a, std::string_view b);
int levenshtein(std::string_view a, std::string_view b);

// Dense N*K distance matrices, row-major (row i is a[i] vs every b[k]). Parallel;
// threads <= 0 => available CPU budget. hamming_matrix throws on any length-mismatched pair.
std::vector<int32_t> hamming_matrix(const std::vector<std::string>& a,
                                    const std::vector<std::string>& b, int threads = 0);
std::vector<int32_t> levenshtein_matrix(const std::vector<std::string>& a,
                                        const std::vector<std::string>& b, int threads = 0);

// Generic label histograms over one prepared native query/reference population.
struct GroupedGeometry {
    std::vector<PairedFeatureOffsets> additional_features;
    std::optional<std::vector<int32_t>> gap_prior;
    std::string alignment_mode = "gapblock";
    int32_t alignment_scale = 1;
    uint32_t trim_left = 0, trim_right = 0;
    int32_t score_quantum = 1;
    bool return_hits = false;
    size_t max_hit_bytes = 64*1024*1024;
    std::optional<std::vector<int32_t>> query_target_ids;
    std::optional<std::vector<int64_t>> query_reference_ids;
    std::vector<std::vector<int32_t>> weights, group_distances, query_weights, reference_weights;
    std::vector<int32_t> query_groups, reference_groups;
    std::optional<PairedFeatureOffsets> features;
};
struct CountBuffer {
    std::vector<uint32_t> data;
    std::vector<size_t> shape;
};
struct GroupedCountChunk {
    size_t start = 0, stop = 0;
    CountBuffer counts, sizes, exact_counts;
    CountBuffer hits{{},{0,3}};
};
class GroupedCountBatch {
public:
    GroupedCountBatch(const std::vector<std::string>& qa, const std::vector<std::string>& ra,
        std::vector<int32_t> thresholds, std::vector<uint32_t> target_offsets,
        std::vector<uint32_t> target_ids, std::vector<int32_t> target_peptides,
        std::vector<int32_t> exclusive_peptides, std::vector<int32_t> query_genes,
        std::vector<int32_t> reference_genes, size_t chunk_size, Alphabet alphabet,
        const SubstitutionMatrix* matrix, int32_t go, int32_t ge, std::vector<int32_t> prior,
        uint32_t width, int threads, bool exclude_exact, const GroupedGeometry& alpha,
        const std::optional<GroupedGeometry>& beta = std::nullopt,
        const std::vector<std::string>& qb = {}, const std::vector<std::string>& rb = {},
        std::optional<std::vector<int32_t>> query_withheld_peptides = std::nullopt,
        std::vector<uint32_t> reference_peptide_offsets = {},
        std::vector<int32_t> reference_peptide_ids = {});
    std::optional<GroupedCountChunk> next();
    size_t query_count() const;
    size_t reference_count() const;
    size_t target_count() const;
    size_t peptide_count() const;
private:
    struct Impl;
    std::shared_ptr<Impl> impl_;
};

}  // namespace seqtree
