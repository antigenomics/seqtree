#pragma once
#include "seqtree/types.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace seqtree {

struct Trie;     // internal, src/trie.hpp
struct Scratch;  // internal, src/engines.hpp

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
    // threads <= 0 => hardware_concurrency. Results aligned with input order.
    std::vector<std::vector<Hit>> search_batch(const std::vector<std::string>& queries,
                                               const SearchParams&, int threads = 0) const;

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
    void search_into(std::string_view query, const SearchParams&, std::vector<Hit>& out);
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
// threads <= 0 => hardware_concurrency.
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
// the diagonal zero. `matrix` may be null for unit cost. threads <= 0 => hardware_concurrency.
std::vector<int32_t> gapblock_matrix(const std::vector<std::string>& queries,
                                     const std::vector<std::string>& refs, Alphabet,
                                     const SubstitutionMatrix* matrix, int32_t gap_open,
                                     int32_t gap_extend, const std::vector<int32_t>& prior,
                                     uint32_t prior_width, int threads = 0);

}  // namespace seqtree
