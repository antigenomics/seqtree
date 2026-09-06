#pragma once
#include "seqtree/seqtree.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace seqtree {

struct TextStore;  // internal, src/text_index.cpp

// One row of the optional group fold: hits collapsed onto the build-time group ids.
struct TextGroup {
    uint32_t group_id = 0;
    uint16_t min_subs = 0;  // best (smallest) substitution count reaching this group
    uint32_t n_hits   = 0;  // how many hits fell in it, at any distance
};

struct TextQueryOpts {
    uint16_t max_subs      = 0;
    // Insertions + deletions allowed, capped independently of max_subs. 0 keeps the pure
    // Hamming predicate and its verification path untouched.
    //
    // The seed table does not change: with b = max_subs + max_indels + 1 disjoint blocks the
    // total edit count is below the block count, so by pigeonhole SOME block still matches
    // EXACTLY -- the existing exact lookup finds every occurrence, indels included. What that
    // costs is verification, which becomes a banded alignment over the 2*max_indels+1 possible
    // starts rather than one Hamming scan at a fixed start.
    //
    // Because every block must be exact, this needs L >= (max_subs + max_indels + 1) * k;
    // search_batch refuses a shorter query rather than answering incompletely.
    uint16_t max_indels    = 0;
    bool     exclude_exact = false;  // drop 0-mismatch hits (a candidate that IS the text)
    bool     best_only     = false;  // keep only the minimum-n_subs shell (still ALL of it)
    bool     group_by      = false;  // also fold hits onto the build-time group ids
    uint32_t max_hits      = 0;      // 0 = unlimited; a cap is reported, never silent
    // Optional. Scores a hit the Hamming predicate has ALREADY accepted; it never changes
    // which hits are returned. score = sum of similarity(query_aa, text_aa) over the
    // mismatched positions, so 0 is an exact match and less-negative is more conservative.
    const SubstitutionMatrix* matrix = nullptr;
};

// Results as flat parallel arrays, never a vector<vector<TextHit>> with a vector<Mismatch>
// inside each hit: a 445k-query run returns a handful of arrays rather than millions of
// objects, and every array is a zero-copy view at the Python boundary.
//
// Hits for query q are [query_begin[q], query_begin[q+1]); mismatches for hit h are
// [mm_begin[h], mm_begin[h+1]). Hits are ordered (n_subs, ref_id, offset) -- stable across
// runs and thread counts, so a downstream digest is stable.
struct TextResult {
    std::vector<uint32_t> query_begin;  // n_queries + 1, CSR into the hit arrays
    std::vector<uint32_t> ref_id, offset;
    std::vector<uint16_t> n_subs;
    // Indels in the alignment, and the matched substring's length in the text. With
    // max_indels = 0 these are 0, 0 and the query length for every hit. With indels the match
    // is no longer len(query) residues wide, so `length` is what bounds it:
    // ref_seq(ref_id)[offset : offset + length].
    std::vector<uint16_t> n_ins, n_dels, length;
    std::vector<int32_t>  score;
    std::vector<uint32_t> mm_begin;     // n_hits + 1, CSR into the mismatch arrays
    std::vector<uint16_t> mm_pos;
    std::vector<char>     mm_query_aa, mm_text_aa;
    std::vector<uint8_t>  truncated;    // n_queries; 1 if max_hits capped this query

    // Group fold; all empty unless TextQueryOpts::group_by.
    std::vector<uint32_t> group_begin;  // n_queries + 1, CSR into the group arrays
    std::vector<uint32_t> group_id, group_n_hits;
    std::vector<uint16_t> group_min_subs;

    size_t num_queries() const { return query_begin.empty() ? 0 : query_begin.size() - 1; }
    size_t num_hits() const { return ref_id.size(); }
};

// Exact k-mismatch (Hamming) search over a CONCATENATED reference text.
//
// Index builds a trie over reference *strings*, so asking a text question with it means
// enumerating every length-L window as its own string -- and again for every distinct query
// length. TextIndex keys on a k-mer seed table instead, so `k` is a property of the index and
// not of the query: ONE build serves every length and every max_subs.
//
// The predicate is Hamming only: full length, no gaps, no scoring matrix in the decision, no
// E-value. The answer is exact, not a heuristic -- every position within max_subs is returned.
// Immutable after build; share one TextIndex across threads.
class TextIndex {
public:
    // `refs` are whole records, NOT windows -- a proteome, a genome, a set of transcripts.
    // `k` is the seed width; the table is direct-addressed, so alphabet_size^k buckets are
    // allocated (24^4 = 331,776 for amino acids) and k is capped accordingly.
    // `group_ids` is optional and arbitrary -- gene ids, species, cluster labels. seqtree does
    // not know what they mean, only that hits can be folded onto them.
    // Residues outside the alphabet are kept as holes and counted in num_unknown(); refusing
    // the build over them would make the class unusable on a real proteome.
    static std::unique_ptr<TextIndex> build(const std::vector<std::string>& refs, Alphabet,
                                            uint8_t k = 4,
                                            const std::vector<uint32_t>& group_ids = {});
    ~TextIndex();

    uint32_t num_refs() const;
    uint64_t num_residues() const;  // excludes the inter-record sentinels
    // Text residues the codec does not name (U in a UniProt proteome, say). They are kept as
    // holes -- no seed spans one, no hit crosses one -- so they are reported rather than
    // silently dropped. A QUERY containing one is refused instead; see search_batch.
    uint64_t num_unknown() const;
    uint8_t  k() const;
    Alphabet alphabet() const;
    std::string ref_seq(uint32_t ref_id) const;

    // Every query must be at least k long: seeds are only indexed where a whole k-mer fits
    // inside a record, so a shorter query could match within k-L of a record end and be
    // missed. Raises rather than answering incompletely. threads <= 0 => hardware_concurrency.
    TextResult search_batch(const std::vector<std::string>& queries, const TextQueryOpts&,
                            int threads = 0) const;

    // Flat and mmap-able: the arrays sit at known, aligned offsets, so load(mmap=true) maps
    // the file instead of reading it and several processes share one copy of the pages.
    void save(const std::string& path) const;
    static std::unique_ptr<TextIndex> load(const std::string& path, bool mmap = true);

private:
    TextIndex();
    std::unique_ptr<TextStore> store_;
};

}  // namespace seqtree
