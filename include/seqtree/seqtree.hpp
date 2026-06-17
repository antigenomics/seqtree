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

// Non-negative substitution penalties indexed by alphabet code; penalty(a,a)==0.
class SubstitutionMatrix {
public:
    // Unit cost: 0 for a match, 1 for a mismatch.
    static SubstitutionMatrix unit(uint8_t size);
    // BLOSUM62 converted to penalties; valid only for the AminoAcid codec order.
    static SubstitutionMatrix blosum62();
    // Convert a similarity matrix (row-major, size*size) to penalties via
    // pen[a][b] = max(sim[a][a], sim[b][b]) - sim[a][b].
    static SubstitutionMatrix from_similarity(uint8_t size, const int32_t* sim);

    uint8_t size() const { return size_; }
    int32_t penalty(uint8_t a, uint8_t b) const { return pen_[a * size_ + b]; }

private:
    uint8_t size_ = 0;
    std::vector<int32_t> pen_;
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

    const Trie& trie() const { return *trie_; }

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

}  // namespace seqtree
