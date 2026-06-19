#pragma once
#include "seqtree/seqtree.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace seqtree {

// One scored candidate peptide from a seed-and-gather query.
struct Candidate {
    uint32_t peptide_id   = 0;
    uint32_t shared_kmers = 0;  // number of distinct query k-mers that hit this peptide
    int32_t  best_score   = 0;  // best (min) per-k-mer penalty among the hits
};

// Seed-and-extend k-mer index for homology search over (variable-length) peptides.
// Each unique k-mer maps to the peptides containing it (CSR posting lists) with a
// per-peptide allele tag. seed_and_gather fuzzy-matches query k-mers with seqtm and
// merges the posting lists into ranked candidates entirely in C++ (GIL released by the
// binding), so candidate generation for large mimic scans never crosses into Python.
// Anchor-aware k-mer extraction is done by the caller (the Python pMHC layer), which
// owns the parametrized AnchorSpec; this class is the hot query-time path.
class KmerIndex {
public:
    // kmers_per_peptide[i] = the (already canonicalized / anchor-masked) k-mers of peptide i;
    // allele_ids[i] = its allele tag (0 if unused). k-mers may be of mixed length; a query
    // k-mer only matches reference k-mers of equal length under the Hamming seed search.
    static std::unique_ptr<KmerIndex> build(
        const std::vector<std::vector<std::string>>& kmers_per_peptide,
        Alphabet alphabet, const std::vector<uint32_t>& allele_ids);
    ~KmerIndex();

    uint32_t num_peptides() const { return num_peptides_; }
    uint32_t num_kmers() const;
    Alphabet alphabet() const { return kmers_->alphabet(); }

    // For each query (its k-mer list), return candidates with >= min_shared shared k-mers,
    // ranked by (shared_kmers desc, best_score asc). `params` drives the per-k-mer fuzzy
    // match (max_subs / pos_matrix / engine). allele_filter >= 0 keeps only peptides with
    // that allele tag; < 0 means no filter. threads <= 0 => hardware_concurrency.
    std::vector<std::vector<Candidate>> seed_and_gather(
        const std::vector<std::vector<std::string>>& query_kmers, const SearchParams& params,
        uint32_t min_shared = 1, int64_t allele_filter = -1, int threads = 0) const;

    void save(const std::string& path) const;
    static std::unique_ptr<KmerIndex> load(const std::string& path);

private:
    KmerIndex();
    std::unique_ptr<Index> kmers_;      // unique k-mers; ref_id == kmer_id
    std::vector<uint32_t>  post_begin_; // size num_kmers + 1
    std::vector<uint32_t>  post_ids_;   // peptide ids grouped per k-mer (CSR)
    std::vector<uint32_t>  allele_;     // size num_peptides
    uint32_t num_peptides_ = 0;
};

}  // namespace seqtree
