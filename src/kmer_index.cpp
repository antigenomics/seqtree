#include "seqtree/kmer_index.hpp"

#include "atomic_write.hpp"

#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <unordered_map>

namespace seqtree {

KmerIndex::KmerIndex() = default;
KmerIndex::~KmerIndex() = default;
uint32_t KmerIndex::num_kmers() const { return kmers_ ? kmers_->size() : 0; }

std::unique_ptr<KmerIndex> KmerIndex::build(
    const std::vector<std::vector<std::string>>& kmers_per_peptide, Alphabet alphabet,
    const std::vector<uint32_t>& allele_ids) {
    const uint32_t n = static_cast<uint32_t>(kmers_per_peptide.size());
    if (!allele_ids.empty() && allele_ids.size() != n)
        throw std::invalid_argument("KmerIndex::build: allele_ids size mismatch");

    // Assign a stable id to each unique k-mer.
    std::unordered_map<std::string, uint32_t> kmer_id;
    std::vector<std::string> unique_kmers;
    std::vector<std::vector<uint32_t>> postings;  // kmer_id -> peptide ids
    for (uint32_t p = 0; p < n; ++p) {
        // dedup this peptide's k-mers so a repeated k-mer posts once
        std::vector<uint32_t> ids;
        for (const auto& km : kmers_per_peptide[p]) {
            auto it = kmer_id.find(km);
            uint32_t id;
            if (it == kmer_id.end()) {
                id = static_cast<uint32_t>(unique_kmers.size());
                kmer_id.emplace(km, id);
                unique_kmers.push_back(km);
                postings.emplace_back();
            } else {
                id = it->second;
            }
            ids.push_back(id);
        }
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        for (uint32_t id : ids) postings[id].push_back(p);
    }

    auto ki = std::unique_ptr<KmerIndex>(new KmerIndex());
    ki->num_peptides_ = n;
    ki->allele_ = allele_ids.empty() ? std::vector<uint32_t>(n, 0) : allele_ids;
    ki->kmers_ = Index::build(unique_kmers, alphabet);  // ref_id == kmer_id (no dups)

    // flatten postings to CSR
    ki->post_begin_.assign(unique_kmers.size() + 1, 0);
    for (size_t i = 0; i < postings.size(); ++i)
        ki->post_begin_[i + 1] = ki->post_begin_[i] + uint32_t(postings[i].size());
    ki->post_ids_.reserve(ki->post_begin_.back());
    for (auto& lst : postings)
        ki->post_ids_.insert(ki->post_ids_.end(), lst.begin(), lst.end());
    return ki;
}

namespace {
struct Acc {
    uint32_t shared = 0;
    int32_t  best = 0;
    int32_t  last_q = -1;  // index of the query k-mer that last bumped `shared` (dedup)
};
}  // namespace

std::vector<std::vector<Candidate>> KmerIndex::seed_and_gather(
    const std::vector<std::vector<std::string>>& query_kmers, const SearchParams& params,
    uint32_t min_shared, int64_t allele_filter, int threads) const {
    const size_t nq = query_kmers.size();
    std::vector<std::vector<Candidate>> results(nq);
    if (nq == 0) return results;

    unsigned nt = threads > 0 ? unsigned(threads)
                              : std::max(1u, std::thread::hardware_concurrency());
    nt = std::min<unsigned>(nt, std::max<size_t>(1, nq));

    std::atomic<size_t> next{0};
    const size_t chunk = std::clamp<size_t>(nq / (size_t(nt) * 8), size_t(1), size_t(256));
    std::exception_ptr err;
    std::mutex emu;

    auto worker = [&] {
        Searcher s(*kmers_);
        std::unordered_map<uint32_t, Acc> acc;
        std::vector<Hit> hits;
        for (;;) {
            size_t start = next.fetch_add(chunk);
            if (start >= nq) break;
            size_t end = std::min(nq, start + chunk);
            for (size_t qi = start; qi < end; ++qi) {
                try {
                    acc.clear();
                    int qk = 0;
                    for (const auto& km : query_kmers[qi]) {
                        s.search_into(km, params, hits);
                        for (const Hit& h : hits) {  // h.ref_id == matched kmer_id
                            uint32_t kid = h.ref_id;
                            for (uint32_t pi = post_begin_[kid]; pi < post_begin_[kid + 1]; ++pi) {
                                uint32_t pid = post_ids_[pi];
                                if (allele_filter >= 0 && allele_[pid] != uint32_t(allele_filter))
                                    continue;
                                Acc& a = acc[pid];
                                if (a.last_q != qk) {  // first hit of THIS query k-mer for this peptide
                                    a.shared += 1;
                                    a.last_q = qk;
                                    if (a.shared == 1 || h.score < a.best) a.best = h.score;
                                } else if (h.score < a.best) {
                                    a.best = h.score;
                                }
                            }
                        }
                        ++qk;
                    }
                    auto& out = results[qi];
                    for (const auto& [pid, a] : acc)
                        if (a.shared >= min_shared)
                            out.push_back(Candidate{pid, a.shared, a.best});
                    std::sort(out.begin(), out.end(), [](const Candidate& x, const Candidate& y) {
                        if (x.shared_kmers != y.shared_kmers) return x.shared_kmers > y.shared_kmers;
                        if (x.best_score != y.best_score) return x.best_score < y.best_score;
                        return x.peptide_id < y.peptide_id;
                    });
                } catch (...) {
                    std::lock_guard<std::mutex> lk(emu);
                    if (!err) err = std::current_exception();
                    return;
                }
            }
        }
    };

    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    if (err) std::rethrow_exception(err);
    return results;
}

namespace {
constexpr char kMagic[4] = {'S', 'Q', 'K', 'M'};
constexpr uint32_t kVersion = 1;
template <class T>
void wvec(std::ostream& os, const std::vector<T>& v) {
    uint64_t n = v.size();
    os.write(reinterpret_cast<const char*>(&n), sizeof n);
    if (n) os.write(reinterpret_cast<const char*>(v.data()), std::streamsize(n * sizeof(T)));
}
template <class T>
void rvec(std::istream& is, std::vector<T>& v) {
    uint64_t n = 0;
    is.read(reinterpret_cast<char*>(&n), sizeof n);
    v.resize(n);
    if (n) is.read(reinterpret_cast<char*>(v.data()), std::streamsize(n * sizeof(T)));
}
}  // namespace

void KmerIndex::save(const std::string& path) const {
    // Temporary + rename, so a concurrent reader never sees a half-written index.
    detail::atomic_write(path, [this](std::ostream& os) {
        os.write(kMagic, 4);
        os.write(reinterpret_cast<const char*>(&kVersion), sizeof kVersion);
        uint8_t alpha = static_cast<uint8_t>(kmers_->alphabet());
        os.write(reinterpret_cast<const char*>(&alpha), 1);
        os.write(reinterpret_cast<const char*>(&num_peptides_), sizeof num_peptides_);
        // unique k-mer strings (rebuild the Index on load)
        uint32_t nk = kmers_->size();
        os.write(reinterpret_cast<const char*>(&nk), sizeof nk);
        std::vector<uint32_t> off(nk + 1, 0);
        std::string blob;
        for (uint32_t i = 0; i < nk; ++i) {
            blob += kmers_->ref_seq(i);
            off[i + 1] = uint32_t(blob.size());
        }
        wvec(os, off);
        uint64_t bn = blob.size();
        os.write(reinterpret_cast<const char*>(&bn), sizeof bn);
        if (bn) os.write(blob.data(), std::streamsize(bn));
        wvec(os, post_begin_);
        wvec(os, post_ids_);
        wvec(os, allele_);
    });
}

std::unique_ptr<KmerIndex> KmerIndex::load(const std::string& path) {
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("KmerIndex: cannot open '" + path + "'");
    char magic[4] = {0, 0, 0, 0};
    uint32_t version = 0;
    is.read(magic, 4);
    is.read(reinterpret_cast<char*>(&version), sizeof version);
    if (std::memcmp(magic, kMagic, 4) != 0 || version != kVersion)
        throw std::runtime_error("KmerIndex: '" + path + "' is not a v1 k-mer index");
    uint8_t alpha = 0;
    is.read(reinterpret_cast<char*>(&alpha), 1);
    auto ki = std::unique_ptr<KmerIndex>(new KmerIndex());
    is.read(reinterpret_cast<char*>(&ki->num_peptides_), sizeof ki->num_peptides_);
    uint32_t nk = 0;
    is.read(reinterpret_cast<char*>(&nk), sizeof nk);
    std::vector<uint32_t> off;
    rvec(is, off);
    uint64_t bn = 0;
    is.read(reinterpret_cast<char*>(&bn), sizeof bn);
    std::string blob(bn, '\0');
    if (bn) is.read(&blob[0], std::streamsize(bn));
    std::vector<std::string> kmers(nk);
    for (uint32_t i = 0; i < nk; ++i) kmers[i] = blob.substr(off[i], off[i + 1] - off[i]);
    ki->kmers_ = Index::build(std::move(kmers), static_cast<Alphabet>(alpha));
    rvec(is, ki->post_begin_);
    rvec(is, ki->post_ids_);
    rvec(is, ki->allele_);
    if (!is) throw std::runtime_error("KmerIndex: truncated/corrupt index '" + path + "'");
    return ki;
}

}  // namespace seqtree
