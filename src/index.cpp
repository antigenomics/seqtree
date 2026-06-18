#include "seqtree/seqtree.hpp"
#include "trie.hpp"

#include <algorithm>
#include <atomic>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace seqtree {

Index::Index() = default;
Index::~Index() = default;

std::unique_ptr<Index> Index::build(std::vector<std::string> refs, Alphabet a) {
    std::unique_ptr<Index> idx(new Index());
    idx->trie_ = std::make_unique<Trie>(Trie::build(std::move(refs), a));
    return idx;
}

uint32_t Index::size() const { return trie_->size(); }
Alphabet Index::alphabet() const { return trie_->codec.alphabet(); }
const Codec& Index::codec() const { return trie_->codec; }
std::string_view Index::ref_seq(uint32_t ref_id) const { return trie_->ref_seq(ref_id); }

Alignment Index::align(std::string_view query, uint32_t ref_id, const SearchParams& p) const {
    std::string_view ref = ref_seq(ref_id);
    const Codec& cod = trie_->codec;
    const bool unit = (p.matrix == nullptr);
    const int32_t gap = p.gap_open;
    auto sub_pen = [&](char a, char b) -> int32_t {
        if (unit) return a == b ? 0 : 1;
        return p.matrix->penalty(cod.encode(a), cod.encode(b));
    };

    const int m = int(query.size()), n = int(ref.size());
    const int W = n + 1;
    std::vector<int32_t> dp(size_t(m + 1) * W);
    std::vector<uint8_t> bt(size_t(m + 1) * W);  // 0 diag, 1 up(del query), 2 left(ins ref)
    dp[0] = 0;
    for (int j = 1; j <= n; ++j) { dp[j] = j * gap; bt[j] = 2; }
    for (int i = 1; i <= m; ++i) { dp[size_t(i) * W] = i * gap; bt[size_t(i) * W] = 1; }
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            int32_t diag = dp[size_t(i - 1) * W + (j - 1)] + sub_pen(query[i - 1], ref[j - 1]);
            int32_t up = dp[size_t(i - 1) * W + j] + gap;       // query char unmatched -> deletion
            int32_t left = dp[size_t(i) * W + (j - 1)] + gap;   // ref char unmatched -> insertion
            int32_t best = diag; uint8_t b = 0;
            if (up < best) { best = up; b = 1; }
            if (left < best) { best = left; b = 2; }
            dp[size_t(i) * W + j] = best;
            bt[size_t(i) * W + j] = b;
        }
    }

    Alignment al;
    al.score = dp[size_t(m) * W + n];
    int i = m, j = n;
    while (i > 0 || j > 0) {
        uint8_t b = bt[size_t(i) * W + j];
        if (i > 0 && j > 0 && b == 0) {
            char q = query[i - 1], r = ref[j - 1];
            al.aligned_query += q; al.aligned_ref += r;
            al.ops += (q == r ? 'M' : 'S');
            --i; --j;
        } else if (i > 0 && b == 1) {
            al.aligned_query += query[i - 1]; al.aligned_ref += '-'; al.ops += 'D';
            --i;
        } else {
            al.aligned_query += '-'; al.aligned_ref += ref[j - 1]; al.ops += 'I';
            --j;
        }
    }
    std::reverse(al.aligned_query.begin(), al.aligned_query.end());
    std::reverse(al.aligned_ref.begin(), al.aligned_ref.end());
    std::reverse(al.ops.begin(), al.ops.end());
    return al;
}

std::vector<std::vector<Hit>> Index::search_batch(const std::vector<std::string>& queries,
                                                  const SearchParams& p, int threads) const {
    const size_t n = queries.size();
    std::vector<std::vector<Hit>> results(n);
    if (n == 0) return results;

    unsigned nt = threads > 0 ? unsigned(threads)
                              : std::max(1u, std::thread::hardware_concurrency());
    nt = std::min<unsigned>(nt, std::max<size_t>(1, n));

    std::atomic<size_t> next{0};
    // Adaptive chunk: large batches keep ~1024 for low atomic contention, but small
    // batches must split fine enough that every worker gets work (~8 chunks/thread).
    const size_t chunk = std::clamp<size_t>(n / (size_t(nt) * 8), size_t(1), size_t(1024));
    std::exception_ptr err;
    std::mutex emu;

    auto worker = [&] {
        Searcher s(*this);
        for (;;) {
            size_t start = next.fetch_add(chunk);
            if (start >= n) break;
            size_t end = std::min(n, start + chunk);
            for (size_t i = start; i < end; ++i) {
                try {
                    s.search_into(queries[i], p, results[i]);
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

std::vector<uint64_t> Index::collisions_batch(const std::vector<std::string>& queries,
                                              const SearchParams& p, int threads) const {
    const size_t n = queries.size();
    std::vector<uint64_t> out(n, 0);
    if (n == 0) return out;

    unsigned nt = threads > 0 ? unsigned(threads)
                              : std::max(1u, std::thread::hardware_concurrency());
    nt = std::min<unsigned>(nt, std::max<size_t>(1, n));

    std::atomic<size_t> next{0};
    const size_t chunk = std::clamp<size_t>(n / (size_t(nt) * 8), size_t(1), size_t(1024));
    std::exception_ptr err;
    std::mutex emu;

    auto worker = [&] {
        Searcher s(*this);
        std::vector<Hit> scratch;
        for (;;) {
            size_t start = next.fetch_add(chunk);
            if (start >= n) break;
            size_t end = std::min(n, start + chunk);
            for (size_t i = start; i < end; ++i) {
                try {
                    s.search_into(queries[i], p, scratch);
                    out[i] = s.last_collisions();
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
    return out;
}

std::vector<std::vector<Hit>> pairwise_batch(const std::vector<std::string>& a,
                                             const std::vector<std::string>& b,
                                             Alphabet alphabet, const SearchParams& p, int threads) {
    const bool index_a = a.size() >= b.size();  // index larger, stream smaller
    const std::vector<std::string>& refs = index_a ? a : b;
    const std::vector<std::string>& qry = index_a ? b : a;

    auto idx = Index::build(std::vector<std::string>(refs), alphabet);
    auto res = idx->search_batch(qry, p, threads);

    if (!index_a) return res;  // qry == a, ref_id already indexes b

    // ref == a, qry == b: transpose to a-major with ref_id pointing into b.
    std::vector<std::vector<Hit>> out(a.size());
    for (uint32_t bi = 0; bi < res.size(); ++bi)
        for (const Hit& h : res[bi])
            out[h.ref_id].push_back(Hit{bi, h.score, h.n_subs, h.n_ins, h.n_dels});
    return out;
}

}  // namespace seqtree
