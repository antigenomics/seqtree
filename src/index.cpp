#include "seqtree/seqtree.hpp"
#include "trie.hpp"

#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace seqtree {

Index::Index() = default;
Index::~Index() = default;

namespace {
constexpr char kMagic[4] = {'S', 'Q', 'T', 'R'};
constexpr uint32_t kVersion = 1;

template <class T>
void write_vec(std::ostream& os, const std::vector<T>& v) {
    uint64_t n = v.size();
    os.write(reinterpret_cast<const char*>(&n), sizeof n);
    if (n) os.write(reinterpret_cast<const char*>(v.data()), std::streamsize(n * sizeof(T)));
}
template <class T>
void read_vec(std::istream& is, std::vector<T>& v) {
    uint64_t n = 0;
    is.read(reinterpret_cast<char*>(&n), sizeof n);
    v.resize(n);
    if (n) is.read(reinterpret_cast<char*>(v.data()), std::streamsize(n * sizeof(T)));
}
}  // namespace

void Index::save(const std::string& path) const {
    std::ofstream os(path, std::ios::binary);
    if (!os) throw std::runtime_error("seqtree: cannot open '" + path + "' for writing");
    const Trie& t = *trie_;
    os.write(kMagic, 4);
    os.write(reinterpret_cast<const char*>(&kVersion), sizeof kVersion);
    uint8_t alpha = static_cast<uint8_t>(t.codec.alphabet());
    os.write(reinterpret_cast<const char*>(&alpha), 1);
    os.write(reinterpret_cast<const char*>(&t.max_depth), sizeof(uint32_t));
    write_vec(os, t.nodes);
    write_vec(os, t.edge_code);
    write_vec(os, t.edge_child);
    write_vec(os, t.ref_ids);
    write_vec(os, t.str_off);
    uint64_t sd = t.str_data.size();
    os.write(reinterpret_cast<const char*>(&sd), sizeof sd);
    if (sd) os.write(t.str_data.data(), std::streamsize(sd));
    if (!os) throw std::runtime_error("seqtree: write failed for '" + path + "'");
}

std::unique_ptr<Index> Index::load(const std::string& path) {
    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("seqtree: cannot open '" + path + "'");
    char magic[4] = {0, 0, 0, 0};
    uint32_t version = 0;
    is.read(magic, 4);
    is.read(reinterpret_cast<char*>(&version), sizeof version);
    if (std::memcmp(magic, kMagic, 4) != 0 || version != kVersion)
        throw std::runtime_error("seqtree: '" + path + "' is not a v" +
                                 std::to_string(kVersion) + " seqtree index");
    uint8_t alpha = 0;
    is.read(reinterpret_cast<char*>(&alpha), 1);

    std::unique_ptr<Index> idx(new Index());
    idx->trie_ = std::make_unique<Trie>();
    Trie& t = *idx->trie_;
    t.codec = Codec(static_cast<Alphabet>(alpha));
    is.read(reinterpret_cast<char*>(&t.max_depth), sizeof(uint32_t));
    read_vec(is, t.nodes);
    read_vec(is, t.edge_code);
    read_vec(is, t.edge_child);
    read_vec(is, t.ref_ids);
    read_vec(is, t.str_off);
    uint64_t sd = 0;
    is.read(reinterpret_cast<char*>(&sd), sizeof sd);
    t.str_data.resize(sd);
    if (sd) is.read(&t.str_data[0], std::streamsize(sd));
    if (!is) throw std::runtime_error("seqtree: truncated or corrupt index '" + path + "'");
    return idx;
}

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
