#include "seqtree/seqtree.hpp"
#include "seqtree/parallel.hpp"
#include "atomic_write.hpp"
#include "trie.hpp"
#include "count_unit.hpp"

#include <algorithm>
#include <atomic>
#include <cstring>
#include <fstream>
#include <iterator>
#include <limits>
#include <mutex>
#include <stdexcept>
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
    // Written to a temporary and renamed into place, so a concurrent reader never sees a
    // half-written index -- see atomic_write.hpp.
    detail::atomic_write(path, [this](std::ostream& os) {
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
    });
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

// Global affine-gap alignment (Gotoh). A gap run of length L costs
// gap_open + (L-1)*gap_extend; with gap_open == gap_extend this is the linear-gap NW it
// replaces. States: M (match/sub), X (gap in ref -> 'D'), Y (gap in query -> 'I').
Alignment Index::align(std::string_view query, uint32_t ref_id, const SearchParams& p) const {
    std::string_view ref = ref_seq(ref_id);
    const Codec& cod = trie_->codec;
    const bool unit = (p.matrix == nullptr);
    for (char c : query)  // the ref comes from the index and is already validated
        if (cod.encode(c) == Codec::kInvalid)
            throw std::invalid_argument(std::string("invalid symbol '") + c + "' in query");
    if (p.gap_open < 0 || p.gap_extend < 0)
        throw std::invalid_argument("seqtree: gap_open and gap_extend must be >= 0");

    const int32_t go = p.gap_open, ge = p.gap_extend;
    // Compare through the codec, never raw chars. The codec is case-insensitive, so 'a' and 'A'
    // are the SAME residue -- and search() has always treated them that way. align() did not:
    // unit cost called them a mismatch (penalty 1) while the matrix path called them identical
    // (penalty 0), so the two modes disagreed, and the ops string labelled them 'S'.
    auto same = [&](char a, char b) { return cod.encode(a) == cod.encode(b); };
    auto sub_pen = [&](char a, char b) -> int32_t {
        if (unit) return same(a, b) ? 0 : 1;
        return p.matrix->penalty(cod.encode(a), cod.encode(b));
    };

    constexpr int32_t kInf = std::numeric_limits<int32_t>::max() / 4;
    const int m = int(query.size()), n = int(ref.size());
    const int W = n + 1;
    const size_t sz = size_t(m + 1) * W;
    std::vector<int32_t> M(sz, kInf), X(sz, kInf), Y(sz, kInf);
    std::vector<uint8_t> bM(sz, 0), bX(sz, 0), bY(sz, 0);  // predecessor state: 0=M 1=X 2=Y

    M[0] = 0;
    for (int i = 1; i <= m; ++i) { X[size_t(i) * W] = go + (i - 1) * ge; bX[size_t(i) * W] = (i > 1) ? 1 : 0; }
    for (int j = 1; j <= n; ++j) { Y[j] = go + (j - 1) * ge; bY[j] = (j > 1) ? 2 : 0; }

    // best of (M, X, Y) candidates, returning the winning state in `which` (0=M 1=X 2=Y)
    auto best3 = [](int32_t mv, int32_t xv, int32_t yv, uint8_t& which) {
        int32_t v = mv; which = 0;
        if (xv < v) { v = xv; which = 1; }
        if (yv < v) { v = yv; which = 2; }
        return v;
    };
    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            const size_t k = size_t(i) * W + j, kd = size_t(i - 1) * W + (j - 1);
            const size_t ku = size_t(i - 1) * W + j, kl = size_t(i) * W + (j - 1);
            M[k] = sub_pen(query[i - 1], ref[j - 1]) + best3(M[kd], X[kd], Y[kd], bM[k]);
            X[k] = best3(M[ku] + go, X[ku] + ge, Y[ku] + go, bX[k]);  // consume query char
            Y[k] = best3(M[kl] + go, X[kl] + go, Y[kl] + ge, bY[k]);  // consume ref char
        }
    }

    Alignment al;
    const size_t end = size_t(m) * W + n;
    uint8_t st;
    al.score = best3(M[end], X[end], Y[end], st);

    int i = m, j = n;
    while (i > 0 || j > 0) {
        const size_t k = size_t(i) * W + j;
        if (st == 0) {  // M: consumed one query and one ref char
            char q = query[i - 1], r = ref[j - 1];
            al.aligned_query += q; al.aligned_ref += r;
            al.ops += (same(q, r) ? 'M' : 'S');
            st = bM[k]; --i; --j;
        } else if (st == 1) {  // X: query char unmatched -> deletion
            al.aligned_query += query[i - 1]; al.aligned_ref += '-'; al.ops += 'D';
            st = bX[k]; --i;
        } else {  // Y: ref char unmatched -> insertion
            al.aligned_query += '-'; al.aligned_ref += ref[j - 1]; al.ops += 'I';
            st = bY[k]; --j;
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

    parallel_for(n, threads, [&] { return Searcher(*this); },
                 [&](size_t i, Searcher& s) { s.search_into(queries[i], p, results[i]); });
    return results;
}

void Index::search_batch_consume(const std::vector<std::string>& queries,
        const SearchParams& p,
        const std::function<void(size_t, size_t, CountBuffer&&)>& consume,
        int threads, size_t chunk_size, size_t max_hit_bytes) const {
    constexpr size_t columns = 6;
    const size_t max_rows = max_hit_bytes / (columns * sizeof(uint32_t));
    if (!chunk_size || !max_rows || !consume)
        throw std::invalid_argument("positive chunk_size/max_hit_bytes and a consumer required");
    if (queries.size() > uint64_t(std::numeric_limits<uint32_t>::max()) + 1)
        throw std::overflow_error("query indices exceed the uint32 output buffer");
    // Bound both the primitive output and the per-worker Hit vectors. A single
    // query may still exceed the explicit output cap; detect that without truncation.
    chunk_size = std::min(chunk_size, std::max(size_t(1), max_rows / std::max(size_t(1), size_t(size()))));
    const size_t blocks = queries.size() / chunk_size + (queries.size() % chunk_size != 0);
    std::atomic<bool> cancelled{false};
    std::mutex callback_mutex;
    struct Local { Searcher searcher; std::vector<std::vector<Hit>> hits; };
    parallel_for(blocks, threads,
        [&] { return Local{Searcher(*this), std::vector<std::vector<Hit>>(chunk_size)}; },
        [&](size_t block, Local& worker) {
            if (cancelled.load()) return;
            try {
                const size_t start = block * chunk_size;
                const size_t stop = start + std::min(chunk_size, queries.size() - start);
                size_t rows = 0;
                for (size_t i = start; i < stop; ++i) {
                    if (cancelled.load()) return;
                    auto& hits = worker.hits[i-start];
                    worker.searcher.search_into(queries[i], p, hits);
                    if (hits.size() > max_rows - rows)
                        throw std::overflow_error("complete query hits exceed max_hit_bytes; increase the explicit buffer limit");
                    rows += hits.size();
                }
                CountBuffer buffer;
                buffer.shape = {rows, columns};
                buffer.data.resize(rows * columns);
                size_t offset = 0;
                for (size_t i = start; i < stop; ++i)
                    for (const auto& hit : worker.hits[i-start]) {
                        const uint32_t values[] = {uint32_t(i), hit.ref_id, uint32_t(hit.score),
                            hit.n_subs, hit.n_ins, hit.n_dels};
                        std::copy(std::begin(values), std::end(values), buffer.data.begin() + offset);
                        offset += columns;
                    }
                std::lock_guard<std::mutex> lock(callback_mutex);
                if (!cancelled.load()) {
                    try { consume(start, stop, std::move(buffer)); }
                    catch (...) { cancelled.store(true); throw; }
                }
            } catch (...) {
                cancelled.store(true);
                throw;
            }
        });
}

void Index::search_paired_batch_consume(const Index& other,
        const std::vector<std::string>& queries, const std::vector<std::string>& other_queries,
        const SearchParams& p,
        const std::function<void(size_t, size_t, CountBuffer&&, CountBuffer&&)>& consume,
        int threads, size_t chunk_size, size_t max_hit_bytes) const {
    if (queries.size() != other_queries.size())
        throw std::invalid_argument("paired query populations must have equal lengths");
    if (alphabet() != other.alphabet())
        throw std::invalid_argument("paired indexes must use the same alphabet");
    constexpr size_t columns = 6;
    const size_t max_rows = max_hit_bytes / (columns * sizeof(uint32_t));
    if (!chunk_size || !max_rows || !consume)
        throw std::invalid_argument("positive chunk_size/max_hit_bytes and a consumer required");
    if (queries.size() > uint64_t(std::numeric_limits<uint32_t>::max()) + 1)
        throw std::overflow_error("query indices exceed the uint32 output buffer");
    // Bound both the primitive output and the per-worker Hit vectors. A single
    // query may still exceed the explicit output cap; detect that without truncation.
    chunk_size = std::min(chunk_size, std::max(size_t(1), max_rows / std::max(size_t(1), size_t(size()) + size_t(other.size()))));
    const size_t blocks = queries.size() / chunk_size + (queries.size() % chunk_size != 0);
    std::atomic<bool> cancelled{false};
    std::mutex callback_mutex;
    struct Local { Searcher searcher, other_searcher; std::vector<std::vector<Hit>> hits, other_hits; };
    parallel_for(blocks, threads,
        [&] { return Local{Searcher(*this), Searcher(other), std::vector<std::vector<Hit>>(chunk_size), std::vector<std::vector<Hit>>(chunk_size)}; },
        [&](size_t block, Local& worker) {
            if (cancelled.load()) return;
            try {
                const size_t start = block * chunk_size;
                const size_t stop = start + std::min(chunk_size, queries.size() - start);
                size_t rows = 0;
                for (size_t i = start; i < stop; ++i) {
                    if (cancelled.load()) return;
                    auto& hits = worker.hits[i-start];
                    worker.searcher.search_into(queries[i], p, hits);
                    auto& other_hits = worker.other_hits[i-start];
                    worker.other_searcher.search_into(other_queries[i], p, other_hits);
                    if (hits.size() > max_rows - rows || other_hits.size() > max_rows - rows - hits.size())
                        throw std::overflow_error("complete query hits exceed max_hit_bytes; increase the explicit buffer limit");
                    rows += hits.size() + other_hits.size();
                }
                auto pack = [&](const auto& groups) {
                    CountBuffer buffer;
                    size_t count = 0;
                    for (size_t i = start; i < stop; ++i) count += groups[i-start].size();
                    buffer.shape = {count, columns};
                    buffer.data.reserve(count * columns);
                    for (size_t i = start; i < stop; ++i)
                        for (const auto& hit : groups[i-start]) {
                            const uint32_t values[] = {uint32_t(i), hit.ref_id, uint32_t(hit.score),
                                hit.n_subs, hit.n_ins, hit.n_dels};
                            buffer.data.insert(buffer.data.end(), std::begin(values), std::end(values));
                        }
                    return buffer;
                };
                auto buffer = pack(worker.hits), other_buffer = pack(worker.other_hits);
                std::lock_guard<std::mutex> lock(callback_mutex);
                if (!cancelled.load()) {
                    try { consume(start, stop, std::move(buffer), std::move(other_buffer)); }
                    catch (...) { cancelled.store(true); throw; }
                }
            } catch (...) {
                cancelled.store(true);
                throw;
            }
        });
}

std::vector<std::vector<uint64_t>> Index::edit_histogram_batch(
        const std::vector<std::string>& queries, const SearchParams& p,
        int threads, bool exclude_exact, int max_indels) const {
    if (max_indels < -1 || max_indels > 65535)
        throw std::invalid_argument("max_indels must be in [-1,65535]");
    if (!p.max_total_edits)
        throw std::invalid_argument("edit histogram requires explicit positive max_total_edits");
    if (p.engine == Engine::SeqTrie || p.mode != Mode::AllHits || p.max_hits)
        throw std::invalid_argument("edit histogram requires uncapped seqtm/auto all-hit search");
    std::vector<std::vector<uint64_t>> out(queries.size(),
        std::vector<uint64_t>(size_t(p.max_total_edits) + 1));
    struct Local { Searcher s; std::vector<Hit> hits; };
    parallel_for(queries.size(), threads, [&] { return Local{Searcher(*this), {}}; },
        [&](size_t i, Local& worker) {
            if (count_unit::fill_histogram(*this,queries[i],p,exclude_exact,out[i],max_indels)) return;
            worker.s.search_into(queries[i], p, worker.hits, max_indels);
            for (const auto& h : worker.hits) {
                const uint32_t edits = uint32_t(h.n_subs) + h.n_ins + h.n_dels;
                if (edits > p.max_total_edits)
                    throw std::logic_error("search hit exceeds max_total_edits");
                if (!exclude_exact || edits) ++out[i][edits];
            }
        });
    return out;
}

std::vector<uint64_t> Index::collisions_batch(const std::vector<std::string>& queries,
                                              const SearchParams& p, int threads) const {
    const size_t n = queries.size();
    std::vector<uint64_t> out(n, 0);
    if (n == 0) return out;

    struct Local { Searcher s; std::vector<Hit> scratch; };
    parallel_for(n, threads, [&] { return Local{Searcher(*this), {}}; },
                 [&](size_t i, Local& L) {
                     L.s.search_into(queries[i], p, L.scratch);
                     out[i] = L.s.last_collisions();
                 });
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

    // ref == a, qry == b: transpose to a-major with ref_id pointing into b. The hit describes
    // b[bi] -> a[ref_id]; reversing the direction swaps insertions and deletions.
    std::vector<std::vector<Hit>> out(a.size());
    for (uint32_t bi = 0; bi < res.size(); ++bi)
        for (const Hit& h : res[bi])
            out[h.ref_id].push_back(Hit{bi, h.score, h.n_subs, h.n_dels, h.n_ins});
    return out;
}

}  // namespace seqtree
