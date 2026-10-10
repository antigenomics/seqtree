#include "seqtree/seqtree.hpp"
#include "seqtree/parallel.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <tuple>

namespace seqtree {
namespace {

std::vector<uint8_t> encode(const Codec& c, const std::string& s, const char* label, size_t idx) {
    std::vector<uint8_t> out(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        uint8_t v = c.encode(s[i]);
        if (v == Codec::kInvalid)
            throw std::invalid_argument(label + std::string("[") + std::to_string(idx) + "] ('" + s +
                                        "'): symbol '" + s[i] + "' is not in the alphabet");
        out[i] = v;
    }
    return out;
}

// Penalty table flattened to [a * A + b]; unit cost when no matrix is given.
std::vector<int32_t> pen_table(const SubstitutionMatrix* m, uint8_t A) {
    std::vector<int32_t> pen(size_t(A) * A);
    for (uint8_t a = 0; a < A; ++a)
        for (uint8_t b = 0; b < A; ++b)
            pen[size_t(a) * A + b] = m ? m->penalty(a, b) : int32_t(a != b);
    return pen;
}

// One (query, ref) cell. `suf` is caller-owned scratch of length >= min(m,n)+1.
//
// prefix[i] + suffix[i] is the substitution cost of the layout whose gap block opens at
// column i of the shorter sequence, so a single forward sweep visits every block position.
// At d == 0 both diagonals coincide, prefix[i] + suffix[i] is constant, and the prior is
// contractually zero -- hence the early return, which also keeps s(q, q) == 0.
inline int32_t cell(const uint8_t* q, uint32_t m, const uint8_t* r, uint32_t n,
                    const int32_t* pen, uint8_t A, int32_t gap_open, int32_t gap_extend,
                    const int32_t* prior_row, int32_t* suf, const int32_t* weights = nullptr) {
    const uint32_t L = std::min(m, n);
    const uint32_t d = (m > n ? m - n : n - m);
    const bool q_longer = m >= n;

    suf[L] = 0;
    for (uint32_t j = L; j-- > 0;) {
        const uint8_t a = q_longer ? q[j + d] : q[j];
        const uint8_t b = q_longer ? r[j] : r[j + d];
        suf[j] = suf[j + 1] + pen[size_t(a) * A + b] * (weights ? weights[j + d] : 1);
    }
    if (d == 0) return suf[0];

    int32_t pre = 0;
    int32_t best = std::numeric_limits<int32_t>::max();
    for (uint32_t i = 0; i <= L; ++i) {
        const int32_t cand = pre + suf[i] + (prior_row ? prior_row[i] : 0);
        if (cand < best) best = cand;
        if (i < L) pre += pen[size_t(q[i]) * A + r[i]] * (weights ? weights[i] : 1);
    }
    return best + gap_open + int32_t(d - 1) * gap_extend;
}

// Shared preparation for bounded reductions: encoded inputs and prior remain
// shared read-only; no Q*N allocation. Validate before starting worker threads.
struct GapInput {
    std::vector<std::vector<uint8_t>> q, r;
    std::vector<int32_t> pen;
    std::vector<std::vector<int32_t>> weights;
    uint32_t longest = 0;
    uint8_t A = 0;
    size_t W1 = 0;
};

GapInput prepare(const std::vector<std::string>& queries, const std::vector<std::string>& refs,
                 Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
                 const std::vector<int32_t>& prior, uint32_t width, int threads,
                 const std::vector<std::vector<int32_t>>& weights = {}) {
    if (go < 0 || ge < 0 || threads < 0)
        throw std::invalid_argument("gap costs and threads must be >= 0");
    if (refs.size() > std::numeric_limits<uint32_t>::max())
        throw std::invalid_argument("too many reference IDs");
    Codec codec(alphabet);
    GapInput in;
    in.A = codec.size();
    if (matrix && matrix->size() != in.A)
        throw std::invalid_argument("matrix size does not match the alphabet");
    in.pen = pen_table(matrix, in.A);
    auto encode_all = [&](const auto& strings, auto& out, const char* label) {
        out.reserve(strings.size());
        for (size_t i = 0; i < strings.size(); ++i) {
            if (strings[i].size() > std::numeric_limits<uint32_t>::max())
                throw std::invalid_argument("sequence length exceeds uint32");
            out.push_back(encode(codec, strings[i], label, i));
            in.longest = std::max(in.longest, uint32_t(strings[i].size()));
        }
    };
    encode_all(queries, in.q, "queries");
    encode_all(refs, in.r, "refs");
    in.W1 = size_t(width) + 1;
    int32_t largest_prior = 0;
    if (!prior.empty()) {
        if (in.W1 > std::numeric_limits<size_t>::max() / in.W1 ||
            in.W1 * in.W1 > std::numeric_limits<size_t>::max() / in.W1 ||
            prior.size() != in.W1 * in.W1 * in.W1)
            throw std::invalid_argument("prior table must have (prior_width + 1)^3 entries");
        if (in.longest > width)
            throw std::invalid_argument("sequence exceeds prior table width");
        for (auto value : prior) {
            if (value < 0) throw std::invalid_argument("gap prior must be nonnegative");
            largest_prior = std::max(largest_prior, value);
        }
    }
    int32_t largest_weight = 1;
    if (!weights.empty()) {
        in.weights = weights; // one immutable shared table, never copied per worker
        for (size_t length = 0; length < weights.size(); ++length) {
            if (weights[length].empty()) continue;
            if (weights[length].size() != length)
                throw std::invalid_argument("positional weight length must match its frame");
            for (const auto value : weights[length]) {
                if (value < 0) throw std::invalid_argument("positional weights must be nonnegative");
                largest_weight = std::max(largest_weight, value);
            }
        }
        for (const auto* axis : {&in.q, &in.r}) for (const auto& sequence : *axis) {
            const auto length = sequence.size();
            if (length && (length >= weights.size() || weights[length].size() != length))
                throw std::invalid_argument("missing positional weights for an input frame");
        }
    }
    const auto largest_pen = *std::max_element(in.pen.begin(), in.pen.end());
    if (*std::min_element(in.pen.begin(), in.pen.end()) < 0)
        throw std::invalid_argument("matrix penalties must be nonnegative");
    const uint64_t weighted_pen = uint64_t(uint32_t(largest_pen)) * uint32_t(largest_weight);
    if (in.longest && weighted_pen > uint64_t(std::numeric_limits<int32_t>::max()) / in.longest)
        throw std::invalid_argument("gapblock weighted score can overflow int32");
    const uint64_t bound = uint64_t(in.longest) * weighted_pen +
                           uint32_t(largest_prior) + uint32_t(go) +
                           uint64_t(in.longest ? in.longest - 1 : 0) * uint32_t(ge);
    if (bound > std::numeric_limits<int32_t>::max())
        throw std::invalid_argument("gapblock score can overflow int32");
    return in;
}

int32_t score_cell(const GapInput& in, size_t i, size_t j, int32_t go, int32_t ge,
                   const std::vector<int32_t>& prior, int32_t* scratch) {
    const uint32_t m = uint32_t(in.q[i].size()), n = uint32_t(in.r[j].size());
    const uint32_t M = std::max(m, n), d = m > n ? m - n : n - m;
    const int32_t* prow = prior.empty() ? nullptr : prior.data() + (size_t(M) * in.W1 + d) * in.W1;
    return cell(in.q[i].data(), m, in.r[j].data(), n, in.pen.data(), in.A, go, ge, prow, scratch, in.weights.empty() || !M ? nullptr : in.weights[M].data());
}

struct Ranked {
    uint32_t ref_id;
    int32_t a, b;
};
bool better(const Ranked& a, const Ranked& b) {
    return std::tuple(std::max(a.a, a.b), a.ref_id) < std::tuple(std::max(b.a, b.b), b.ref_id);
}

struct Worker {
    std::vector<int32_t> scratch;
    std::vector<Ranked> heap;
};

std::vector<std::vector<Ranked>> ranked(const GapInput& a, const GapInput* b, uint32_t k,
                                      int32_t go, int32_t ge, const std::vector<int32_t>& prior,
                                      int threads, bool exclude_exact) {
    if (!k) throw std::invalid_argument("k must be positive");
    if (b && (b->q.size() != a.q.size() || b->r.size() != a.r.size()))
        throw std::invalid_argument("paired query/reference axes must agree");
    k = uint32_t(std::min<size_t>(k, a.r.size()));
    std::vector<std::vector<Ranked>> out(a.q.size());
    if (!k) return out;
    const auto longest = b ? std::max(a.longest, b->longest) : a.longest;
    parallel_for(a.q.size(), threads, [&] {
        Worker w{std::vector<int32_t>(size_t(longest) + 1), {}};
        w.heap.reserve(k);
        return w;
    }, [&](size_t i, Worker& w) {
        w.heap.clear();
        for (size_t j = 0; j < a.r.size(); ++j) {
            if (exclude_exact && a.q[i] == a.r[j] && (!b || b->q[i] == b->r[j])) continue;
            const auto sa = score_cell(a, i, j, go, ge, prior, w.scratch.data());
            const auto sb = b ? score_cell(*b, i, j, go, ge, prior, w.scratch.data()) : sa;
            const Ranked h{uint32_t(j), sa, sb};
            if (w.heap.size() < k) {
                w.heap.push_back(h);
                std::push_heap(w.heap.begin(), w.heap.end(), better);
            } else if (better(h, w.heap.front())) {
                std::pop_heap(w.heap.begin(), w.heap.end(), better);
                w.heap.back() = h;
                std::push_heap(w.heap.begin(), w.heap.end(), better);
            }
        }
        std::sort_heap(w.heap.begin(), w.heap.end(), better);
        out[i] = w.heap;
    });
    return out;
}

std::vector<std::vector<uint64_t>> counted(const GapInput& a, const GapInput* b,
                  const std::vector<std::vector<int32_t>>& thresholds, int32_t go, int32_t ge,
                  const std::vector<int32_t>& prior, int threads, bool exclude_exact) {
    if (thresholds.size() != a.q.size())
        throw std::invalid_argument("thresholds must have one row per query");
    if (b && (b->q.size() != a.q.size() || b->r.size() != a.r.size()))
        throw std::invalid_argument("paired query/reference axes must agree");
    std::vector<std::vector<uint64_t>> out(a.q.size());
    const auto longest = b ? std::max(a.longest, b->longest) : a.longest;
    parallel_for(a.q.size(), threads, [&] { return std::vector<int32_t>(size_t(longest) + 1); },
        [&](size_t i, std::vector<int32_t>& scratch) {
            std::vector<std::pair<int32_t, size_t>> sorted;
            sorted.reserve(thresholds[i].size());
            for (size_t t = 0; t < thresholds[i].size(); ++t) sorted.emplace_back(thresholds[i][t], t);
            std::sort(sorted.begin(), sorted.end());
            std::vector<uint64_t> bins(sorted.size());
            out[i].resize(sorted.size());
            if (sorted.empty()) return;
            for (size_t j = 0; j < a.r.size(); ++j) {
                if (exclude_exact && a.q[i] == a.r[j] && (!b || b->q[i] == b->r[j])) continue;
                auto score = score_cell(a, i, j, go, ge, prior, scratch.data());
                if (b) score = std::max(score, score_cell(*b, i, j, go, ge, prior, scratch.data()));
                const auto pos = std::lower_bound(sorted.begin(), sorted.end(), score,
                    [](const auto& threshold, int32_t value) { return threshold.first < value; });
                if (pos != sorted.end()) ++bins[size_t(pos - sorted.begin())];
            }
            uint64_t cumulative = 0;
            for (size_t t = 0; t < sorted.size(); ++t) {
                cumulative += bins[t];
                out[i][sorted[t].second] = cumulative;
            }
        });
    return out;
}

}  // namespace

std::vector<int32_t> gapblock_matrix(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, Alphabet alphabet,
        const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads,
        const std::vector<std::vector<int32_t>>& weights) {
    const auto in = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, std::max(threads, 0), weights);
    if (!refs.empty() && queries.size() > std::numeric_limits<size_t>::max() / refs.size())
        throw std::invalid_argument("gapblock matrix size overflows size_t");
    std::vector<int32_t> out(queries.size() * refs.size());
    parallel_for(queries.size(), threads, [&] { return std::vector<int32_t>(size_t(in.longest) + 1); },
        [&](size_t i, std::vector<int32_t>& scratch) {
            for (size_t j = 0; j < refs.size(); ++j)
                out[i * refs.size() + j] = score_cell(in, i, j, go, ge, prior, scratch.data());
        });
    return out;
}

std::vector<std::vector<Hit>> gapblock_ball_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, const std::vector<int32_t>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights) {
    if (thresholds.size() != queries.size())
        throw std::invalid_argument("ball thresholds must have one value per query");
    const auto in = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, threads, weights);
    std::vector<std::vector<Hit>> out(queries.size());
    parallel_for(queries.size(), threads, [&] { return std::vector<int32_t>(size_t(in.longest) + 1); },
        [&](size_t i, std::vector<int32_t>& scratch) {
            if (thresholds[i] < 0) return;
            for (size_t j = 0; j < refs.size(); ++j) {
                if (exclude_exact && in.q[i] == in.r[j]) continue;
                const auto score = score_cell(in, i, j, go, ge, prior, scratch.data());
                if (score <= thresholds[i]) out[i].push_back(Hit{uint32_t(j), score, 0, 0, 0});
            }
        });
    return out;
}

std::vector<std::vector<Hit>> gapblock_topk_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, uint32_t k, Alphabet alphabet,
        const SubstitutionMatrix* matrix, int32_t go, int32_t ge, const std::vector<int32_t>& prior,
        uint32_t width, int threads, bool exclude_exact) {
    const auto input = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, threads);
    const auto rows = ranked(input, nullptr, k, go, ge, prior, threads, exclude_exact);
    std::vector<std::vector<Hit>> out(rows.size());
    for (size_t i = 0; i < rows.size(); ++i) {
        out[i].reserve(rows[i].size());
        for (auto h : rows[i]) out[i].push_back(Hit{h.ref_id, h.a});
    }
    return out;
}

std::vector<std::vector<std::tuple<uint32_t, int32_t, int32_t>>> gapblock_paired_topk_batch(
        const std::vector<std::string>& qa, const std::vector<std::string>& qb,
        const std::vector<std::string>& ra, const std::vector<std::string>& rb,
        uint32_t k, Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact) {
    const auto a = prepare(qa, ra, alphabet, matrix, go, ge, prior, width, threads);
    const auto b = prepare(qb, rb, alphabet, matrix, go, ge, prior, width, threads);
    const auto rows = ranked(a, &b, k, go, ge, prior, threads, exclude_exact);
    std::vector<std::vector<std::tuple<uint32_t, int32_t, int32_t>>> out(rows.size());
    for (size_t i = 0; i < rows.size(); ++i)
        for (auto h : rows[i]) out[i].emplace_back(h.ref_id, h.a, h.b);
    return out;
}

std::vector<std::vector<uint64_t>> gapblock_count_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, const std::vector<std::vector<int32_t>>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights) {
    const auto a = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, threads, weights);
    return counted(a, nullptr, thresholds, go, ge, prior, threads, exclude_exact);
}

std::vector<std::vector<uint64_t>> gapblock_paired_count_batch(
        const std::vector<std::string>& qa, const std::vector<std::string>& qb,
        const std::vector<std::string>& ra, const std::vector<std::string>& rb,
        const std::vector<std::vector<int32_t>>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact) {
    const auto a = prepare(qa, ra, alphabet, matrix, go, ge, prior, width, threads);
    const auto b = prepare(qb, rb, alphabet, matrix, go, ge, prior, width, threads);
    return counted(a, &b, thresholds, go, ge, prior, threads, exclude_exact);
}

}  // namespace seqtree
