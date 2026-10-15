#include "seqtree/seqtree.hpp"
#include "seqtree/parallel.hpp"
#include "pairwise_detail.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

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
                    const int32_t* prior_row, int32_t* suf, const int32_t* weights = nullptr,
                    const int32_t* qw = nullptr, const int32_t* rw = nullptr) {
    const uint32_t L = std::min(m, n);
    const uint32_t d = (m > n ? m - n : n - m);
    const bool q_longer = m >= n;

    suf[L] = 0;
    for (uint32_t j = L; j-- > 0;) {
        const uint8_t a = q_longer ? q[j + d] : q[j];
        const uint8_t b = q_longer ? r[j] : r[j + d];
        suf[j] = suf[j + 1] + pen[size_t(a) * A + b] * (qw && rw ? std::max(qw[j + (q_longer ? d : 0)], rw[j + (q_longer ? 0 : d)]) : (weights ? weights[j + d] : 1));
    }
    if (d == 0) return suf[0];

    int32_t pre = 0;
    int32_t best = std::numeric_limits<int32_t>::max();
    for (uint32_t i = 0; i <= L; ++i) {
        const int32_t cand = pre + suf[i] + (prior_row ? prior_row[i] : 0);
        if (cand < best) best = cand;
        if (i < L) pre += pen[size_t(q[i]) * A + r[i]] * (qw && rw ? std::max(qw[i], rw[i]) : (weights ? weights[i] : 1));
    }
    return best + gap_open + int32_t(d - 1) * gap_extend;
}

// Shared preparation for bounded reductions: encoded inputs and prior remain
// shared read-only; no Q*N allocation. Validate before starting worker threads.
struct GapInput {
    std::vector<std::vector<uint8_t>> q, r;
    std::vector<int32_t> pen;
    std::vector<std::vector<int32_t>> weights, query_weights, reference_weights;
    std::vector<int32_t> query_groups, reference_groups;
    std::vector<std::vector<int32_t>> group_distances;
    std::optional<PairedFeatureOffsets> features;
    uint32_t longest = 0;
    uint8_t A = 0;
    size_t W1 = 0;
};

// Explicit ties-even rounding is independent of the process floating rounding mode.
inline double round_even(double value) {
    const double base = std::floor(value), fraction = value - base;
    return base + (fraction > 0.5 || (fraction == 0.5 && std::fmod(base, 2.0) != 0.0));
}

uint64_t validate_feature(const PairedFeatureOffsets& f, size_t query_count, size_t reference_count) {
    if (!std::isfinite(f.weight) || f.weight < 0 || f.scale < 0)
        throw std::invalid_argument("feature weight and scale must be finite and nonnegative");
    if (f.query_codes.size() != query_count || f.reference_codes.size() != reference_count)
        throw std::invalid_argument("feature code rows must match query/reference axes");
    if (f.costs.empty() || f.costs.front().empty())
        throw std::invalid_argument("feature costs must be a nonempty rectangular matrix");
    const size_t columns = f.costs.front().size();
    int32_t maximum = 0;
    for (const auto& row : f.costs) {
        if (row.size() != columns) throw std::invalid_argument("feature costs must be rectangular");
        for (auto cost : row) {
            if (cost < 0) throw std::invalid_argument("feature costs must be nonnegative");
            maximum = std::max(maximum, cost);
        }
    }
    const size_t width = !f.query_codes.empty() ? f.query_codes.front().size() :
                         (!f.reference_codes.empty() ? f.reference_codes.front().size() : 0);
    if (width % 2 || (!width && (query_count!=0 || reference_count!=0)))
        throw std::invalid_argument("feature rows must contain a positive number of code pairs");
    auto validate = [&](const auto& rows, size_t bound) {
        for (const auto& row : rows) {
            if (row.size() != width) throw std::invalid_argument("feature rows must have equal depth");
            for (auto code : row)
                if (code < 0 || size_t(code) >= bound)
                    throw std::invalid_argument("feature code outside costs");
        }
    };
    validate(f.query_codes, f.costs.size());
    validate(f.reference_codes, columns);
    const size_t depth = width / 2;
    if (depth && uint64_t(maximum) > std::numeric_limits<uint64_t>::max() / depth)
        throw std::invalid_argument("feature cost sum can overflow uint64");
    const double bound = depth && f.weight && f.scale ? round_even(f.weight * double(uint64_t(maximum) * depth) / double(depth)) * f.scale : 0;
    if (!std::isfinite(bound) || bound > std::numeric_limits<int32_t>::max())
        throw std::invalid_argument("feature offset can overflow int32");
    return uint64_t(bound);
}

GapInput prepare(const std::vector<std::string>& queries, const std::vector<std::string>& refs,
                 Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
                 const std::vector<int32_t>& prior, uint32_t width, int threads,
                 const std::vector<std::vector<int32_t>>& weights = {},
                 const std::vector<int32_t>& query_groups = {},
                 const std::vector<int32_t>& reference_groups = {},
                 const std::vector<std::vector<int32_t>>& group_distances = {},
                 const std::vector<std::vector<int32_t>>& query_weights = {},
                 const std::vector<std::vector<int32_t>>& reference_weights = {},
                 const std::optional<PairedFeatureOffsets>& feature_offsets = std::nullopt,
                 uint64_t* score_bound = nullptr) {
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
    int32_t largest_weight = weights.empty() ? 1 : 0;
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
    if (!query_weights.empty() || !reference_weights.empty()) {
        if (!weights.empty()) throw std::invalid_argument("length and per-input weights cannot be combined");
        auto validate = [&](const auto& rows, const auto& axis) {
            if (rows.size() != axis.size())
                throw std::invalid_argument("per-input weight rows must match their sequence axis");
            for (size_t i = 0; i < axis.size(); ++i) {
                if (rows[i].size() != axis[i].size())
                    throw std::invalid_argument("per-input weight length must match its sequence");
                for (auto value : rows[i]) {
                    if (value < 0) throw std::invalid_argument("positional weights must be nonnegative");
                    largest_weight = std::max(largest_weight, value);
                }
            }
        };
        largest_weight = 0;
        validate(query_weights, in.q);
        validate(reference_weights, in.r);
        in.query_weights = query_weights;
        in.reference_weights = reference_weights;
    }
    int32_t largest_offset = 0;
    if (!query_groups.empty() || !reference_groups.empty() || !group_distances.empty()) {
        if (query_groups.size() != queries.size() || reference_groups.size() != refs.size())
            throw std::invalid_argument("group IDs must match their query/reference axes");
        if (group_distances.empty() || group_distances.front().empty())
            throw std::invalid_argument("group distances must be a nonempty rectangular matrix");
        const auto columns = group_distances.front().size();
        for (const auto& row : group_distances) {
            if (row.size() != columns)
                throw std::invalid_argument("group distances must be rectangular");
            for (const auto value : row) {
                if (value < 0) throw std::invalid_argument("group distances must be nonnegative");
                largest_offset = std::max(largest_offset, value);
            }
        }
        for (const auto id : query_groups)
            if (id < 0 || size_t(id) >= group_distances.size())
                throw std::invalid_argument("query group ID outside group distances");
        for (const auto id : reference_groups)
            if (id < 0 || size_t(id) >= columns)
                throw std::invalid_argument("reference group ID outside group distances");
        in.query_groups = query_groups;
        in.reference_groups = reference_groups;
        in.group_distances = group_distances;
    }
    const uint64_t largest_feature_offset = feature_offsets ?
        validate_feature(*feature_offsets,queries.size(),refs.size()) : 0;
    if (feature_offsets) in.features=*feature_offsets;
    const auto largest_pen = *std::max_element(in.pen.begin(), in.pen.end());
    if (*std::min_element(in.pen.begin(), in.pen.end()) < 0)
        throw std::invalid_argument("matrix penalties must be nonnegative");
    const uint64_t weighted_pen = uint64_t(uint32_t(largest_pen)) * uint32_t(largest_weight);
    if (in.longest && weighted_pen > uint64_t(std::numeric_limits<int32_t>::max()) / in.longest)
        throw std::invalid_argument("gapblock weighted score can overflow int32");
    const uint64_t bound = uint64_t(in.longest) * weighted_pen +
                           uint32_t(largest_prior) + uint32_t(largest_offset) + largest_feature_offset + uint32_t(go) +
                           uint64_t(in.longest ? in.longest - 1 : 0) * uint32_t(ge);
    if (bound > std::numeric_limits<int32_t>::max())
        throw std::invalid_argument("gapblock score can overflow int32");
    if (score_bound) *score_bound=bound;
    return in;
}

inline int32_t feature_offset(const PairedFeatureOffsets& f,size_t i,size_t j) {
    if (f.weight==0 || f.scale==0) return 0;
    const auto& q=f.query_codes[i];const auto& r=f.reference_codes[j];
    uint64_t sum=0;
    for (size_t k=0;k<q.size();k+=2)
        sum+=std::max(f.costs[q[k]][r[k]],f.costs[q[k+1]][r[k+1]]);
    return int32_t(round_even(f.weight*double(sum)/double(q.size()/2))*f.scale);
}

template<bool Factored = true>
inline int32_t group_offset(const GapInput& in, size_t i, size_t j) {
    int32_t offset = in.group_distances.empty() ? 0 : in.group_distances[in.query_groups[i]][in.reference_groups[j]];
    if constexpr (Factored) if (in.features) offset+=feature_offset(*in.features,i,j);
    return offset;
}

int32_t score_cell(const GapInput& in, size_t i, size_t j, int32_t go, int32_t ge,
                   const std::vector<int32_t>& prior, int32_t* scratch, int32_t offset) {
    const uint32_t m = uint32_t(in.q[i].size()), n = uint32_t(in.r[j].size());
    const uint32_t M = std::max(m, n), d = m > n ? m - n : n - m;
    const int32_t* prow = prior.empty() ? nullptr : prior.data() + (size_t(M) * in.W1 + d) * in.W1;
    return offset + cell(in.q[i].data(), m, in.r[j].data(), n, in.pen.data(), in.A, go, ge, prow, scratch, in.weights.empty() || !M ? nullptr : in.weights[M].data(),
        in.query_weights.empty() ? nullptr : in.query_weights[i].data(),
        in.reference_weights.empty() ? nullptr : in.reference_weights[j].data());
}

template<bool Factored = true>
int32_t score_cell(const GapInput& in, size_t i, size_t j, int32_t go, int32_t ge,
                   const std::vector<int32_t>& prior, int32_t* scratch) {
    return score_cell(in,i,j,go,ge,prior,scratch,group_offset<Factored>(in,i,j));
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

template<bool Factored>
std::vector<std::vector<uint64_t>> counted_impl(const GapInput& a, const GapInput* b,
                  const std::vector<std::vector<int32_t>>& thresholds, int32_t go, int32_t ge,
                  const std::vector<int32_t>& prior, int threads, bool exclude_exact,
                  std::vector<std::vector<uint64_t>>* masses = nullptr, bool sum_pair_distance = false) {
    if (thresholds.size() != a.q.size())
        throw std::invalid_argument("thresholds must have one row per query");
    if (b && (b->q.size() != a.q.size() || b->r.size() != a.r.size()))
        throw std::invalid_argument("paired query/reference axes must agree");
    if (sum_pair_distance && !b)
        throw std::invalid_argument("sum distance requires two linked views");
    std::vector<std::vector<uint64_t>> out(a.q.size());
    if (masses) {
        masses->resize(a.q.size());
        if (std::any_of(a.pen.begin(), a.pen.end(), [](int32_t cost) { return cost < 0; }) ||
            (b && std::any_of(b->pen.begin(), b->pen.end(), [](int32_t cost) { return cost < 0; })))
            throw std::invalid_argument("linear mass requires nonnegative scores");
        // Every score admitted to a row is <= its largest positive threshold.
        // Bound both score sums and threshold*count before workers start.
        for (const auto& row : thresholds) for (auto threshold : row)
            if (threshold > 0 && a.r.size() > std::numeric_limits<uint64_t>::max() / uint64_t(threshold))
                throw std::overflow_error("gapblock linear mass exceeds uint64");
    }
    const auto longest = b ? std::max(a.longest, b->longest) : a.longest;
    parallel_for(a.q.size(), threads, [&] { return std::vector<int32_t>(size_t(longest) + 1); },
        [&](size_t i, std::vector<int32_t>& scratch) {
            std::vector<std::pair<int32_t, size_t>> sorted;
            sorted.reserve(thresholds[i].size());
            for (size_t t = 0; t < thresholds[i].size(); ++t) sorted.emplace_back(thresholds[i][t], t);
            std::sort(sorted.begin(), sorted.end());
            std::vector<uint64_t> bins(sorted.size());
            std::vector<uint64_t> score_bins(masses ? sorted.size() : 0);
            out[i].resize(sorted.size());
            if (masses) (*masses)[i].resize(sorted.size());
            if (sorted.empty()) return;
            for (size_t j = 0; j < a.r.size(); ++j) {
                if (exclude_exact && a.q[i] == a.r[j] && (!b || b->q[i] == b->r[j])) continue;
                const int64_t offset = int64_t(group_offset<Factored>(a, i, j)) +
                    (sum_pair_distance ? int64_t(group_offset<Factored>(*b, i, j)) : 0);
                if (offset > sorted.back().first) continue;
                int64_t score = score_cell<Factored>(a, i, j, go, ge, prior, scratch.data());
                if (sum_pair_distance) {
                    if (score > sorted.back().first) continue;
                    score += int64_t(score_cell<Factored>(*b, i, j, go, ge, prior, scratch.data()));
                } else if (b) {
                    score = std::max<int64_t>(score, score_cell<Factored>(*b, i, j, go, ge, prior, scratch.data()));
                }
                const auto pos = std::lower_bound(sorted.begin(), sorted.end(), score,
                    [](const auto& threshold, int64_t value) { return threshold.first < value; });
                if (pos != sorted.end()) {
                    const auto bin = size_t(pos - sorted.begin());
                    ++bins[bin];
                    if (masses) score_bins[bin] += uint64_t(score);
                }
            }
            uint64_t cumulative = 0, score_sum = 0;
            for (size_t t = 0; t < sorted.size(); ++t) {
                cumulative += bins[t];
                out[i][sorted[t].second] = cumulative;
                if (masses) {
                    score_sum += score_bins[t];
                    (*masses)[i][sorted[t].second] = uint64_t(std::max(sorted[t].first, 0)) * cumulative - score_sum;
                }
            }
        });
    return out;
}

template<bool Factored>
std::vector<std::vector<uint64_t>> counted_shared_impl(const GapInput& a, const GapInput* b,
                  const std::vector<std::vector<int32_t>>& thresholds, int32_t go, int32_t ge,
                  const std::vector<int32_t>& prior, int threads, bool exclude_exact,
                  std::vector<std::vector<uint64_t>>* masses = nullptr, bool sum_pair_distance = false) {
    if (thresholds.size() != 1)
        throw std::invalid_argument("shared thresholds must contain exactly one row");
    if (b && (b->q.size() != a.q.size() || b->r.size() != a.r.size()))
        throw std::invalid_argument("paired query/reference axes must agree");
    if (sum_pair_distance && !b)
        throw std::invalid_argument("sum distance requires two linked views");
    std::vector<std::vector<uint64_t>> out(a.q.size());
    if (masses) {
        masses->resize(a.q.size());
        if (std::any_of(a.pen.begin(), a.pen.end(), [](int32_t cost) { return cost < 0; }) ||
            (b && std::any_of(b->pen.begin(), b->pen.end(), [](int32_t cost) { return cost < 0; })))
            throw std::invalid_argument("linear mass requires nonnegative scores");
        // Every score admitted to a row is <= its largest positive threshold.
        // Bound both score sums and threshold*count before workers start.
        for (const auto& row : thresholds) for (auto threshold : row)
            if (threshold > 0 && a.r.size() > std::numeric_limits<uint64_t>::max() / uint64_t(threshold))
                throw std::overflow_error("gapblock linear mass exceeds uint64");
    }
    auto sort_row = [](auto& sorted, const auto& row) {
        sorted.reserve(row.size());
        for (size_t t = 0; t < row.size(); ++t) sorted.emplace_back(row[t], t);
        std::sort(sorted.begin(), sorted.end());
    };
    std::vector<std::pair<int32_t, size_t>> shared_sorted;
    sort_row(shared_sorted, thresholds.front());
    const auto longest = b ? std::max(a.longest, b->longest) : a.longest;
    parallel_for(a.q.size(), threads, [&] { return std::vector<int32_t>(size_t(longest) + 1); },
        [&](size_t i, std::vector<int32_t>& scratch) {
            const auto& sorted = shared_sorted;
            std::vector<uint64_t> bins(sorted.size());
            std::vector<uint64_t> score_bins(masses ? sorted.size() : 0);
            out[i].resize(sorted.size());
            if (masses) (*masses)[i].resize(sorted.size());
            if (sorted.empty()) return;
            for (size_t j = 0; j < a.r.size(); ++j) {
                if (exclude_exact && a.q[i] == a.r[j] && (!b || b->q[i] == b->r[j])) continue;
                const int64_t offset = int64_t(group_offset<Factored>(a, i, j)) +
                    (sum_pair_distance ? int64_t(group_offset<Factored>(*b, i, j)) : 0);
                if (offset > sorted.back().first) continue;
                int64_t score = score_cell<Factored>(a, i, j, go, ge, prior, scratch.data());
                if (sum_pair_distance) {
                    if (score > sorted.back().first) continue;
                    score += int64_t(score_cell<Factored>(*b, i, j, go, ge, prior, scratch.data()));
                } else if (b) {
                    score = std::max<int64_t>(score, score_cell<Factored>(*b, i, j, go, ge, prior, scratch.data()));
                }
                const auto pos = std::lower_bound(sorted.begin(), sorted.end(), score,
                    [](const auto& threshold, int64_t value) { return threshold.first < value; });
                if (pos != sorted.end()) {
                    const auto bin = size_t(pos - sorted.begin());
                    ++bins[bin];
                    if (masses) score_bins[bin] += uint64_t(score);
                }
            }
            uint64_t cumulative = 0, score_sum = 0;
            for (size_t t = 0; t < sorted.size(); ++t) {
                cumulative += bins[t];
                out[i][sorted[t].second] = cumulative;
                if (masses) {
                    score_sum += score_bins[t];
                    (*masses)[i][sorted[t].second] = uint64_t(std::max(sorted[t].first, 0)) * cumulative - score_sum;
                }
            }
        });
    return out;
}

// Dispatch once per batch: absent features retain the original hot loop.
std::vector<std::vector<uint64_t>> counted(const GapInput& a, const GapInput* b,
                  const std::vector<std::vector<int32_t>>& thresholds, int32_t go, int32_t ge,
                  const std::vector<int32_t>& prior, int threads, bool exclude_exact,
                  std::vector<std::vector<uint64_t>>* masses = nullptr, bool sum_pair_distance = false, bool shared_thresholds = false) {
    if (shared_thresholds) {
        if (!a.features && (!b || !b->features))
            return counted_shared_impl<false>(a, b, thresholds, go, ge, prior, threads, exclude_exact, masses, sum_pair_distance);
        return counted_shared_impl<true>(a, b, thresholds, go, ge, prior, threads, exclude_exact, masses, sum_pair_distance);
    }
    if (!a.features && (!b || !b->features))
        return counted_impl<false>(a, b, thresholds, go, ge, prior, threads, exclude_exact, masses, sum_pair_distance);
    return counted_impl<true>(a, b, thresholds, go, ge, prior, threads, exclude_exact, masses, sum_pair_distance);
}

}  // namespace

std::vector<int32_t> gapblock_matrix(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, Alphabet alphabet,
        const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads,
        const std::vector<std::vector<int32_t>>& weights,
        const std::vector<int32_t>& query_groups, const std::vector<int32_t>& reference_groups,
        const std::vector<std::vector<int32_t>>& group_distances, const std::vector<std::vector<int32_t>>& query_weights, const std::vector<std::vector<int32_t>>& reference_weights,
        const std::optional<PairedFeatureOffsets>& feature_offsets) {
    const auto in = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, std::max(threads, 0), weights, query_groups, reference_groups, group_distances, query_weights, reference_weights, feature_offsets);
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
        const std::vector<std::vector<int32_t>>& weights,
        const std::vector<int32_t>& query_groups, const std::vector<int32_t>& reference_groups,
        const std::vector<std::vector<int32_t>>& group_distances, const std::vector<std::vector<int32_t>>& query_weights, const std::vector<std::vector<int32_t>>& reference_weights,
        const std::optional<PairedFeatureOffsets>& feature_offsets) {
    if (thresholds.size() != queries.size())
        throw std::invalid_argument("ball thresholds must have one value per query");
    const auto in = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, threads, weights, query_groups, reference_groups, group_distances, query_weights, reference_weights, feature_offsets);
    std::vector<std::vector<Hit>> out(queries.size());
    parallel_for(queries.size(), threads, [&] { return std::vector<int32_t>(size_t(in.longest) + 1); },
        [&](size_t i, std::vector<int32_t>& scratch) {
            if (thresholds[i] < 0) return;
            for (size_t j = 0; j < refs.size(); ++j) {
                if (exclude_exact && in.q[i] == in.r[j]) continue;
                if (group_offset(in, i, j) > thresholds[i]) continue;
                const auto score = score_cell(in, i, j, go, ge, prior, scratch.data());
                if (score <= thresholds[i]) out[i].push_back(Hit{uint32_t(j), score, 0, 0, 0});
            }
        });
    return out;
}

std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<int32_t>>
gapblock_group_audit_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, const std::vector<int32_t>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights,
        const std::vector<int32_t>& query_groups, const std::vector<int32_t>& reference_groups,
        const std::vector<std::vector<int32_t>>& group_distances,
        const std::vector<std::vector<int32_t>>& query_weights,
        const std::vector<std::vector<int32_t>>& reference_weights,
        const std::optional<PairedFeatureOffsets>& feature_offsets) {
    if (thresholds.size() != queries.size())
        throw std::invalid_argument("group audit thresholds must have one value per query");
    const auto in = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, threads,
                           weights, query_groups, reference_groups, group_distances,
                           query_weights, reference_weights, feature_offsets);
    std::tuple<std::vector<uint64_t>, std::vector<uint64_t>, std::vector<int32_t>> out{
        std::vector<uint64_t>(queries.size()), std::vector<uint64_t>(queries.size()),
        std::vector<int32_t>(queries.size(), -1)};
    auto& [before, after, nearest] = out;
    parallel_for(queries.size(), threads,
        [&] { return std::vector<int32_t>(size_t(in.longest) + 1); },
        [&](size_t i, std::vector<int32_t>& scratch) {
            if (thresholds[i] < 0) return;
            int32_t minimum = std::numeric_limits<int32_t>::max();
            for (size_t j = 0; j < refs.size(); ++j) {
                if (exclude_exact && in.q[i] == in.r[j]) continue;
                const auto score = score_cell(in, i, j, go, ge, prior, scratch.data());
                if (score - group_offset(in, i, j) > thresholds[i]) continue;
                ++before[i];
                after[i] += score <= thresholds[i];
                minimum = std::min(minimum, score);
            }
            if (before[i]) nearest[i] = minimum;
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
        const std::vector<std::vector<int32_t>>& weights,
        const std::vector<int32_t>& query_groups, const std::vector<int32_t>& reference_groups,
        const std::vector<std::vector<int32_t>>& group_distances, const std::vector<std::vector<int32_t>>& query_weights, const std::vector<std::vector<int32_t>>& reference_weights,
        const std::optional<PairedFeatureOffsets>& feature_offsets, bool shared_thresholds) {
    const auto a = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, threads, weights, query_groups, reference_groups, group_distances, query_weights, reference_weights, feature_offsets);
    return counted(a, nullptr, thresholds, go, ge, prior, threads, exclude_exact, nullptr, false, shared_thresholds);
}

std::pair<std::vector<std::vector<uint64_t>>, std::vector<std::vector<uint64_t>>> gapblock_count_mass_batch(const std::vector<std::string>& queries,
        const std::vector<std::string>& refs, const std::vector<std::vector<int32_t>>& thresholds,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        const std::vector<int32_t>& prior, uint32_t width, int threads, bool exclude_exact,
        const std::vector<std::vector<int32_t>>& weights,
        const std::vector<int32_t>& query_groups, const std::vector<int32_t>& reference_groups,
        const std::vector<std::vector<int32_t>>& group_distances, const std::vector<std::vector<int32_t>>& query_weights, const std::vector<std::vector<int32_t>>& reference_weights,
        const std::optional<PairedFeatureOffsets>& feature_offsets, bool shared_thresholds) {
    const auto a = prepare(queries, refs, alphabet, matrix, go, ge, prior, width, threads, weights, query_groups, reference_groups, group_distances, query_weights, reference_weights, feature_offsets);
    std::vector<std::vector<uint64_t>> masses;
    auto counts = counted(a, nullptr, thresholds, go, ge, prior, threads, exclude_exact, &masses, false, shared_thresholds);
    return {std::move(counts), std::move(masses)};
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
        const std::vector<std::vector<int32_t>>& query_weights_alpha,
        const std::vector<std::vector<int32_t>>& query_weights_beta,
        const std::vector<std::vector<int32_t>>& reference_weights_alpha,
        const std::vector<std::vector<int32_t>>& reference_weights_beta,
        const std::optional<PairedFeatureOffsets>& feature_offsets_alpha,
        const std::optional<PairedFeatureOffsets>& feature_offsets_beta, bool shared_thresholds) {
    const auto a = prepare(qa, ra, alphabet, matrix, go, ge, prior, width, threads,
                           weights_alpha, query_groups_alpha, reference_groups_alpha, group_distances_alpha, query_weights_alpha, reference_weights_alpha, feature_offsets_alpha);
    const auto b = prepare(qb, rb, alphabet, matrix, go, ge, prior, width, threads,
                           weights_beta, query_groups_beta, reference_groups_beta, group_distances_beta, query_weights_beta, reference_weights_beta, feature_offsets_beta);
    std::vector<std::vector<uint64_t>> masses;
    auto counts = counted(a, &b, thresholds, go, ge, prior, threads, exclude_exact, &masses, true, shared_thresholds);
    return {std::move(counts), std::move(masses)};
}

struct GroupedCountBatch::Impl {
    GapInput a;
    std::optional<GapInput> b;
    std::vector<int32_t> cuts, prior, target_peptides, exclusive, qgenes, rgenes;
    std::vector<uint32_t> offsets, labels, peptide_offsets;
    std::vector<int32_t> peptide_ids;
    std::optional<std::vector<int32_t>> withheld;
    size_t targets = 0, peptides = 0, channels = 0, chunk_size = 0, cursor = 0;
    int32_t go = 0, ge = 0;
    int32_t quantum_a = 1, quantum_b = 1;
    int threads = 0;
    bool exclude = false;
    std::mutex mutex;
    std::optional<std::vector<int32_t>> beta_prior;
    std::vector<PairedFeatureOffsets> extra_a,extra_b;
    uint64_t extra_bound_a=0,extra_bound_b=0;
    std::optional<SubstitutionMatrix> similarity;
    AlignMode mode = AlignMode::Global;
    int32_t scale = 1;
    uint32_t left = 0, right = 0;
    std::vector<int64_t> qdiag_a, rdiag_a, qdiag_b, rdiag_b;
    int64_t bin_step = 1;
    std::vector<size_t> score_bins;

    void prepare_bins() {
        if (cuts.empty() || cuts.back()<0) return;
        int64_t step=0;
        for (const auto cut:cuts) step=std::gcd(step,cut<0 ? -int64_t(cut) : int64_t(cut));
        bin_step=std::max(int64_t(1),step);
        const uint64_t slots=uint64_t(cuts.back()/bin_step)+1, K=cuts.size();
        // Building more than K squared bins cannot repay even K binary searches.
        if (slots/K+uint64_t(slots%K!=0)>K) return;
        score_bins.resize(size_t(slots));
        size_t bin=size_t(std::lower_bound(cuts.begin(),cuts.end(),0)-cuts.begin());
        for (size_t i=0;i<score_bins.size();++i) {
            while (bin<cuts.size() && cuts[bin]<int64_t(i)*bin_step) ++bin;
            score_bins[i]=bin;
        }
    }

    size_t score_bin(int64_t score) const {
        if (!score_bins.empty() && score>=0) {
            if (score>cuts.back()) return cuts.size();
            return score_bins[size_t((score+bin_step-1)/bin_step)];
        }
        return size_t(std::lower_bound(cuts.begin(),cuts.end(),score)-cuts.begin());
    }

    std::vector<uint32_t> static_sizes;
    std::unordered_map<int32_t,std::vector<std::pair<uint32_t,uint32_t>>> gene_sizes;

    std::unordered_map<std::string,std::vector<uint32_t>> exact_references;

    std::string identity(size_t i, bool query) const {
        const auto& alpha=query ? a.q[i] : a.r[i];
        std::string key(alpha.begin(),alpha.end());
        if (b) {
            // Encoded symbols cannot equal the invalid-code separator.
            key.push_back(char(Codec::kInvalid));
            const auto& beta=query ? b->q[i] : b->r[i];
            key.append(beta.begin(),beta.end());
        }
        return key;
    }

    void prepare_sizes() {
        if (withheld) return;
        static_sizes.resize(channels);
        static_sizes[0]=uint32_t(a.r.size());
        // Sparse gene/target totals: at most one entry per distinct reference
        // membership, independent of the number of query rows.
        std::unordered_map<uint64_t,uint32_t> totals;
        for (size_t j=0;j<a.r.size();++j) {
            for (size_t p=offsets[j];p<offsets[j+1];++p) {
                const auto target=labels[p];
                ++static_sizes[1+target];
                ++totals[(uint64_t(uint32_t(rgenes[j]))<<32)|target];
            }
            exact_references[identity(j,false)].push_back(uint32_t(j));
            if (exclusive[j]>=0) ++static_sizes[1+2*targets+size_t(exclusive[j])];
        }
        for (const auto& [key,count]:totals)
            gene_sizes[int32_t(key>>32)].emplace_back(uint32_t(key),count);
    }

    void prepare_affine(const GroupedGeometry& geometry, const GapInput& lane,
            std::vector<int64_t>& qdiag, std::vector<int64_t>& rdiag,
            const std::vector<int32_t>& lane_prior,uint64_t extra_bound) {
        if (!geometry.weights.empty() || !geometry.query_weights.empty() ||
                !geometry.reference_weights.empty() || !lane_prior.empty())
            throw std::invalid_argument("affine alignment does not accept position weights or gap prior");
        uint64_t max_q=0, max_r=0;
        auto diagonals = [&](const auto& sequences, auto& output, uint64_t& longest) {
            output.reserve(sequences.size());
            for (const auto& sequence : sequences) {
                if (uint64_t(left)+right>sequence.size())
                    throw std::invalid_argument("alignment trims exceed sequence length");
                const size_t length=sequence.size()-left-right;
                longest=std::max(longest,uint64_t(length));
                int64_t diagonal=0;
                for (size_t k=left;k<sequence.size()-right;++k)
                    diagonal+=similarity->similarity(sequence[k],sequence[k]);
                output.push_back(diagonal);
            }
        };
        // Bound before summation/DP: |similarity| and gap magnitudes cover every
        // consumed residue; diagonal sum minus twice DP is <=3*(m+n)*max_cost.
        uint64_t maximum=std::max(go,ge);
        for (uint8_t x=0;x<lane.A;++x) for (uint8_t y=0;y<lane.A;++y) {
            const int64_t value=similarity->similarity(x,y);
            maximum=std::max(maximum,uint64_t(value<0 ? -value : value));
        }
        for (const auto& sequence:lane.q) {
            if (uint64_t(left)+right>sequence.size())
                throw std::invalid_argument("alignment trims exceed sequence length");
            max_q=std::max(max_q,uint64_t(sequence.size()-left-right));
        }
        for (const auto& sequence:lane.r) {
            if (uint64_t(left)+right>sequence.size())
                throw std::invalid_argument("alignment trims exceed sequence length");
            max_r=std::max(max_r,uint64_t(sequence.size()-left-right));
        }
        uint64_t additive=extra_bound;
        for (const auto& row:lane.group_distances)
            for (auto value:row) additive=std::max(additive,extra_bound+uint64_t(value));
        if (lane.features) {
            const auto& features=*lane.features;
            int32_t cost=0;
            for (const auto& row:features.costs) for (auto value:row) cost=std::max(cost,value);
            const size_t depth=(!features.query_codes.empty() ? features.query_codes.front().size() :
                (!features.reference_codes.empty() ? features.reference_codes.front().size() : 0))/2;
            if (depth && features.weight && features.scale)
                additive+=uint64_t(round_even(features.weight*double(uint64_t(cost)*depth)/double(depth))*features.scale);
        }
        const uint64_t limit=(uint64_t(INT32_MAX)-additive)/uint32_t(scale);
        if (maximum && max_q+max_r>limit/(3*maximum))
            throw std::overflow_error("affine grouped distance can overflow int32");
        diagonals(lane.q,qdiag,max_q); diagonals(lane.r,rdiag,max_r);
    }

    int64_t affine_cell(const GapInput& lane, size_t i, size_t j,
            const std::vector<int64_t>& qdiag, const std::vector<int64_t>& rdiag,
            detail::PairwiseScratch& scratch, int32_t offset) const {
        const auto& q=lane.q[i]; const auto& r=lane.r[j];
        // data() can be null for an empty vector; avoid pointer arithmetic then.
        const auto* qp=q.empty() ? q.data() : q.data()+left;
        const auto* rp=r.empty() ? r.data() : r.data()+left;
        const int64_t sim=detail::align_score_encoded(qp,q.size()-left-right,
            rp,r.size()-left-right,*similarity,mode,go,ge,scratch);
        return std::max(int64_t(0),qdiag[i]+rdiag[j]-2*sim)*scale+offset;
    }

    template<bool Additional> int64_t extra_offset(const std::vector<PairedFeatureOffsets>& terms,
            size_t i,size_t j) const {
        int64_t offset=0;
        if constexpr (Additional) for (const auto& term:terms) offset+=feature_offset(term,i,j);
        return offset;
    }

    template<bool Factored, bool Withheld=false, bool Quantized=false, bool Affine=false, bool Additional=false> GroupedCountChunk run(size_t start, size_t stop) {
        GroupedCountChunk out;
        out.start = start; out.stop = stop;
        const size_t rows = stop-start, K = cuts.size();
        auto make_buffer = [](std::vector<size_t> shape) {
            size_t cells = 1;
            for (size_t axis : shape) {
                if (axis && cells > size_t(PTRDIFF_MAX)/sizeof(uint32_t)/axis)
                    throw std::overflow_error("grouped count buffer exceeds addressable memory");
                cells *= axis;
            }
            return CountBuffer{std::vector<uint32_t>(cells),std::move(shape)};
        };
        out.counts = make_buffer({rows,channels,K});
        out.sizes = make_buffer({rows,channels});
        out.exact_counts = make_buffer({rows,targets});
        const auto longest = b ? std::max(a.longest,b->longest) : a.longest;
        const auto& other_prior=beta_prior ? *beta_prior : prior;
        parallel_for(rows,threads,[&] {
                if constexpr (Affine) {
                    detail::PairwiseScratch scratch;
                    for (auto& row:scratch) row.reserve(size_t(longest)+1);
                    return scratch;
                } else return std::vector<int32_t>(size_t(longest)+1);
            }, [&](size_t row, auto& scratch) {
                const size_t i = start+row;
                auto* sizes = out.sizes.data.data()+row*channels;
                auto* hist = K ? out.counts.data.data()+row*channels*K : nullptr;
                auto* exact = targets ? out.exact_counts.data.data()+row*targets : nullptr;
                const std::vector<uint32_t>* identical_refs=nullptr;
                size_t identical_cursor=0;
                if constexpr (!Withheld) {
                    std::copy(static_sizes.begin(),static_sizes.end(),sizes);
                    if (const auto found=gene_sizes.find(qgenes[i]);found!=gene_sizes.end())
                        for (const auto& [target,count]:found->second)
                            sizes[1+targets+target]=count;
                    if (const auto found=exact_references.find(identity(i,true));found!=exact_references.end())
                        identical_refs=&found->second;
                }
                for (size_t j=0;j<a.r.size();++j) {
                    int32_t exclusive_peptide=exclusive[j], withheld_peptide=-1;
                    if constexpr (Withheld) {
                        withheld_peptide=(*withheld)[i];
                        const auto begin=peptide_ids.begin()+peptide_offsets[j];
                        const auto end=peptide_ids.begin()+peptide_offsets[j+1];
                        const bool removed=std::binary_search(begin,end,withheld_peptide);
                        const size_t surviving=size_t(end-begin)-size_t(removed);
                        if (!surviving) continue;
                        exclusive_peptide=-1;
                        if (surviving==1) {
                            const int32_t peptide=*begin==withheld_peptide ? *(begin+1) : *begin;
                            if (size_t(peptide)<peptides) exclusive_peptide=peptide;
                        }
                    }
                    bool identical;
                    if constexpr (Withheld) identical=a.q[i]==a.r[j] && (!b || b->q[i]==b->r[j]);
                    else {
                        identical=identical_refs && identical_cursor<identical_refs->size() &&
                            (*identical_refs)[identical_cursor]==j;
                        if (identical) ++identical_cursor;
                    }
                    const bool same_gene = qgenes[i]==rgenes[j];
                    if (identical)
                        for (size_t p=offsets[j];p<offsets[j+1];++p) {
                            if constexpr (Withheld)
                                if (target_peptides[labels[p]]==withheld_peptide) continue;
                            ++exact[labels[p]];
                        }
                    if constexpr (!Withheld) {
                        if (exclude && identical) {
                            --sizes[0];
                            for (size_t p=offsets[j];p<offsets[j+1];++p) {
                                --sizes[1+labels[p]];
                                if (same_gene) --sizes[1+targets+labels[p]];
                            }
                            if (exclusive_peptide>=0) --sizes[1+2*targets+size_t(exclusive_peptide)];
                            continue;
                        }
                    } else {
                        if (exclude && identical) continue;
                        ++sizes[0];
                        for (size_t p=offsets[j];p<offsets[j+1];++p) {
                            if (target_peptides[labels[p]]==withheld_peptide) continue;
                            ++sizes[1+labels[p]];
                            if (same_gene) ++sizes[1+targets+labels[p]];
                        }
                        if (exclusive_peptide>=0) ++sizes[1+2*targets+size_t(exclusive_peptide)];
                    }
                    if (!K) continue;
                    const int64_t extra_alpha=extra_offset<Additional>(extra_a,i,j);
                    const int64_t extra_beta=b ? extra_offset<Additional>(extra_b,i,j) : 0;
                    const int32_t offset_alpha=group_offset<Factored>(a,i,j);
                    const int32_t offset_beta=b ? group_offset<Factored>(*b,i,j) : 0;
                    const int64_t offset=int64_t(offset_alpha)+extra_alpha+offset_beta+extra_beta;
                    if (offset>cuts.back()) continue;
                    int64_t score;
                    if constexpr (Affine) score=affine_cell(a,i,j,qdiag_a,rdiag_a,scratch,offset_alpha);
                    else score=score_cell(a,i,j,go,ge,prior,scratch.data(),offset_alpha);
                    score+=extra_alpha;
                    if constexpr (Quantized) {
                        const int64_t additive=offset_alpha+extra_alpha;
                        score=additive+(score-additive)/quantum_a*quantum_a;
                    }
                    if (b) {
                        if (score>cuts.back()) continue;
                        int64_t other;
                        if constexpr (Affine) other=affine_cell(*b,i,j,qdiag_b,rdiag_b,scratch,offset_beta);
                        else other=score_cell(*b,i,j,go,ge,other_prior,scratch.data(),offset_beta);
                        other+=extra_beta;
                        if constexpr (Quantized) {
                            const int64_t additive=offset_beta+extra_beta;
                            other=additive+(other-additive)/quantum_b*quantum_b;
                        }
                        score+=other;
                    }
                    const size_t bin=score_bin(score);
                    if (bin==K) continue;
                    ++hist[bin];
                    for (size_t p=offsets[j];p<offsets[j+1];++p) {
                        if constexpr (Withheld)
                            if (target_peptides[labels[p]]==withheld_peptide) continue;
                        ++hist[(1+labels[p])*K+bin];
                        if (same_gene) ++hist[(1+targets+labels[p])*K+bin];
                    }
                    if (exclusive_peptide>=0) ++hist[(1+2*targets+size_t(exclusive_peptide))*K+bin];
                }
                for (size_t channel=0;channel<channels;++channel)
                    for (size_t t=1;t<K;++t) hist[channel*K+t]+=hist[channel*K+t-1];
            });
        return out;
    }
    template<bool Affine,bool Additional> GroupedCountChunk dispatch(size_t start,size_t stop) {
        const bool features=a.features || (b && b->features);
        const bool quantized=quantum_a!=1 || quantum_b!=1;
        return quantized ?
            (withheld ?
                (features ? run<true,true,true,Affine,Additional>(start,stop) : run<false,true,true,Affine,Additional>(start,stop)) :
                (features ? run<true,false,true,Affine,Additional>(start,stop) : run<false,false,true,Affine,Additional>(start,stop))) :
            (withheld ?
                (features ? run<true,true,false,Affine,Additional>(start,stop) : run<false,true,false,Affine,Additional>(start,stop)) :
                (features ? run<true,false,false,Affine,Additional>(start,stop) : run<false,false,false,Affine,Additional>(start,stop)));
    }
};

GroupedCountBatch::GroupedCountBatch(const std::vector<std::string>& qa, const std::vector<std::string>& ra,
        std::vector<int32_t> cuts, std::vector<uint32_t> offsets, std::vector<uint32_t> labels,
        std::vector<int32_t> target_peptides, std::vector<int32_t> exclusive,
        std::vector<int32_t> qgenes, std::vector<int32_t> rgenes, size_t chunk_size,
        Alphabet alphabet, const SubstitutionMatrix* matrix, int32_t go, int32_t ge,
        std::vector<int32_t> prior, uint32_t width, int threads, bool exclude,
        const GroupedGeometry& alpha, const std::optional<GroupedGeometry>& beta,
        const std::vector<std::string>& qb, const std::vector<std::string>& rb,
        std::optional<std::vector<int32_t>> withheld,
        std::vector<uint32_t> peptide_offsets, std::vector<int32_t> peptide_ids)
        : impl_(std::make_shared<Impl>()) {
    auto& in=*impl_;
    if (alpha.score_quantum<1 || (beta && beta->score_quantum<1))
        throw std::invalid_argument("score quantum must be positive");
    in.quantum_a=alpha.score_quantum;
    in.quantum_b=beta ? beta->score_quantum : 1;
    if (!chunk_size) throw std::invalid_argument("chunk_size must be positive");
    if (!std::is_sorted(cuts.begin(),cuts.end()) || std::adjacent_find(cuts.begin(),cuts.end())!=cuts.end())
        throw std::invalid_argument("grouped thresholds must be strictly increasing");
    if (ra.size()>std::numeric_limits<uint32_t>::max())
        throw std::overflow_error("grouped reference count exceeds uint32");
    if (offsets.size()!=ra.size()+1 || offsets.front()!=0 || offsets.back()!=labels.size() ||
            !std::is_sorted(offsets.begin(),offsets.end()))
        throw std::invalid_argument("reference target offsets must form valid CSR");
    if (exclusive.size()!=ra.size() || qgenes.size()!=qa.size() || rgenes.size()!=ra.size())
        throw std::invalid_argument("grouped metadata axes must match queries/references");
    if (std::any_of(qgenes.begin(),qgenes.end(),[](auto x){return x<0;}) ||
        std::any_of(rgenes.begin(),rgenes.end(),[](auto x){return x<0;}))
        throw std::invalid_argument("gene IDs must be nonnegative");
    for (auto peptide:target_peptides) {
        if (peptide<0) throw std::invalid_argument("target peptide IDs must be nonnegative");
        in.peptides=std::max(in.peptides,size_t(peptide)+1);
    }
    in.targets=target_peptides.size();
    for (size_t j=0;j<ra.size();++j) {
        if (exclusive[j]<-1 || (exclusive[j]>=0 && size_t(exclusive[j])>=in.peptides))
            throw std::invalid_argument("exclusive peptide ID outside peptide axis");
        // Empty memberships and -1 exclusive IDs are valid: an included target
        // can share a reference with an unselected label in the whole population.
        uint32_t previous=0;
        for (size_t p=offsets[j];p<offsets[j+1];++p) {
            const auto label=labels[p];
            if (label>=in.targets || (p>offsets[j] && label<=previous))
                throw std::invalid_argument("reference target IDs must be sorted, unique and in range");
            if (exclusive[j]>=0 && target_peptides[label]!=exclusive[j])
                throw std::invalid_argument("exclusive peptide contradicts reference target membership");
            previous=label;
        }
    }
    if (withheld) {
        if (withheld->size()!=qa.size() ||
            std::any_of(withheld->begin(),withheld->end(),[](auto x){return x<-1;}))
            throw std::invalid_argument("withheld peptide IDs must match queries and be >= -1");
        if (peptide_ids.size()>std::numeric_limits<uint32_t>::max())
            throw std::overflow_error("reference peptide memberships exceed uint32");
        if (peptide_offsets.size()!=ra.size()+1 || peptide_offsets.front()!=0 ||
            peptide_offsets.back()!=peptide_ids.size() ||
            !std::is_sorted(peptide_offsets.begin(),peptide_offsets.end()))
            throw std::invalid_argument("reference peptide offsets must form valid CSR");
        for (size_t j=0;j<ra.size();++j) {
            const auto begin=peptide_ids.begin()+peptide_offsets[j];
            const auto end=peptide_ids.begin()+peptide_offsets[j+1];
            if (std::any_of(begin,end,[](auto x){return x<0;}) ||
                !std::is_sorted(begin,end) || std::adjacent_find(begin,end)!=end)
                throw std::invalid_argument("reference peptide IDs must be sorted, unique and nonnegative");
            for (size_t p=offsets[j];p<offsets[j+1];++p)
                if (!std::binary_search(begin,end,target_peptides[labels[p]]))
                    throw std::invalid_argument("reference peptide memberships omit selected target parent");
            const int32_t expected=end-begin==1 && size_t(*begin)<in.peptides ? *begin : -1;
            if (exclusive[j]!=expected)
                throw std::invalid_argument("exclusive peptide contradicts full reference peptide membership");
        }
    } else if (!peptide_offsets.empty() || !peptide_ids.empty())
        throw std::invalid_argument("reference peptide memberships require query withheld peptide IDs");
    in.withheld=std::move(withheld);
    in.peptide_offsets=std::move(peptide_offsets);in.peptide_ids=std::move(peptide_ids);
    const auto& alpha_prior=alpha.gap_prior ? *alpha.gap_prior : prior;
    const auto& beta_prior=beta && beta->gap_prior ? *beta->gap_prior : prior;
    uint64_t bound_a=0,bound_b=0;
    in.a=prepare(qa,ra,alphabet,matrix,go,ge,alpha_prior,width,threads,alpha.weights,
        alpha.query_groups,alpha.reference_groups,alpha.group_distances,
        alpha.query_weights,alpha.reference_weights,alpha.features,&bound_a);
    if (beta) {
        if (qa.size()!=qb.size() || ra.size()!=rb.size())
            throw std::invalid_argument("paired query/reference axes must agree");
        in.b=prepare(qb,rb,alphabet,matrix,go,ge,beta_prior,width,threads,beta->weights,
            beta->query_groups,beta->reference_groups,beta->group_distances,
            beta->query_weights,beta->reference_weights,beta->features,&bound_b);
    } else if (!qb.empty() || !rb.empty()) throw std::invalid_argument("paired axes require beta geometry");
    if (threads<0 || unsigned(threads)>available_cpus())
        throw std::invalid_argument("worker count exceeds the available CPU budget");
    in.cuts=std::move(cuts); in.offsets=std::move(offsets); in.labels=std::move(labels);
    in.prepare_bins();
    in.target_peptides=std::move(target_peptides); in.exclusive=std::move(exclusive);
    in.qgenes=std::move(qgenes);in.rgenes=std::move(rgenes);
    in.prior=alpha.gap_prior ? *alpha.gap_prior : std::move(prior);
    if (beta && (beta->gap_prior || alpha.gap_prior))
        in.beta_prior=beta->gap_prior ? *beta->gap_prior : prior;
    auto additional=[&](const auto& terms,size_t nq,size_t nr,uint64_t base_bound) {
        uint64_t bound=0;
        for (const auto& term:terms) {
            const uint64_t next=validate_feature(term,nq,nr);
            if (next>uint64_t(INT32_MAX)-base_bound-bound)
                throw std::overflow_error("additional factored offsets can overflow int32");
            bound+=next;
        }
        return bound;
    };
    in.extra_bound_a=additional(alpha.additional_features,qa.size(),ra.size(),bound_a);
    in.extra_a=alpha.additional_features;
    if (beta) {
        in.extra_bound_b=additional(beta->additional_features,qb.size(),rb.size(),bound_b);
        in.extra_b=beta->additional_features;
    }
    in.channels=1+2*in.targets+in.peptides;in.chunk_size=chunk_size;
    in.go=go;in.ge=ge;in.threads=threads;in.exclude=exclude;
    if (!in.a.q.empty()) in.prepare_sizes();
    if (alpha.alignment_scale<1 || (alpha.alignment_mode!="gapblock" &&
            alpha.alignment_mode!="global" && alpha.alignment_mode!="local"))
        throw std::invalid_argument("invalid grouped alignment mode or scale");
    if (beta && (beta->alignment_mode!=alpha.alignment_mode ||
            beta->alignment_scale!=alpha.alignment_scale || beta->trim_left!=alpha.trim_left ||
            beta->trim_right!=alpha.trim_right))
        throw std::invalid_argument("paired alignment options must agree");
    if (alpha.alignment_mode=="gapblock") {
        if (alpha.alignment_scale!=1 || alpha.trim_left || alpha.trim_right)
            throw std::invalid_argument("alignment scale/trims require global or local mode");
    } else {
        if (!matrix) throw std::invalid_argument("affine alignment requires a similarity matrix");
        in.similarity=*matrix;
        in.mode=alpha.alignment_mode=="local" ? AlignMode::Local : AlignMode::Global;
        in.scale=alpha.alignment_scale;in.left=alpha.trim_left;in.right=alpha.trim_right;
        in.prepare_affine(alpha,in.a,in.qdiag_a,in.rdiag_a,in.prior,in.extra_bound_a);
        if (beta) in.prepare_affine(*beta,*in.b,in.qdiag_b,in.rdiag_b,in.beta_prior ? *in.beta_prior : in.prior,in.extra_bound_b);
    }
}
std::optional<GroupedCountChunk> GroupedCountBatch::next() {
    auto& in=*impl_;
    std::lock_guard<std::mutex> lock(in.mutex);
    if (in.cursor==in.a.q.size()) return std::nullopt;
    const size_t stop=in.cursor+std::min(in.chunk_size,in.a.q.size()-in.cursor);
    const bool extra=!in.extra_a.empty() || !in.extra_b.empty();
    auto out=in.similarity ?
        (extra ? in.dispatch<true,true>(in.cursor,stop) : in.dispatch<true,false>(in.cursor,stop)) :
        (extra ? in.dispatch<false,true>(in.cursor,stop) : in.dispatch<false,false>(in.cursor,stop));
    in.cursor=stop;
    return out;
}
size_t GroupedCountBatch::query_count() const {return impl_->a.q.size();}
size_t GroupedCountBatch::reference_count() const {return impl_->a.r.size();}
size_t GroupedCountBatch::target_count() const {return impl_->targets;}
size_t GroupedCountBatch::peptide_count() const {return impl_->peptides;}

}  // namespace seqtree
