#include "seqtree/seqtree.hpp"

#include <algorithm>
#include <atomic>
#include <limits>
#include <stdexcept>
#include <thread>

namespace seqtree {
namespace {

std::vector<uint8_t> encode(const Codec& c, const std::string& s) {
    std::vector<uint8_t> out(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        uint8_t v = c.encode(s[i]);
        if (v == Codec::kInvalid)
            throw std::invalid_argument(std::string("symbol '") + s[i] +
                                        "' is not in the alphabet");
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
                    const int32_t* prior_row, int32_t* suf) {
    const uint32_t L = std::min(m, n);
    const uint32_t d = (m > n ? m - n : n - m);
    const bool q_longer = m >= n;

    suf[L] = 0;
    for (uint32_t j = L; j-- > 0;) {
        const uint8_t a = q_longer ? q[j + d] : q[j];
        const uint8_t b = q_longer ? r[j] : r[j + d];
        suf[j] = suf[j + 1] + pen[size_t(a) * A + b];
    }
    if (d == 0) return suf[0];

    int32_t pre = 0;
    int32_t best = std::numeric_limits<int32_t>::max();
    for (uint32_t i = 0; i <= L; ++i) {
        const int32_t cand = pre + suf[i] + (prior_row ? prior_row[i] : 0);
        if (cand < best) best = cand;
        if (i < L) pre += pen[size_t(q[i]) * A + r[i]];
    }
    return best + gap_open + int32_t(d - 1) * gap_extend;
}

}  // namespace

std::vector<int32_t> gapblock_matrix(const std::vector<std::string>& queries,
                                     const std::vector<std::string>& refs, Alphabet alphabet,
                                     const SubstitutionMatrix* matrix, int32_t gap_open,
                                     int32_t gap_extend, const std::vector<int32_t>& prior,
                                     uint32_t prior_width, int threads) {
    if (gap_open < 0 || gap_extend < 0)
        throw std::invalid_argument("gap_open and gap_extend must be >= 0");

    const Codec codec(alphabet);
    const uint8_t A = codec.size();
    if (matrix && matrix->size() != A)
        throw std::invalid_argument("matrix size does not match the alphabet");

    const size_t N = queries.size(), K = refs.size();
    std::vector<int32_t> out(N * K);
    if (N == 0 || K == 0) return out;

    const std::vector<int32_t> pen = pen_table(matrix, A);

    std::vector<std::vector<uint8_t>> qc(N), rc(K);
    uint32_t longest = 0;
    for (size_t i = 0; i < N; ++i) {
        qc[i] = encode(codec, queries[i]);
        longest = std::max<uint32_t>(longest, uint32_t(qc[i].size()));
    }
    for (size_t k = 0; k < K; ++k) {
        rc[k] = encode(codec, refs[k]);
        longest = std::max<uint32_t>(longest, uint32_t(rc[k].size()));
    }

    const size_t W1 = size_t(prior_width) + 1;
    if (!prior.empty()) {
        if (prior.size() != W1 * W1 * W1)
            throw std::invalid_argument("prior table must have (prior_width + 1)^3 entries");
        if (longest > prior_width)
            throw std::invalid_argument("a sequence is longer than the prior table's width");
    }
    const int32_t* P = prior.empty() ? nullptr : prior.data();

    unsigned nt = threads > 0 ? unsigned(threads) : std::max(1u, std::thread::hardware_concurrency());
    nt = std::min<unsigned>(nt, std::max<size_t>(1, N));

    // Every symbol was validated above and the kernel is arithmetic, so no worker throws
    // and rows are disjoint -- no exception plumbing, no locking.
    std::atomic<size_t> next{0};
    auto worker = [&] {
        std::vector<int32_t> suf(size_t(longest) + 1);
        for (;;) {
            const size_t i = next.fetch_add(1);
            if (i >= N) break;
            const uint8_t* q = qc[i].data();
            const uint32_t m = uint32_t(qc[i].size());
            int32_t* row = out.data() + i * K;
            for (size_t k = 0; k < K; ++k) {
                const uint32_t n = uint32_t(rc[k].size());
                const uint32_t M = std::max(m, n), d = (m > n ? m - n : n - m);
                const int32_t* prow = P ? P + (size_t(M) * W1 + d) * W1 : nullptr;
                row[k] = cell(q, m, rc[k].data(), n, pen.data(), A, gap_open, gap_extend, prow,
                              suf.data());
            }
        }
    };

    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    return out;
}

}  // namespace seqtree
