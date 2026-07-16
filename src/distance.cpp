// Plain string edit distances -- Hamming and Levenshtein -- on raw characters, with unit
// costs. Deliberately alphabet-agnostic: no Codec, no substitution matrix, no gap model. This
// is the distance you want when you just need a number, not a scored alignment. Everything
// else in this library scores against a matrix; these two do not, so they live apart from
// pairwise.cpp even though the threading plumbing is the same.
//
// Comparison is byte-for-byte and case-sensitive -- it does NOT fold case the way the search
// Codec does, because these are generic string distances, not sequence-search primitives.

#include "seqtree/seqtree.hpp"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace seqtree {
namespace {

// Shared driver for both N*K matrices. `pair` is the scalar distance applied to every cell;
// it may throw (Hamming does, on a length mismatch), so the worker catches and the exception
// is rethrown after the join -- an exception escaping a std::thread entry calls std::terminate
// and would take the interpreter down with an uncatchable SIGABRT. Same plumbing as
// pairwise.cpp / Index::search_batch.
template <class Pair>
std::vector<int32_t> matrix_impl(const std::vector<std::string>& a,
                                 const std::vector<std::string>& b, int threads, Pair pair) {
    const size_t N = a.size(), K = b.size();
    std::vector<int32_t> out(N * K);
    if (N == 0 || K == 0) return out;

    unsigned nt =
        threads > 0 ? unsigned(threads) : std::max(1u, std::thread::hardware_concurrency());
    nt = std::min<unsigned>(nt, std::max<size_t>(1, N));

    std::atomic<size_t> next{0};
    std::exception_ptr err;
    std::mutex emu;

    auto worker = [&] {
        for (;;) {
            const size_t i = next.fetch_add(1);
            if (i >= N) break;
            try {
                int32_t* row = out.data() + i * K;
                for (size_t k = 0; k < K; ++k) row[k] = pair(a[i], b[k]);
            } catch (...) {
                std::lock_guard<std::mutex> lk(emu);
                if (!err) err = std::current_exception();
                return;
            }
        }
    };

    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    if (err) std::rethrow_exception(err);
    return out;
}

}  // namespace

int hamming(std::string_view a, std::string_view b) {
    if (a.size() != b.size())
        throw std::invalid_argument("seqtree: hamming distance is only defined for equal-length "
                                    "sequences (got " + std::to_string(a.size()) + " and " +
                                    std::to_string(b.size()) + ")");
    int d = 0;
    for (size_t i = 0; i < a.size(); ++i) d += (a[i] != b[i]);
    return d;
}

int levenshtein(std::string_view a, std::string_view b) {
    // Keep the inner row over the shorter sequence: O(min(m,n)) memory. Two rows, no traceback.
    if (b.size() > a.size()) std::swap(a, b);
    const size_t m = a.size(), n = b.size();
    if (n == 0) return int(m);

    std::vector<int> prev(n + 1), cur(n + 1);
    for (size_t j = 0; j <= n; ++j) prev[j] = int(j);
    for (size_t i = 1; i <= m; ++i) {
        cur[0] = int(i);
        for (size_t j = 1; j <= n; ++j) {
            const int cost = (a[i - 1] == b[j - 1]) ? 0 : 1;
            cur[j] = std::min({prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost});
        }
        prev.swap(cur);
    }
    return prev[n];
}

std::vector<int32_t> hamming_matrix(const std::vector<std::string>& a,
                                    const std::vector<std::string>& b, int threads) {
    return matrix_impl(a, b, threads, [](std::string_view x, std::string_view y) {
        return hamming(x, y);
    });
}

std::vector<int32_t> levenshtein_matrix(const std::vector<std::string>& a,
                                        const std::vector<std::string>& b, int threads) {
    return matrix_impl(a, b, threads, [](std::string_view x, std::string_view y) {
        return levenshtein(x, y);
    });
}

}  // namespace seqtree
