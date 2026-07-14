// Needleman-Wunsch (global) and Smith-Waterman (local) with affine gaps, on the raw
// similarity scale. The conventions match Bio.Align.PairwiseAligner exactly, and
// tests/python/test_pairwise.py pins every one of them against BioPython as an oracle.
//
// Gotoh's three-state recurrence, MAXIMISING (everything else in this library minimises a
// non-negative penalty; this does not):
//
//   M[i][j] = s(q_i, r_j) + max(M, X, Y)[i-1][j-1]     both residues aligned
//   X[i][j] = max(M[i-1][j] - open, X[i-1][j] - extend, Y[i-1][j] - open)   query residue vs gap
//   Y[i][j] = max(M[i][j-1] - open, Y[i][j-1] - extend, X[i][j-1] - open)   ref residue vs gap
//
// so a gap run of length L costs open + (L-1)*extend, and open == extend gives linear gaps.
// Local adds a 0 floor on M and takes the best M cell anywhere (Smith-Waterman); the gap
// states are NOT floored -- a local alignment never ends in a gap, you would just trim it.

#include "seqtree/seqtree.hpp"

#include <algorithm>
#include <atomic>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace seqtree {
namespace {

// Low enough to be unreachable, with headroom so that kNegInf - gap does not overflow.
constexpr int32_t kNegInf = std::numeric_limits<int32_t>::min() / 4;

std::vector<uint8_t> encode(const Codec& c, std::string_view s) {
    std::vector<uint8_t> out(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        uint8_t v = c.encode(s[i]);
        if (v == Codec::kInvalid)
            throw std::invalid_argument(std::string("invalid symbol '") + s[i] + "' in sequence");
        out[i] = v;
    }
    return out;
}

void check(const SubstitutionMatrix& mat, const Codec& codec, int32_t open, int32_t extend) {
    if (open < 0 || extend < 0)
        throw std::invalid_argument("seqtree: gap_open and gap_extend must be >= 0 (they are "
                                    "magnitudes, and are subtracted)");
    if (mat.size() != codec.size())
        throw std::invalid_argument("matrix size does not match the alphabet");
}

// Cost of a terminal gap run of length L in global mode.
inline int32_t gap_run(int32_t L, int32_t open, int32_t extend) {
    return L <= 0 ? 0 : open + (L - 1) * extend;
}

// Optimal score, two rows per state. O(n) memory, no traceback.
int32_t score_encoded(const uint8_t* q, size_t m, const uint8_t* r, size_t n,
                      const SubstitutionMatrix& mat, AlignMode mode, int32_t go, int32_t ge) {
    const bool local = (mode == AlignMode::Local);
    if (m == 0 || n == 0) {
        if (local) return 0;
        return -gap_run(int32_t(m ? m : n), go, ge);
    }

    std::vector<int32_t> Mp(n + 1), Xp(n + 1), Yp(n + 1);
    std::vector<int32_t> Mc(n + 1), Xc(n + 1), Yc(n + 1);

    // Row 0. Global: the only way to consume ref[0..j) is a leading gap (state Y). Local: an
    // alignment may start anywhere, which the 0 floor on M expresses; the gap states cannot be
    // entered before the alignment has started, so they stay unreachable.
    Mp[0] = 0;
    Xp[0] = kNegInf;
    Yp[0] = kNegInf;
    for (size_t j = 1; j <= n; ++j) {
        Mp[j] = local ? 0 : kNegInf;
        Xp[j] = kNegInf;
        Yp[j] = local ? kNegInf : -gap_run(int32_t(j), go, ge);
    }

    int32_t best = local ? 0 : kNegInf;

    for (size_t i = 1; i <= m; ++i) {
        Mc[0] = local ? 0 : kNegInf;
        Xc[0] = local ? kNegInf : -gap_run(int32_t(i), go, ge);
        Yc[0] = kNegInf;

        for (size_t j = 1; j <= n; ++j) {
            const int32_t s = mat.similarity(q[i - 1], r[j - 1]);
            const int32_t diag = std::max(Mp[j - 1], std::max(Xp[j - 1], Yp[j - 1]));
            int32_t mv = diag + s;
            if (local && mv < 0) mv = 0;
            Mc[j] = mv;

            Xc[j] = std::max(Mp[j] - go, std::max(Xp[j] - ge, Yp[j] - go));
            Yc[j] = std::max(Mc[j - 1] - go, std::max(Yc[j - 1] - ge, Xc[j - 1] - go));

            if (local && mv > best) best = mv;
        }
        Mp.swap(Mc);
        Xp.swap(Xc);
        Yp.swap(Yc);
    }

    if (local) return best;
    return std::max(Mp[n], std::max(Xp[n], Yp[n]));
}

enum : uint8_t { kM = 0, kX = 1, kY = 2, kStop = 3 };

}  // namespace

int32_t align_score(std::string_view query, std::string_view ref, const SubstitutionMatrix& mat,
                    Alphabet alphabet, AlignMode mode, int32_t gap_open, int32_t gap_extend) {
    const Codec codec(alphabet);
    check(mat, codec, gap_open, gap_extend);
    const auto q = encode(codec, query);
    const auto r = encode(codec, ref);
    return score_encoded(q.data(), q.size(), r.data(), r.size(), mat, mode, gap_open, gap_extend);
}

Alignment align_pair(std::string_view query, std::string_view ref, const SubstitutionMatrix& mat,
                     Alphabet alphabet, AlignMode mode, int32_t gap_open, int32_t gap_extend) {
    const Codec codec(alphabet);
    check(mat, codec, gap_open, gap_extend);
    const auto q = encode(codec, query);
    const auto r = encode(codec, ref);
    const int m = int(q.size()), n = int(r.size());
    const bool local = (mode == AlignMode::Local);
    const int32_t go = gap_open, ge = gap_extend;

    Alignment al;
    if (m == 0 || n == 0) {
        al.score = local ? 0 : -gap_run(m ? m : n, go, ge);
        if (!local) {
            for (int i = 0; i < m; ++i) {
                al.aligned_query += query[size_t(i)];
                al.aligned_ref += '-';
                al.ops += 'D';
            }
            for (int j = 0; j < n; ++j) {
                al.aligned_query += '-';
                al.aligned_ref += ref[size_t(j)];
                al.ops += 'I';
            }
        }
        return al;
    }

    const int W = n + 1;
    const size_t sz = size_t(m + 1) * W;
    std::vector<int32_t> M(sz, kNegInf), X(sz, kNegInf), Y(sz, kNegInf);
    std::vector<uint8_t> bM(sz, kStop), bX(sz, kM), bY(sz, kM);

    M[0] = 0;
    for (int j = 1; j <= n; ++j) {
        M[j] = local ? 0 : kNegInf;
        Y[j] = local ? kNegInf : -gap_run(j, go, ge);
        bY[j] = (j > 1) ? kY : kM;
    }
    for (int i = 1; i <= m; ++i) {
        const size_t k = size_t(i) * W;
        M[k] = local ? 0 : kNegInf;
        X[k] = local ? kNegInf : -gap_run(i, go, ge);
        bX[k] = (i > 1) ? kX : kM;
    }

    // argmax over the three states, reporting which one won.
    auto pick = [](int32_t vm, int32_t vx, int32_t vy, uint8_t& who) {
        int32_t v = vm;
        who = kM;
        if (vx > v) { v = vx; who = kX; }
        if (vy > v) { v = vy; who = kY; }
        return v;
    };

    int32_t best = local ? 0 : kNegInf;
    int bi = 0, bj = 0;

    for (int i = 1; i <= m; ++i) {
        for (int j = 1; j <= n; ++j) {
            const size_t k = size_t(i) * W + j;
            const size_t kd = size_t(i - 1) * W + (j - 1);
            const size_t ku = size_t(i - 1) * W + j;
            const size_t kl = size_t(i) * W + (j - 1);

            uint8_t who;
            int32_t mv = pick(M[kd], X[kd], Y[kd], who) + mat.similarity(q[i - 1], r[j - 1]);
            if (local && mv < 0) { mv = 0; who = kStop; }
            M[k] = mv;
            bM[k] = who;

            X[k] = pick(M[ku] - go, X[ku] - ge, Y[ku] - go, bX[k]);
            Y[k] = pick(M[kl] - go, X[kl] - go, Y[kl] - ge, bY[k]);

            if (local && mv > best) { best = mv; bi = i; bj = j; }
        }
    }

    int i, j;
    uint8_t st;
    if (local) {
        al.score = best;
        i = bi;
        j = bj;
        st = kM;
        if (best == 0) return al;  // no positive-scoring local alignment exists
    } else {
        const size_t end = size_t(m) * W + n;
        al.score = pick(M[end], X[end], Y[end], st);
        i = m;
        j = n;
    }

    while (i > 0 || j > 0) {
        const size_t k = size_t(i) * W + j;
        if (st == kM) {
            if (i == 0 || j == 0) break;
            if (local && M[k] == 0) break;  // Smith-Waterman stops where the score resets
            const char a = codec.decode(q[i - 1]), b = codec.decode(r[j - 1]);
            al.aligned_query += a;
            al.aligned_ref += b;
            al.ops += (a == b ? 'M' : 'S');
            st = bM[k];
            --i;
            --j;
            if (st == kStop) break;
        } else if (st == kX) {
            if (i == 0) break;
            al.aligned_query += codec.decode(q[i - 1]);
            al.aligned_ref += '-';
            al.ops += 'D';
            st = bX[k];
            --i;
        } else {
            if (j == 0) break;
            al.aligned_query += '-';
            al.aligned_ref += codec.decode(r[j - 1]);
            al.ops += 'I';
            st = bY[k];
            --j;
        }
    }
    std::reverse(al.aligned_query.begin(), al.aligned_query.end());
    std::reverse(al.aligned_ref.begin(), al.aligned_ref.end());
    std::reverse(al.ops.begin(), al.ops.end());
    return al;
}

namespace {

// Shared driver for the two N*K matrices. `dist` applies the sequence-level Gram transform
// d = s(a,a) + s(b,b) - 2*s(a,b); the self-scores are computed once per sequence, not per pair.
std::vector<int32_t> matrix_impl(const std::vector<std::string>& queries,
                                 const std::vector<std::string>& refs,
                                 const SubstitutionMatrix& mat, Alphabet alphabet, AlignMode mode,
                                 int32_t go, int32_t ge, int threads, bool dist) {
    const Codec codec(alphabet);
    check(mat, codec, go, ge);

    const size_t N = queries.size(), K = refs.size();
    std::vector<int32_t> out(N * K);
    if (N == 0 || K == 0) return out;

    std::vector<std::vector<uint8_t>> qc(N), rc(K);
    for (size_t i = 0; i < N; ++i) qc[i] = encode(codec, queries[i]);
    for (size_t k = 0; k < K; ++k) rc[k] = encode(codec, refs[k]);

    std::vector<int32_t> q_self(N), r_self(K);
    if (dist) {
        for (size_t i = 0; i < N; ++i)
            q_self[i] = score_encoded(qc[i].data(), qc[i].size(), qc[i].data(), qc[i].size(), mat,
                                      mode, go, ge);
        for (size_t k = 0; k < K; ++k)
            r_self[k] = score_encoded(rc[k].data(), rc[k].size(), rc[k].data(), rc[k].size(), mat,
                                      mode, go, ge);
    }

    unsigned nt =
        threads > 0 ? unsigned(threads) : std::max(1u, std::thread::hardware_concurrency());
    nt = std::min<unsigned>(nt, std::max<size_t>(1, N));

    std::atomic<size_t> next{0};
    auto worker = [&] {
        for (;;) {
            const size_t i = next.fetch_add(1);
            if (i >= N) break;
            int32_t* row = out.data() + i * K;
            for (size_t k = 0; k < K; ++k) {
                const int32_t s = score_encoded(qc[i].data(), qc[i].size(), rc[k].data(),
                                                rc[k].size(), mat, mode, go, ge);
                row[k] = dist ? (q_self[i] + r_self[k] - 2 * s) : s;
            }
        }
    };

    std::vector<std::thread> pool;
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    return out;
}

}  // namespace

std::vector<int32_t> align_score_matrix(const std::vector<std::string>& queries,
                                        const std::vector<std::string>& refs,
                                        const SubstitutionMatrix& mat, Alphabet alphabet,
                                        AlignMode mode, int32_t gap_open, int32_t gap_extend,
                                        int threads) {
    return matrix_impl(queries, refs, mat, alphabet, mode, gap_open, gap_extend, threads, false);
}

std::vector<int32_t> align_dist_matrix(const std::vector<std::string>& queries,
                                       const std::vector<std::string>& refs,
                                       const SubstitutionMatrix& mat, Alphabet alphabet,
                                       AlignMode mode, int32_t gap_open, int32_t gap_extend,
                                       int threads) {
    return matrix_impl(queries, refs, mat, alphabet, mode, gap_open, gap_extend, threads, true);
}

}  // namespace seqtree
