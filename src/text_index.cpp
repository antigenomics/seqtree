// Exact k-mismatch (Hamming) search over a concatenated reference text.
//
// The structure is one direct-addressed k-mer seed table over a flat text, so `k` belongs to
// the INDEX and not to the query: one build answers every query length and every max_subs.
//
// ONE SEARCH SCHEME, not two paths. Split the query into `b` disjoint blocks of width
// bw = L / b >= k, give block j an error budget c_j, and probe the radius-c_j neighbourhood of
// that block's leading k-mer.
//
//   LOSSLESS  <=>  sum_j c_j >= m - b + 1.
//
//   Proof. The scheme misses an occurrence only if its error vector (e_0 .. e_{b-1}), with
//   sum e_j <= m, has e_j > c_j for EVERY j. The cheapest such vector is e_j = c_j + 1, whose
//   sum is sum(c_j) + b. So no missing vector exists iff sum(c_j) + b > m. Errors landing past
//   b*bw belong to no block and only shrink sum e_j, so an uncovered tail is free.
//
// Probing block j at radius c costs N(c) = sum_{i<=c} C(k,i)(A-1)^i, whose increments
// C(k,c)(A-1)^c rise steeply in c and are the same for every block -- so the cheapest legal
// scheme takes b as large as it can and spreads the budget as evenly as possible:
//
//   b = min(m + 1, L / k);  r = max(0, m - b + 1);  c_j = r/b + (j < r%b)
//
// Both of the two paths this replaced are special cases: b = m+1 gives r = 0 and all-exact
// seeds (pigeonhole), b = 1 gives c_0 = m and a single ball over q[0:k]. What was missing was
// everything between, which is where the short queries live: at L = 8, m = 2, k = 4 two
// disjoint 4-mers plainly fit, and the scheme (1,0) probes 94 variants where a one-block ball
// probes 3,267 -- 35x fewer candidates on the human proteome, and 186-279x at m = 3.
// The framing is the "search scheme" of Kianfar, Pockrandt, Torkamandi, Luo & Reinert,
// arXiv:1711.02035; only the partition-and-budget half transfers, since a direct-addressed
// table has no bidirectional extension to optimise a search ORDER over.
//
// Under proper substitutions (replacement != the original residue) each block's enumeration is
// duplicate-free by construction, so it is a direct write into a reused buffer -- no sort, no
// hash set, no allocation per query. Across blocks a start CAN repeat, and it is emitted by the
// lowest-indexed block that could have produced it; that test is read off the mismatch positions
// verification already computed, so deduplication costs no extra memory traffic and no sort.
//
// Codec::kInvalid between records is load-bearing: a k-mer containing it is never inserted, so
// a seed cannot span a record boundary, and verification refuses it outright, so a hit cannot
// either. Refusing rather than counting it as a mismatch is the point -- counting would let a
// query with max_subs >= 2 straddle two records by "paying" for the sentinel.

#include "seqtree/text_index.hpp"
#include "seqtree/parallel.hpp"

#include "atomic_write.hpp"

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <stdexcept>

#ifdef _WIN32
// NOMINMAX or windows.h defines min/max as macros and breaks every std::min / std::max /
// std::clamp below; WIN32_LEAN_AND_MEAN keeps the rest of the header out.
#define NOMINMAX
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace seqtree {
namespace {

constexpr char kMagic[4] = {'S', 'Q', 'T', 'X'};
constexpr uint32_t kVersion = 2;
constexpr size_t kHeaderBytes = 64;

// Direct addressing means alphabet_size^k uint32 buckets are allocated up front, so k has to be
// capped or a k=8 amino-acid index would ask for 44 GB. 2^27 buckets is 512 MB, already far more
// than any measured configuration needs (24^4 = 331,776; 24^5 = 7,962,624).
constexpr uint64_t kMaxBuckets = uint64_t(1) << 27;

// TextResult::mm_pos is uint16_t, so a longer query would wrap a mismatch position silently.
constexpr size_t kMaxQueryLen = 65535;

// The first read-only memory map in the tree. Only load(mmap=true) uses it: a 280 MB seed table
// is then paged in on demand and SHARED across processes rather than copied into each.
struct Mapping {
    void*  addr = nullptr;
    size_t len  = 0;
#ifdef _WIN32
    HANDLE file = INVALID_HANDLE_VALUE;
    HANDLE mapping = nullptr;
#endif

    Mapping() = default;
    Mapping(const Mapping&) = delete;
    Mapping& operator=(const Mapping&) = delete;

    ~Mapping() {
        if (!addr) return;
#ifdef _WIN32
        UnmapViewOfFile(addr);
        if (mapping) CloseHandle(mapping);
        if (file != INVALID_HANDLE_VALUE) CloseHandle(file);
#else
        ::munmap(addr, len);
#endif
    }

    void open(const std::string& path) {
#ifdef _WIN32
        file = CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING,
                           FILE_ATTRIBUTE_NORMAL, nullptr);
        if (file == INVALID_HANDLE_VALUE)
            throw std::runtime_error("seqtree: cannot open '" + path + "'");
        LARGE_INTEGER sz;
        if (!GetFileSizeEx(file, &sz)) {
            CloseHandle(file);
            file = INVALID_HANDLE_VALUE;
            throw std::runtime_error("seqtree: cannot size '" + path + "'");
        }
        len = size_t(sz.QuadPart);
        mapping = CreateFileMappingA(file, nullptr, PAGE_READONLY, 0, 0, nullptr);
        if (mapping) addr = MapViewOfFile(mapping, FILE_MAP_READ, 0, 0, 0);
        if (!addr) {
            if (mapping) CloseHandle(mapping);
            CloseHandle(file);
            mapping = nullptr;
            file = INVALID_HANDLE_VALUE;
            throw std::runtime_error("seqtree: cannot map '" + path + "'");
        }
#else
        int fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) throw std::runtime_error("seqtree: cannot open '" + path + "'");
        struct stat st{};
        if (::fstat(fd, &st) != 0 || st.st_size <= 0) {
            ::close(fd);
            throw std::runtime_error("seqtree: cannot size '" + path + "'");
        }
        len = size_t(st.st_size);
        addr = ::mmap(nullptr, len, PROT_READ, MAP_SHARED, fd, 0);
        ::close(fd);  // the mapping keeps its own reference; the descriptor is not needed
        if (addr == MAP_FAILED) {
            addr = nullptr;
            throw std::runtime_error("seqtree: cannot map '" + path + "'");
        }
#endif
    }
};

}  // namespace

// Owns (or maps) the flat arrays and points the hot path at them. Exactly one of the two
// ownership modes is populated; the const views below are what search actually reads, so the
// query code never branches on where the bytes came from.
struct TextStore {
    Alphabet alphabet = Alphabet::AminoAcid;
    uint8_t  A = 0;  // codec size
    uint8_t  k = 4;
    uint32_t num_refs = 0;
    uint64_t text_len = 0;      // includes one sentinel after every record
    uint64_t num_residues = 0;  // excludes them
    uint64_t num_buckets = 0;   // A^k
    uint64_t num_post = 0;
    uint64_t num_unknown = 0;  // text residues outside the codec; holes, never matches
    bool     has_group = false;

    const uint8_t*  text = nullptr;
    const uint32_t* starts = nullptr;      // num_refs + 1; 32-bit, build caps text at 2^32
    const uint32_t* post_begin = nullptr;  // num_buckets + 1
    const uint32_t* post_ids = nullptr;    // num_post
    const uint32_t* group = nullptr;       // num_refs, null when unused

    std::vector<uint8_t>  text_own;
    std::vector<uint32_t> starts_own, post_begin_own, post_ids_own, group_own;
    Mapping map;

    void point_at_owned() {
        text = text_own.data();
        starts = starts_own.data();
        post_begin = post_begin_own.data();
        post_ids = post_ids_own.data();
        group = group_own.empty() ? nullptr : group_own.data();
    }

    uint32_t ref_of(uint32_t abs_pos) const {
        return uint32_t(std::upper_bound(starts, starts + num_refs + 1, abs_pos) - starts - 1);
    }
};

namespace {

// Saturates rather than wrapping: a wrapped A^k could slip under the bucket cap and allocate
// a table that does not address the alphabet.
uint64_t pow_checked(uint8_t A, uint8_t k) {
    uint64_t n = 1;
    for (uint8_t i = 0; i < k; ++i) {
        if (n > (uint64_t(1) << 60) / A) return ~uint64_t(0);
        n *= A;
    }
    return n;
}

// Encode with the gapblock error form (0.6.1): name the side, the index, and the string, so a
// single malformed row in a batch of hundreds of thousands is findable.
void encode_into(const Codec& c, const std::string& s, const char* label, size_t idx,
                 std::vector<uint8_t>& out) {
    for (size_t i = 0; i < s.size(); ++i) {
        uint8_t v = c.encode(s[i]);
        if (v == Codec::kInvalid)
            throw std::invalid_argument(label + std::string("[") + std::to_string(idx) + "] ('" +
                                        s + "'): symbol '" + s[i] +
                                        "' is not in the alphabet");
        out.push_back(v);
    }
}

// The TEXT is encoded with a different rule from a query. A real proteome carries residues the
// codec does not name -- 36 U (selenocysteine) in human, 33 in mouse -- and refusing the whole
// build over them would make the class unusable on its own reference data. They are kept as
// kInvalid instead, which makes them behave exactly like the inter-record separator: no seed
// spans one and no hit crosses one. Counted, not silent: build reports how many it found.
uint64_t encode_text_into(const Codec& c, const std::string& s, std::vector<uint8_t>& out) {
    uint64_t unknown = 0;
    for (char ch : s) {
        uint8_t v = c.encode(ch);
        unknown += (v == Codec::kInvalid);
        out.push_back(v);
    }
    return unknown;
}

// Every <= m proper-substitution variant of q[0..k), as packed k-mer codes. Duplicate-free by
// construction, so this appends straight into `out` with no sort and no set.
void enumerate_ball(const uint8_t* q, uint8_t k, uint8_t A, uint16_t m, const uint64_t* pw,
                    uint8_t pos, uint32_t code, uint16_t used, std::vector<uint32_t>& out) {
    if (pos == k) {
        out.push_back(code);
        return;
    }
    const uint32_t place = uint32_t(pw[k - 1 - pos]);
    enumerate_ball(q, k, A, m, pw, pos + 1, code + q[pos] * place, used, out);
    if (used < m)
        for (uint8_t a = 0; a < A; ++a)
            if (a != q[pos])
                enumerate_ball(q, k, A, m, pw, pos + 1, code + a * place, used + 1, out);
}

uint32_t kmer_code(const uint8_t* p, uint8_t k, uint8_t A) {
    uint32_t code = 0;
    for (uint8_t i = 0; i < k; ++i) code = code * A + p[i];
    return code;
}

// A verified hit, before the CSR flatten. `mm_off` indexes a per-query mismatch buffer that is
// appended to in discovery order, so sorting these carries the detail along without touching it.
struct RawHit {
    uint32_t ref_id, offset;
    uint16_t n_subs;
    uint16_t n_ins, n_dels, length;
    int32_t  score;
    uint32_t mm_off;
    // How many mismatch entries this hit owns at mm_off. Equal to n_subs on the substitution
    // path; 0 on the indel path, which reports counts without per-column detail because
    // recovering the substituted columns would need an alignment traceback. Reading n_subs
    // entries instead would walk off the end of a buffer the indel path never wrote to.
    uint16_t mm_len;
    uint32_t abs;
};

// Reused across every query a worker handles, so the probe buffer is allocated once per thread
// rather than once per query. There is no candidate buffer: a posting is verified the moment it
// is read, which is also what removes the sort the old two-path dispatch needed to deduplicate.
struct Work {
    std::vector<uint32_t> probes;  // k-mer codes to look up
    // Banded-DP scratch for the indel path, allocated once per worker. Indexed
    // [d + max_indels][g]; see verify_banded.
    std::vector<uint16_t> dp_cur, dp_next;
};

// Per-query results. Mismatch positions are appended in discovery order and referenced by
// RawHit::mm_off, so sorting the hits carries their detail along without reordering it.
struct QueryOut {
    std::vector<RawHit>    hits;
    std::vector<uint16_t>  mm_pos;
    std::vector<TextGroup> groups;
    uint8_t                truncated = 0;
};

// Banded alignment of the whole query against the text starting at `s`, under two independent
// caps: at most `ms` substitutions and at most `mi` insertions+deletions.
//
// State is (i, d, g): i query residues consumed, d = (text consumed) - (query consumed), and g
// indels spent. d is confined to [-mi, mi] and g to [0, mi], and every move changes them
// together -- an insertion (a query residue with no text opposite it) is d-1, g+1; a deletion
// (a text residue with none opposite) is d+1, g+1 -- so g >= |d| and g == |d| (mod 2). The cell
// holds the fewest substitutions reaching that state, or kInf.
//
// The alignment must BEGIN and END on an aligned pair -- no leading or trailing gap of either
// kind. Two reasons, and they agree. A gap at either edge extends the reported interval over a
// residue that matches nothing, so `offset` and `length` would no longer bound the part of the
// text the query actually explains. And a leading gap is just the same alignment starting one
// residue over, which the caller already reaches by sweeping every start in the seed's window,
// so allowing it would report one occurrence two or three times over.
//
// Every hit that survives is appended via `emit(end_offset_from_s, n_subs, n_ins, n_dels)`. One
// start can yield several ends -- they are genuinely different intervals -- and the caller
// deduplicates on (ref, offset, length).
template <class Emit>
void verify_banded(const uint8_t* q, size_t L, const uint8_t* text, uint64_t text_len, uint64_t s,
                   uint16_t ms, uint16_t mi, std::vector<uint16_t>& cur,
                   std::vector<uint16_t>& next, Emit emit) {
    constexpr uint16_t kInf = 0xFFFF;
    const size_t W = size_t(2) * mi + 1, G = size_t(mi) + 1, N = W * G;
    const auto at = [mi, G](int d, size_t g) { return size_t(d + int(mi)) * G + g; };

    cur.assign(N, kInf);
    next.assign(N, kInf);
    cur[at(0, 0)] = 0;  // no leading gap: the alignment starts flush at s

    for (size_t i = 0; i < L; ++i) {
        // Deletions consume text without consuming a query residue, so they move within this
        // same i. Relax in increasing d: a run of them chains, each costing one unit of g.
        // Not at i == 0: that would be a leading gap.
        for (int d = -int(mi); i > 0 && d < int(mi); ++d) {
            for (size_t g = 0; g + 1 < G; ++g) {
                const uint16_t v = cur[at(d, g)];
                if (v == kInf) continue;
                const int64_t tpos = int64_t(s) + int64_t(i) + d;
                if (tpos < 0 || uint64_t(tpos) >= text_len) continue;
                if (text[tpos] == Codec::kInvalid) continue;  // never match through a hole
                uint16_t& dst = cur[at(d + 1, g + 1)];
                if (v < dst) dst = v;
            }
        }

        std::fill(next.begin(), next.end(), kInf);
        bool any = false;
        for (int d = -int(mi); d <= int(mi); ++d) {
            for (size_t g = 0; g < G; ++g) {
                const uint16_t v = cur[at(d, g)];
                if (v == kInf) continue;
                // match / substitution: consume q[i] against text[s + i + d]
                const int64_t tpos = int64_t(s) + int64_t(i) + d;
                if (tpos >= 0 && uint64_t(tpos) < text_len && text[tpos] != Codec::kInvalid) {
                    const uint16_t nv = uint16_t(v + (q[i] != text[tpos] ? 1 : 0));
                    if (nv <= ms) {
                        uint16_t& dst = next[at(d, g)];
                        if (nv < dst) dst = nv;
                        any = true;
                    }
                }
                // Insertion: consume q[i] with no text opposite it. Never at the first or last
                // query residue -- that is a leading or trailing gap.
                if (i > 0 && i + 1 < L && d > -int(mi) && g + 1 < G) {
                    uint16_t& dst = next[at(d - 1, g + 1)];
                    if (v < dst) dst = v;
                    any = true;
                }
            }
        }
        // Every surviving state already exceeds one of the caps, so no suffix can rescue this
        // start. Without this the DP walked all L rows on a candidate the Hamming path would
        // have rejected after two residues, and most candidates are rejected: on a 5 M-residue
        // corpus it is the difference between 0.107 and 0.019 ms/query at L=16.
        if (!any) return;
        cur.swap(next);
    }

    // The query is fully consumed. Every finite state is a complete alignment; d fixes the
    // matched length and g splits into insertions and deletions, since d = dels - ins and
    // g = dels + ins.
    for (int d = -int(mi); d <= int(mi); ++d) {
        for (size_t g = 0; g < G; ++g) {
            const uint16_t v = cur[at(d, g)];
            if (v == kInf || v > ms) continue;
            if (size_t(std::abs(d)) > g || ((g - size_t(std::abs(d))) & 1)) continue;
            const size_t dels = (g + size_t(d)) / 2, ins = (g - size_t(d)) / 2;
            if (int64_t(L) + d < 0) continue;
            emit(size_t(int64_t(L) + d), v, uint16_t(ins), uint16_t(dels));
        }
    }
}

// One query against the text at a fixed max_subs. Appends to out.hits / out.mm_pos.
// The scheme and why it is lossless are at the top of this file.
// The indel path. Seeding is unchanged and every block is exact -- with b = m + mi + 1 blocks
// the total edit count cannot reach the block count, so pigeonhole still puts a 0-error block
// somewhere and the existing exact lookup finds it. Only verification differs: the seed pins a
// text position, not a start, because the prefix before the block may have gained or lost up to
// mi residues, so every start in [pos - base - mi, pos - base + mi] is aligned.
//
// Hits are appended unsorted and with duplicates; the caller sorts by (n_subs, ref_id, offset,
// length) and drops adjacent equals, which is the same sort the substitution path already needs
// for its stable output order.
void search_indels(const TextStore& st, const uint8_t* q, size_t L, uint16_t m, uint16_t mi,
                   const SubstitutionMatrix* matrix, Work& s, QueryOut& out) {
    const uint8_t k = st.k;
    const size_t b = size_t(m) + size_t(mi) + 1;  // search_batch enforces L >= b * k
    const size_t bw = L / b;

    for (size_t j = 0; j < b; ++j) {
        const size_t base = j * bw;
        uint32_t code = 0;
        for (uint8_t x = 0; x < k; ++x) code = uint32_t(code * st.A + q[base + x]);
        const uint32_t lo = st.post_begin[code], hi = st.post_begin[code + 1];
        for (uint32_t i = lo; i < hi; ++i) {
            const uint32_t pos = st.post_ids[i];
            const int64_t anchor = int64_t(pos) - int64_t(base);
            for (int64_t start = anchor - mi; start <= anchor + mi; ++start) {
                if (start < 0 || uint64_t(start) >= st.text_len) continue;
                verify_banded(
                    q, L, st.text, st.text_len, uint64_t(start), m, mi, s.dp_cur, s.dp_next,
                    [&](size_t len, uint16_t n_subs, uint16_t n_ins, uint16_t n_dels) {
                        if (uint64_t(start) + len > st.text_len) return;
                        const uint32_t c = uint32_t(start);
                        const uint32_t ref = st.ref_of(c);
                        // A hit must lie wholly inside one record. The sentinels already stop
                        // the alignment, but a zero-length tail would slip through.
                        if (uint64_t(c) + len > st.starts[ref + 1]) return;
                        int32_t score = 0;
                        if (matrix) {
                            // No mismatch positions on this path (see the comment on mm_pos in
                            // search_batch), so the matrix scores what it can: the substituted
                            // columns of the best ungapped prefix are not recoverable without a
                            // traceback, and a partial score would be worse than none.
                            score = 0;
                        }
                        out.hits.push_back(RawHit{ref, c - st.starts[ref], n_subs, n_ins, n_dels,
                                                  uint16_t(len), score,
                                                  uint32_t(out.mm_pos.size()), 0, c});
                    });
            }
        }
    }
}

void search_one(const TextStore& st, const uint8_t* q, size_t L, uint16_t m,
                const SubstitutionMatrix* matrix, const uint64_t* pw, Work& s, QueryOut& out) {
    const uint8_t k = st.k, A = st.A;
    const size_t b  = std::min<size_t>(size_t(m) + 1, L / k);  // >= 1: search_batch enforces L >= k
    const size_t bw = L / b;
    const size_t r  = size_t(m) + 1 > b ? size_t(m) + 1 - b : 0;
    // Spare budget units go to the LAST blocks, and that is worth 1.4x on the human proteome
    // (L=8, m=2, k=4: 1.51 -> 1.09 ms/query). The lemma is symmetric in j, so this is free to
    // choose; verification is not symmetric. A candidate produced by block j matches its own
    // k-mer, so scanning left-to-right hits the UNCONSTRAINED prefix [0, j*bw) first and
    // early-exits after a residue or two -- and it is the high-budget block that contributes
    // almost all the candidates, so that is the one whose rejects must be cheap.
    auto budget = [r, b](size_t j) { return uint16_t(r / b + (j >= b - r % b ? 1 : 0)); };

    for (size_t j = 0; j < b; ++j) {
        const size_t base = j * bw;
        s.probes.clear();
        enumerate_ball(q + base, k, A, budget(j), pw, 0, 0, 0, s.probes);
        for (uint32_t code : s.probes) {
            const uint32_t lo = st.post_begin[code], hi = st.post_begin[code + 1];
            for (uint32_t i = lo; i < hi; ++i) {
                const uint32_t pos = st.post_ids[i];
                if (pos < base) continue;
                const uint32_t c = uint32_t(pos - base);
                if (uint64_t(c) + L > st.text_len) continue;
                const uint8_t* t = st.text + c;
                const uint32_t mm_off = uint32_t(out.mm_pos.size());
                uint16_t n = 0;
                bool ok = true;
                for (size_t x = 0; x < L; ++x) {
                    const uint8_t ch = t[x];
                    if (ch == Codec::kInvalid) { ok = false; break; }  // never match through a hole
                    if (ch == q[x]) continue;
                    if (n == m) { ok = false; break; }  // early exit: already too far
                    ++n;
                    out.mm_pos.push_back(uint16_t(x));
                }
                // Deduplicate across blocks without a sort: block i produces this same start
                // exactly when its leading k-mer carries <= c_i mismatches, and the positions
                // just written already name them. Emit only from the lowest such block.
                for (size_t i = 0; ok && i < j; ++i) {
                    uint16_t e = 0;
                    for (uint16_t x = 0; x < n; ++x)
                        e += size_t(out.mm_pos[mm_off + x]) - i * bw < k;  // unsigned wrap: p < i*bw
                    ok = e > budget(i);
                }
                if (!ok) {
                    out.mm_pos.resize(mm_off);
                    continue;
                }
                const uint32_t ref = st.ref_of(c);
                int32_t score = 0;
                if (matrix)
                    for (uint32_t x = mm_off; x < mm_off + n; ++x)
                        score += matrix->similarity(q[out.mm_pos[x]], t[out.mm_pos[x]]);
                out.hits.push_back(RawHit{ref, c - st.starts[ref], n, 0, 0, uint16_t(L), score,
                                          mm_off, n, c});
            }
        }
    }
}

}  // namespace

TextIndex::TextIndex() = default;
TextIndex::~TextIndex() = default;

uint32_t TextIndex::num_refs() const { return store_->num_refs; }
uint64_t TextIndex::num_residues() const { return store_->num_residues; }
uint64_t TextIndex::num_unknown() const { return store_->num_unknown; }
uint8_t TextIndex::k() const { return store_->k; }
Alphabet TextIndex::alphabet() const { return store_->alphabet; }

// Residues the codec never named come back as 'X'. That is LOSSY and 'X' is a real symbol in
// this codec, so a record containing U round-trips to a string that would search differently
// from the one that was indexed; ref_seq reports what the index holds, not what was handed in.
std::string TextIndex::ref_seq(uint32_t ref_id) const {
    const TextStore& st = *store_;
    if (ref_id >= st.num_refs)
        throw std::invalid_argument("seqtree: ref_id " + std::to_string(ref_id) +
                                    " out of range (" + std::to_string(st.num_refs) + " records)");
    Codec c(st.alphabet);
    const uint64_t b = st.starts[ref_id], e = st.starts[ref_id + 1] - 1;  // -1: drop the sentinel
    std::string out;
    out.reserve(size_t(e - b));
    for (uint64_t i = b; i < e; ++i)
        out += st.text[i] == Codec::kInvalid ? 'X' : c.decode(st.text[i]);
    return out;
}

std::unique_ptr<TextIndex> TextIndex::build(const std::vector<std::string>& refs, Alphabet alphabet,
                                            uint8_t k, const std::vector<uint32_t>& group_ids) {
    Codec codec(alphabet);
    const uint8_t A = codec.size();
    if (k == 0) throw std::invalid_argument("seqtree: k must be >= 1");
    const uint64_t buckets = pow_checked(A, k);
    if (buckets > kMaxBuckets)
        throw std::invalid_argument(
            "seqtree: k=" + std::to_string(k) + " over a " + std::to_string(A) +
            "-symbol alphabet needs " + std::to_string(buckets) +
            " seed buckets, above the " + std::to_string(kMaxBuckets) + " cap; use a smaller k");
    if (!group_ids.empty() && group_ids.size() != refs.size())
        throw std::invalid_argument("seqtree: group_ids has " + std::to_string(group_ids.size()) +
                                    " entries for " + std::to_string(refs.size()) + " records");

    auto ix = std::unique_ptr<TextIndex>(new TextIndex());
    ix->store_ = std::make_unique<TextStore>();
    TextStore& st = *ix->store_;
    st.alphabet = alphabet;
    st.A = A;
    st.k = k;
    st.num_refs = uint32_t(refs.size());
    st.num_buckets = buckets;

    uint64_t total = 0;
    for (const auto& r : refs) total += r.size() + 1;  // + one sentinel per record
    if (total >= uint64_t(1) << 32)
        throw std::invalid_argument("seqtree: text is " + std::to_string(total) +
                                    " symbols; TextIndex postings are 32-bit and cap it at 2^32");

    st.text_own.reserve(size_t(total));
    st.starts_own.reserve(refs.size() + 1);
    for (size_t i = 0; i < refs.size(); ++i) {
        st.starts_own.push_back(uint32_t(st.text_own.size()));
        st.num_unknown += encode_text_into(codec, refs[i], st.text_own);
        st.text_own.push_back(Codec::kInvalid);
    }
    st.starts_own.push_back(uint32_t(st.text_own.size()));
    st.text_len = st.text_own.size();
    st.num_residues = st.text_len - refs.size();
    if (!group_ids.empty()) {
        st.group_own = group_ids;
        st.has_group = true;
    }

    // Counting sort in two passes over the text: count per bucket, prefix-sum, then scatter.
    // No comparison sort, no hashing, no argsort -- the argsort is what dominated the 42 s
    // NumPy prototype build.
    const uint64_t mod = buckets;  // A^k; the rolling code is kept modulo this
    st.post_begin_own.assign(size_t(buckets) + 1, 0);
    auto roll = [&](auto&& emit) {
        uint32_t code = 0;
        uint32_t run = 0;
        for (uint64_t p = 0; p < st.text_len; ++p) {
            const uint8_t ch = st.text_own[size_t(p)];
            if (ch == Codec::kInvalid) {
                run = 0;
                code = 0;
                continue;
            }
            code = uint32_t((uint64_t(code) * A + ch) % mod);
            if (++run >= k) emit(code, uint32_t(p - k + 1));
        }
    };
    roll([&](uint32_t code, uint32_t) { ++st.post_begin_own[code + 1]; });
    for (uint64_t b = 0; b < buckets; ++b) st.post_begin_own[size_t(b) + 1] += st.post_begin_own[size_t(b)];
    st.num_post = st.post_begin_own.back();
    st.post_ids_own.resize(size_t(st.num_post));
    std::vector<uint32_t> cursor(st.post_begin_own.begin(), st.post_begin_own.end() - 1);
    roll([&](uint32_t code, uint32_t pos) { st.post_ids_own[cursor[code]++] = pos; });
    // The two roll() passes must visit exactly the same positions or the scatter overruns its
    // bucket and silently corrupts a neighbour's postings. Nothing else ties them together.
    for (uint64_t bkt = 0; bkt < buckets; ++bkt)
        if (cursor[size_t(bkt)] != st.post_begin_own[size_t(bkt) + 1])
            throw std::logic_error("seqtree: TextIndex seed-table scatter disagreed with its count");

    st.point_at_owned();
    return ix;
}

TextResult TextIndex::search_batch(const std::vector<std::string>& queries,
                                   const TextQueryOpts& opts, int threads) const {
    const TextStore& st = *store_;
    Codec codec(st.alphabet);
    if (opts.matrix && opts.matrix->size() != st.A)
        throw std::invalid_argument("seqtree: matrix size does not match the alphabet");

    // Encode and validate every query up front, in this thread. The workers then cannot throw,
    // which is why they need no exception plumbing (same reasoning as gapblock.cpp).
    const size_t n = queries.size();
    std::vector<uint8_t> qbuf;
    std::vector<uint64_t> qoff(n + 1, 0);
    for (size_t i = 0; i < n; ++i) {
        if (queries[i].size() < st.k)
            throw std::invalid_argument(
                "queries[" + std::to_string(i) + "] ('" + queries[i] + "') is " +
                std::to_string(queries[i].size()) + " long; TextIndex indexes seeds of k=" +
                std::to_string(st.k) + " and cannot answer a shorter query completely");
        if (opts.max_indels > 0) {
            const size_t need = (size_t(opts.max_subs) + opts.max_indels + 1) * st.k;
            if (queries[i].size() < need)
                throw std::invalid_argument(
                    "queries[" + std::to_string(i) + "] ('" + queries[i] + "') is " +
                    std::to_string(queries[i].size()) + " long; searching with max_indels=" +
                    std::to_string(opts.max_indels) + " and max_subs=" +
                    std::to_string(opts.max_subs) + " needs every one of the " +
                    std::to_string(size_t(opts.max_subs) + opts.max_indels + 1) +
                    " seed blocks to be exact, so a query must be at least " +
                    std::to_string(need) + " long at k=" + std::to_string(st.k) +
                    ". Rebuild the index with a smaller k, or lower max_subs/max_indels");
        }
        if (queries[i].size() > kMaxQueryLen)
            throw std::invalid_argument(
                "queries[" + std::to_string(i) + "] is " + std::to_string(queries[i].size()) +
                " long; TextIndex reports mismatch positions as 16-bit and caps a query at " +
                std::to_string(kMaxQueryLen));
        encode_into(codec, queries[i], "queries", i, qbuf);
        qoff[i + 1] = qbuf.size();
    }

    std::vector<uint64_t> pw(size_t(st.k) + 1);
    pw[0] = 1;
    for (uint8_t i = 1; i <= st.k; ++i) pw[i] = pw[i - 1] * st.A;

    std::vector<QueryOut> per_query(n);

    auto run_query = [&](size_t i, Work& w) {
        QueryOut& s = per_query[i];
        const uint8_t* q = qbuf.data() + qoff[i];
        const size_t L = size_t(qoff[i + 1] - qoff[i]);
        const uint16_t top = opts.max_subs;
        // best_only walks m upward and stops at the first shell that yields anything, so a query
        // resolved at m = 0 never enters the m = 1 pass. No second index, no second table.
        for (uint16_t m = opts.best_only ? 0 : top; m <= top; ++m) {
            s.hits.clear();
            s.mm_pos.clear();
            if (opts.max_indels > 0)
                search_indels(st, q, L, m, opts.max_indels, opts.matrix, w, s);
            else
                search_one(st, q, L, m, opts.matrix, pw.data(), w, s);
            if (opts.exclude_exact)
                s.hits.erase(std::remove_if(s.hits.begin(), s.hits.end(),
                                            [](const RawHit& h) { return h.n_subs == 0; }),
                             s.hits.end());
            if (!opts.best_only || !s.hits.empty()) break;
        }
        // Deterministic and independent of the order candidates were discovered in. The indel
        // path adds length to the key and then deduplicates: one text interval is reachable
        // from several blocks and several starts, and the cheapest alignment of that interval
        // is the one to keep, so sorting by (n_subs, n_ins + n_dels) first puts it in front.
        std::sort(s.hits.begin(), s.hits.end(), [](const RawHit& a, const RawHit& b) {
            if (a.n_subs != b.n_subs) return a.n_subs < b.n_subs;
            if (a.n_ins + a.n_dels != b.n_ins + b.n_dels)
                return a.n_ins + a.n_dels < b.n_ins + b.n_dels;
            if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
            if (a.offset != b.offset) return a.offset < b.offset;
            return a.length < b.length;
        });
        if (opts.max_indels > 0) {
            std::vector<RawHit> keep;
            keep.reserve(s.hits.size());
            for (const RawHit& h : s.hits) {
                const bool seen = std::any_of(keep.begin(), keep.end(), [&](const RawHit& x) {
                    return x.ref_id == h.ref_id && x.offset == h.offset && x.length == h.length;
                });
                if (!seen) keep.push_back(h);
            }
            s.hits.swap(keep);
            // The dedup pass above is a linear scan per hit, which is fine while the kept set
            // is the handful a real query returns; max_hits bounds it in the pathological case.
            std::sort(s.hits.begin(), s.hits.end(), [](const RawHit& a, const RawHit& b) {
                if (a.n_subs != b.n_subs) return a.n_subs < b.n_subs;
                if (a.ref_id != b.ref_id) return a.ref_id < b.ref_id;
                if (a.offset != b.offset) return a.offset < b.offset;
                return a.length < b.length;
            });
        }
        // A cap keeps the BEST hits, because it is applied after the sort -- and it is always
        // reported. A cap that is invisible is a recall bug wearing a performance costume.
        bool capped = opts.max_hits && s.hits.size() > opts.max_hits;
        if (capped) s.hits.resize(opts.max_hits);

        s.groups.clear();
        if (opts.group_by && st.group) {
            for (const RawHit& h : s.hits) {
                const uint32_t g = st.group[h.ref_id];
                auto it = std::find_if(s.groups.begin(), s.groups.end(),
                                       [g](const TextGroup& r) { return r.group_id == g; });
                if (it == s.groups.end())
                    s.groups.push_back(TextGroup{g, h.n_subs, 1});
                else {
                    it->min_subs = std::min(it->min_subs, h.n_subs);
                    ++it->n_hits;
                }
            }
            std::sort(s.groups.begin(), s.groups.end(),
                      [](const TextGroup& a, const TextGroup& b) { return a.group_id < b.group_id; });
        }
        s.truncated = capped ? 1 : 0;
    };

    if (n == 0) return TextResult{};

    // Chunk floor 8 rather than parallel_for's default 1: QueryOut is ~80 bytes, so a chunk of
    // 1 has neighbouring threads writing into the same cache line for the whole batch. `Work`
    // is the per-worker probe buffer, allocated once per thread.
    parallel_for(n, threads, [] { return Work{}; },
                 [&](size_t i, Work& w) { run_query(i, w); }, size_t(8), size_t(1024));

    // The flatten is serial, and after the search scheme cut per-query time it is a real share
    // of the batch -- so size every array once from what the workers actually produced rather
    // than growing nine vectors by doubling.
    size_t n_hits = 0, n_mm = 0, n_groups = 0;
    for (const QueryOut& s : per_query) {
        n_hits += s.hits.size();
        n_groups += s.groups.size();
        for (const RawHit& h : s.hits) n_mm += h.mm_len;
    }
    TextResult out;
    out.query_begin.reserve(n + 1);
    out.truncated.reserve(n);
    out.ref_id.reserve(n_hits);
    out.offset.reserve(n_hits);
    out.n_subs.reserve(n_hits);
    out.n_ins.reserve(n_hits);
    out.n_dels.reserve(n_hits);
    out.length.reserve(n_hits);
    out.score.reserve(n_hits);
    out.mm_begin.reserve(n_hits + 1);
    out.mm_pos.reserve(n_mm);
    out.mm_query_aa.reserve(n_mm);
    out.mm_text_aa.reserve(n_mm);
    if (opts.group_by) {
        out.group_begin.reserve(n + 1);
        out.group_id.reserve(n_groups);
        out.group_min_subs.reserve(n_groups);
        out.group_n_hits.reserve(n_groups);
    }
    out.query_begin.push_back(0);
    out.mm_begin.push_back(0);
    if (opts.group_by) out.group_begin.push_back(0);
    for (size_t i = 0; i < n; ++i) {
        QueryOut& s = per_query[i];
        for (const RawHit& h : s.hits) {
            out.ref_id.push_back(h.ref_id);
            out.offset.push_back(h.offset);
            out.n_subs.push_back(h.n_subs);
            out.n_ins.push_back(h.n_ins);
            out.n_dels.push_back(h.n_dels);
            out.length.push_back(h.length);
            out.score.push_back(h.score);
            const uint8_t* t = st.text + h.abs;
            for (uint16_t j = 0; j < h.mm_len; ++j) {
                const uint16_t p = s.mm_pos[h.mm_off + j];
                out.mm_pos.push_back(p);
                out.mm_query_aa.push_back(codec.decode(qbuf[qoff[i] + p]));
                out.mm_text_aa.push_back(codec.decode(t[p]));
            }
            out.mm_begin.push_back(uint32_t(out.mm_pos.size()));
        }
        out.query_begin.push_back(uint32_t(out.ref_id.size()));
        out.truncated.push_back(s.truncated);
        if (opts.group_by) {
            for (const TextGroup& g : s.groups) {
                out.group_id.push_back(g.group_id);
                out.group_min_subs.push_back(g.min_subs);
                out.group_n_hits.push_back(g.n_hits);
            }
            out.group_begin.push_back(uint32_t(out.group_id.size()));
        }
    }
    return out;
}

void TextIndex::save(const std::string& path) const {
    const TextStore& st = *store_;
    detail::atomic_write(path, [&](std::ostream& os) {
        char head[kHeaderBytes] = {};
        std::memcpy(head, kMagic, 4);
        const uint32_t version = kVersion, num_refs = st.num_refs;
        const uint8_t alph = uint8_t(st.alphabet), kk = st.k, hg = st.has_group ? 1 : 0;
        std::memcpy(head + 4, &version, 4);
        std::memcpy(head + 8, &alph, 1);
        std::memcpy(head + 9, &kk, 1);
        std::memcpy(head + 10, &hg, 1);
        std::memcpy(head + 12, &num_refs, 4);
        std::memcpy(head + 16, &st.text_len, 8);
        std::memcpy(head + 24, &st.num_buckets, 8);
        std::memcpy(head + 32, &st.num_post, 8);
        std::memcpy(head + 40, &st.num_residues, 8);
        std::memcpy(head + 48, &st.num_unknown, 8);
        os.write(head, kHeaderBytes);
        auto put = [&os](const void* p, size_t bytes) {
            os.write(reinterpret_cast<const char*>(p), std::streamsize(bytes));
        };
        put(st.starts, (size_t(st.num_refs) + 1) * 4);
        put(st.post_begin, (size_t(st.num_buckets) + 1) * 4);
        put(st.post_ids, size_t(st.num_post) * 4);
        if (st.has_group) put(st.group, size_t(st.num_refs) * 4);
        put(st.text, size_t(st.text_len));
    });
}

std::unique_ptr<TextIndex> TextIndex::load(const std::string& path, bool mmap) {
    auto ix = std::unique_ptr<TextIndex>(new TextIndex());
    ix->store_ = std::make_unique<TextStore>();
    TextStore& st = *ix->store_;

    auto parse_header = [&](const char* head) {
        uint32_t version = 0;
        if (std::memcmp(head, kMagic, 4) != 0) version = 0;
        else std::memcpy(&version, head + 4, 4);
        if (version != kVersion)
            throw std::runtime_error("seqtree: '" + path + "' is not a v" +
                                     std::to_string(kVersion) + " seqtree text index");
        uint8_t alph = 0, kk = 0, hg = 0;
        std::memcpy(&alph, head + 8, 1);
        std::memcpy(&kk, head + 9, 1);
        std::memcpy(&hg, head + 10, 1);
        std::memcpy(&st.num_refs, head + 12, 4);
        std::memcpy(&st.text_len, head + 16, 8);
        std::memcpy(&st.num_buckets, head + 24, 8);
        std::memcpy(&st.num_post, head + 32, 8);
        std::memcpy(&st.num_residues, head + 40, 8);
        std::memcpy(&st.num_unknown, head + 48, 8);
        st.alphabet = Alphabet(alph);
        st.k = kk;
        st.has_group = hg != 0;
        st.A = Codec(st.alphabet).size();
        // A .sti is untrusted input: every size below is read from the file and then used for
        // resize() or, on the mmap path, for raw pointer arithmetic into the mapping. Cross-check
        // them against each other before anything is dereferenced.
        auto bad = [&](const char* why) {
            throw std::runtime_error("seqtree: corrupt text index '" + path + "': " +
                                     std::string(why));
        };
        if (st.A == 0 || st.k == 0) bad("alphabet or k is zero");
        if (st.num_buckets != pow_checked(st.A, st.k) || st.num_buckets > kMaxBuckets)
            bad("bucket count does not match the alphabet and k in the header");
        if (st.text_len >= (uint64_t(1) << 32)) bad("text length exceeds the 32-bit posting cap");
        if (st.num_post > st.text_len) bad("more postings than text positions");
        if (st.num_residues > st.text_len || st.num_unknown > st.text_len)
            bad("residue or unknown count exceeds the text length");
        if (uint64_t(st.num_refs) + 1 > st.text_len + 1) bad("more records than text positions");
    };
    // Checked once the arrays are readable, on either path.
    auto validate_arrays = [&] {
        if (st.starts[0] != 0 || st.starts[st.num_refs] != st.text_len)
            throw std::runtime_error("seqtree: corrupt text index '" + path +
                                     "': record offsets do not span the text");
        if (st.post_begin[0] != 0 || st.post_begin[st.num_buckets] != st.num_post)
            throw std::runtime_error("seqtree: corrupt text index '" + path +
                                     "': seed table does not span its postings");
    };
    // Sizes are all derived from the header, so a truncated file has to be caught explicitly
    // rather than by a short read -- with mmap there is no read to come up short.
    auto expected_bytes = [&] {
        return kHeaderBytes + (size_t(st.num_refs) + 1) * 4 + (size_t(st.num_buckets) + 1) * 4 +
               size_t(st.num_post) * 4 + (st.has_group ? size_t(st.num_refs) * 4 : 0) +
               size_t(st.text_len);
    };

    if (mmap) {
        st.map.open(path);
        if (st.map.len < kHeaderBytes)
            throw std::runtime_error("seqtree: '" + path + "' is too short to be a text index");
        const char* base = static_cast<const char*>(st.map.addr);
        parse_header(base);
        if (st.map.len < expected_bytes())
            throw std::runtime_error("seqtree: truncated or corrupt text index '" + path + "'");
        size_t at = kHeaderBytes;
        st.starts = reinterpret_cast<const uint32_t*>(base + at);
        at += (size_t(st.num_refs) + 1) * 4;
        st.post_begin = reinterpret_cast<const uint32_t*>(base + at);
        at += (size_t(st.num_buckets) + 1) * 4;
        st.post_ids = reinterpret_cast<const uint32_t*>(base + at);
        at += size_t(st.num_post) * 4;
        if (st.has_group) {
            st.group = reinterpret_cast<const uint32_t*>(base + at);
            at += size_t(st.num_refs) * 4;
        }
        st.text = reinterpret_cast<const uint8_t*>(base + at);
        validate_arrays();
        return ix;
    }

    std::ifstream is(path, std::ios::binary);
    if (!is) throw std::runtime_error("seqtree: cannot open '" + path + "'");
    char head[kHeaderBytes] = {};
    is.read(head, kHeaderBytes);
    if (!is) throw std::runtime_error("seqtree: '" + path + "' is too short to be a text index");
    parse_header(head);
    auto get = [&is](void* p, size_t bytes) {
        is.read(reinterpret_cast<char*>(p), std::streamsize(bytes));
    };
    st.starts_own.resize(size_t(st.num_refs) + 1);
    get(st.starts_own.data(), st.starts_own.size() * 4);
    st.post_begin_own.resize(size_t(st.num_buckets) + 1);
    get(st.post_begin_own.data(), st.post_begin_own.size() * 4);
    st.post_ids_own.resize(size_t(st.num_post));
    get(st.post_ids_own.data(), st.post_ids_own.size() * 4);
    if (st.has_group) {
        st.group_own.resize(st.num_refs);
        get(st.group_own.data(), st.group_own.size() * 4);
    }
    st.text_own.resize(size_t(st.text_len));
    get(st.text_own.data(), st.text_own.size());
    if (!is) throw std::runtime_error("seqtree: truncated or corrupt text index '" + path + "'");
    st.point_at_owned();
    validate_arrays();
    return ix;
}

}  // namespace seqtree
