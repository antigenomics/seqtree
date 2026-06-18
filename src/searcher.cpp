#include "engines.hpp"
#include "seqtree/seqtree.hpp"

#include <algorithm>
#include <limits>
#include <stdexcept>

namespace seqtree {

Limits resolve_limits(const SearchParams& p) {
    Limits L{};
    L.mat = p.matrix;
    L.unit = (p.matrix == nullptr);
    L.gap = p.gap_open;
    L.max_sub = p.max_substitutions;
    L.max_ins = p.max_insertions;
    L.max_del = p.max_deletions;
    // max_total_edits is an independent total cap (0 => no total cap, fall back
    // to the per-type sum). It is NOT clamped by the per-type sum, so seqtrie --
    // which ignores per-type caps -- can be driven by max_total_edits alone.
    int sum = int(p.max_substitutions) + int(p.max_insertions) + int(p.max_deletions);
    L.max_tot = (p.max_total_edits > 0) ? int(p.max_total_edits) : sum;
    if (p.max_score_penalty > 0) L.budget = p.max_score_penalty;
    else if (L.unit) L.budget = L.max_tot;  // unit cost: one penalty per edit
    else L.budget = std::numeric_limits<int32_t>::max() / 4;  // matrix mode: rely on count caps
    L.max_hits = p.max_hits;
    return L;
}

namespace {
Engine pick_engine(const SearchParams& p) {
    if (p.engine != Engine::Auto) return p.engine;
    if (p.max_insertions == 0 && p.max_deletions == 0) return Engine::SeqTm;  // Hamming fast path
    if (p.matrix != nullptr) return Engine::SeqTrie;                          // matrix-weighted budget
    return Engine::SeqTm;                                                     // small-k indel
}

bool by_score(const Hit& a, const Hit& b) {
    if (a.score != b.score) return a.score < b.score;
    return a.ref_id < b.ref_id;
}
}  // namespace

Searcher::Searcher(const Index& idx) : idx_(idx), scratch_(std::make_unique<Scratch>()) {}
Searcher::~Searcher() = default;

void Searcher::search_into(std::string_view query, const SearchParams& p, std::vector<Hit>& out) {
    const Codec& cod = idx_.codec();
    qcodes_.resize(query.size());
    for (size_t i = 0; i < query.size(); ++i) {
        uint8_t c = cod.encode(query[i]);
        if (c == Codec::kInvalid)
            throw std::invalid_argument(std::string("invalid symbol '") + query[i] + "' in query");
        qcodes_[i] = c;
    }

    Limits L = resolve_limits(p);
    out.clear();
    scratch_->collisions = 0;  // seqtrie leaves it 0; seqtm accumulates per re-reached ref
    if (pick_engine(p) == Engine::SeqTm)
        search_seqtm(idx_.trie(), qcodes_.data(), int(query.size()), L, p.mode, *scratch_, out);
    else
        search_seqtrie(idx_.trie(), qcodes_.data(), int(query.size()), L, *scratch_, out);

    if (p.mode == Mode::TopHit) {
        uint32_t keep = p.max_hits > 0 ? p.max_hits : 1;
        if (out.size() > keep) {
            std::partial_sort(out.begin(), out.begin() + keep, out.end(), by_score);
            out.resize(keep);
        } else {
            std::sort(out.begin(), out.end(), by_score);
        }
    } else if (p.max_hits > 0 && out.size() > p.max_hits) {
        out.resize(p.max_hits);
    }
}

std::vector<Hit> Searcher::search(std::string_view query, const SearchParams& p) {
    std::vector<Hit> out;
    search_into(query, p, out);
    return out;
}

uint64_t Searcher::last_collisions() const { return scratch_->collisions; }

bool Searcher::search_top(std::string_view query, const SearchParams& p, Hit& out) {
    SearchParams q = p;
    q.mode = Mode::TopHit;
    q.max_hits = 1;
    std::vector<Hit> hits;
    search_into(query, q, hits);
    if (hits.empty()) return false;
    out = hits[0];
    return true;
}

}  // namespace seqtree
