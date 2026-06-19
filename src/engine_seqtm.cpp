#include "engines.hpp"

#include <limits>

namespace seqtree {

namespace {

struct TmCtx {
    const Trie& trie;
    const uint8_t* q;
    int qlen;
    Limits lim;
    bool hamming;          // no indels allowed
    bool posmode;          // per-position penalties active (Hamming, width == qlen)
    Mode mode;
    std::unordered_map<uint32_t, uint32_t>* seen;
    uint64_t* collisions;
    std::vector<Hit>* out;

    int32_t sub_pen(int qpos, uint8_t a, uint8_t b) const {
        if (posmode) return lim.posmat->penalty(uint16_t(qpos), a, b);
        return lim.unit ? (a == b ? 0 : 1) : lim.mat->penalty(a, b);
    }
    // A masked position (anchor) is free: it neither adds penalty nor counts as a substitution.
    bool pos_masked(int qpos) const { return posmode && lim.posmat->masked(uint16_t(qpos)); }
};

void emit(TmCtx& c, uint32_t node, int32_t pen, int ns, int ni, int nd) {
    const Trie::Node& N = c.trie.nodes[node];
    for (uint32_t k = 0; k < N.ref_count; ++k) {
        uint32_t ref_id = c.trie.ref_ids[N.ref_begin + k];
        auto it = c.seen->find(ref_id);
        if (it == c.seen->end()) {
            c.seen->emplace(ref_id, static_cast<uint32_t>(c.out->size()));
            c.out->push_back(Hit{ref_id, pen, uint16_t(ns), uint16_t(ni), uint16_t(nd)});
        } else {  // re-reached via a different edit path: a collision
            ++*c.collisions;
            if (pen < (*c.out)[it->second].score)
                (*c.out)[it->second] = Hit{ref_id, pen, uint16_t(ns), uint16_t(ni), uint16_t(nd)};
        }
    }
}

// qpos query chars consumed; (ns,ni,nd) edits so far; pen accumulated penalty.
// Callers guarantee the limits already hold for this state.
void recurse(TmCtx& c, uint32_t node, int qpos, int ns, int ni, int nd, int32_t pen) {
    const Trie::Node& N = c.trie.nodes[node];
    if (qpos == c.qlen && N.ref_count > 0) emit(c, node, pen, ns, ni, nd);

    const int tot = ns + ni + nd;
    for (uint8_t e = 0; e < N.child_count; ++e) {
        uint32_t ei = N.child_begin + e;
        uint8_t code = c.trie.edge_code[ei];
        uint32_t child = c.trie.edge_child[ei];

        if (qpos < c.qlen) {  // substitution / match: consume query + ref char
            bool counts = (c.q[qpos] != code) && !c.pos_masked(qpos);  // masked anchors are free
            int nns = ns + (counts ? 1 : 0);
            int32_t npen = pen + c.sub_pen(qpos, c.q[qpos], code);
            if (nns <= c.lim.max_sub && tot + (counts ? 1 : 0) <= c.lim.max_tot && npen <= c.lim.budget)
                recurse(c, child, qpos + 1, nns, ni, nd, npen);
        }
        if (!c.hamming) {  // insertion: consume ref char only
            int nni = ni + 1;
            int32_t npen = pen + c.lim.gap;
            if (nni <= c.lim.max_ins && tot + 1 <= c.lim.max_tot && npen <= c.lim.budget)
                recurse(c, child, qpos, ns, nni, nd, npen);
        }
    }
    if (!c.hamming && qpos < c.qlen) {  // deletion: consume query char, stay at node
        int nnd = nd + 1;
        int32_t npen = pen + c.lim.gap;
        if (nnd <= c.lim.max_del && tot + 1 <= c.lim.max_tot && npen <= c.lim.budget)
            recurse(c, node, qpos + 1, ns, ni, nnd, npen);
    }
}

}  // namespace

void search_seqtm(const Trie& trie, const uint8_t* qcodes, int qlen, const Limits& lim,
                  Mode mode, Scratch& s, std::vector<Hit>& out) {
    s.seen.clear();
    bool hamming = lim.max_ins == 0 && lim.max_del == 0;
    bool posmode = lim.posmat != nullptr && hamming && lim.posmat->width() == uint16_t(qlen);
    TmCtx c{trie, qcodes, qlen, lim, hamming, posmode, mode,
            &s.seen, &s.collisions, &out};
    recurse(c, /*node=*/0, /*qpos=*/0, 0, 0, 0, 0);
}

}  // namespace seqtree
