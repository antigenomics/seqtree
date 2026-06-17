#include "engines.hpp"

#include <limits>

namespace seqtree {

namespace {

constexpr int32_t kInf = std::numeric_limits<int32_t>::max() / 4;

struct TrieCtx {
    const Trie& trie;
    const uint8_t* q;
    int qlen;
    int W;                 // qlen + 1
    Limits lim;
    int32_t* buf;          // (max_depth+2) * W scratch rows
    std::vector<Hit>* out;

    int32_t sub_pen(uint8_t a, uint8_t b) const {
        return lim.unit ? (a == b ? 0 : 1) : lim.mat->penalty(a, b);
    }
};

// row at `depth` is already filled. Emit terminal, prune, descend.
void dfs(TrieCtx& c, uint32_t node, int depth) {
    int32_t* row = c.buf + size_t(depth) * c.W;
    const Trie::Node& N = c.trie.nodes[node];

    if (N.ref_count > 0 && row[c.qlen] <= c.lim.budget) {
        for (uint32_t k = 0; k < N.ref_count; ++k)
            c.out->push_back(Hit{c.trie.ref_ids[N.ref_begin + k], row[c.qlen], 0, 0, 0});
    }

    int32_t mn = kInf;  // best achievable below: every edit only adds cost
    for (int j = 0; j <= c.qlen; ++j) mn = mn < row[j] ? mn : row[j];
    if (mn > c.lim.budget) return;

    int32_t* nrow = c.buf + size_t(depth + 1) * c.W;
    const int32_t gap = c.lim.gap;
    for (uint8_t e = 0; e < N.child_count; ++e) {
        uint32_t ei = N.child_begin + e;
        uint8_t code = c.trie.edge_code[ei];
        nrow[0] = row[0] + gap;
        for (int j = 1; j <= c.qlen; ++j) {
            int32_t sub = row[j - 1] + c.sub_pen(c.q[j - 1], code);
            int32_t del = row[j] + gap;       // ref char unmatched (query side gap)
            int32_t ins = nrow[j - 1] + gap;  // query char unmatched
            int32_t v = sub < del ? sub : del;
            nrow[j] = v < ins ? v : ins;
        }
        dfs(c, c.trie.edge_child[ei], depth + 1);
    }
}

}  // namespace

void search_seqtrie(const Trie& trie, const uint8_t* qcodes, int qlen, const Limits& lim,
                    Scratch& s, std::vector<Hit>& out) {
    const int W = qlen + 1;
    s.dp.assign(size_t(trie.max_depth + 2) * W, 0);
    int32_t* root = s.dp.data();
    for (int j = 0; j <= qlen; ++j) root[j] = j * lim.gap;  // delete j query chars vs empty ref

    TrieCtx c{trie, qcodes, qlen, W, lim, s.dp.data(), &out};
    dfs(c, /*node=*/0, /*depth=*/0);
}

}  // namespace seqtree
