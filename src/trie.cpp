#include "trie.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace seqtree {

namespace {
// Mutable build node; children kept sorted by code for deterministic freeze.
struct BuildNode {
    std::vector<std::pair<uint8_t, uint32_t>> children;
    std::vector<uint32_t> refs;
};

uint32_t get_or_add_child(std::vector<BuildNode>& arena, uint32_t parent, uint8_t code) {
    {
        auto& ch = arena[parent].children;
        size_t lo = 0, hi = ch.size();
        while (lo < hi) {  // binary search on sorted (code, idx)
            size_t mid = (lo + hi) / 2;
            if (ch[mid].first < code) lo = mid + 1;
            else hi = mid;
        }
        if (lo < ch.size() && ch[lo].first == code) return ch[lo].second;
    }
    // emplace_back may reallocate arena, so re-find the insert position after.
    uint32_t idx = static_cast<uint32_t>(arena.size());
    arena.emplace_back();
    auto& ch = arena[parent].children;
    size_t lo = 0, hi = ch.size();
    while (lo < hi) {
        size_t mid = (lo + hi) / 2;
        if (ch[mid].first < code) lo = mid + 1;
        else hi = mid;
    }
    ch.insert(ch.begin() + lo, {code, idx});
    return idx;
}
}  // namespace

Trie Trie::build(std::vector<std::string> refs, Alphabet alphabet) {
    Trie t;
    t.codec = Codec(alphabet);
    const Codec& codec = t.codec;

    std::vector<BuildNode> arena;
    arena.emplace_back();  // root

    t.str_off.reserve(refs.size() + 1);
    t.str_off.push_back(0);
    for (auto& s : refs) t.str_data += s;  // store originals; offsets below
    {
        uint32_t off = 0;
        for (auto& s : refs) { off += static_cast<uint32_t>(s.size()); t.str_off.push_back(off); }
    }

    for (uint32_t id = 0; id < refs.size(); ++id) {
        const std::string& s = refs[id];
        uint32_t node = 0;
        for (char ch : s) {
            uint8_t code = codec.encode(ch);
            if (code == Codec::kInvalid)
                throw std::invalid_argument(std::string("invalid symbol '") + ch +
                                            "' in reference #" + std::to_string(id));
            node = get_or_add_child(arena, node, code);
        }
        arena[node].refs.push_back(id);
        if (s.size() > t.max_depth) t.max_depth = static_cast<uint32_t>(s.size());
    }

    // Freeze: final node index == arena index, so child links stay valid.
    const size_t n = arena.size();
    t.nodes.resize(n);
    size_t total_edges = 0, total_refs = 0;
    for (auto& bn : arena) { total_edges += bn.children.size(); total_refs += bn.refs.size(); }
    t.edge_code.reserve(total_edges);
    t.edge_child.reserve(total_edges);
    t.ref_ids.reserve(total_refs);

    for (size_t i = 0; i < n; ++i) {
        BuildNode& bn = arena[i];
        Node& node = t.nodes[i];
        node.child_begin = static_cast<uint32_t>(t.edge_code.size());
        node.child_count = static_cast<uint8_t>(bn.children.size());
        for (auto& [code, child] : bn.children) {
            t.edge_code.push_back(code);
            t.edge_child.push_back(child);
        }
        node.ref_begin = static_cast<uint32_t>(t.ref_ids.size());
        node.ref_count = static_cast<uint32_t>(bn.refs.size());
        for (uint32_t r : bn.refs) t.ref_ids.push_back(r);
    }
    return t;
}

}  // namespace seqtree
