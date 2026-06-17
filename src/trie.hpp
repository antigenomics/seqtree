#pragma once
#include "seqtree/seqtree.hpp"

#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace seqtree {

// Arena trie addressed by uint32 node index (no pointers). Children of a node
// occupy a contiguous slice of the parallel edge arrays, sorted by code.
struct Trie {
    struct Node {
        uint32_t child_begin = 0;  // offset into edge_code / edge_child
        uint32_t ref_begin   = 0;  // offset into ref_ids (terminal only)
        uint32_t ref_count   = 0;  // refs ending here; terminal iff > 0
        uint8_t  child_count = 0;
        uint8_t  _pad[3]     = {0, 0, 0};
    };

    std::vector<Node>     nodes;       // node 0 is the root
    std::vector<uint8_t>  edge_code;   // symbol on the edge to a child
    std::vector<uint32_t> edge_child;  // destination node index
    std::vector<uint32_t> ref_ids;     // ref ids grouped per terminal node

    std::string           str_data;    // concatenated reference strings
    std::vector<uint32_t> str_off;     // size n+1; ref i is [str_off[i], str_off[i+1])

    Codec    codec{Alphabet::AminoAcid};
    uint32_t max_depth = 0;

    uint32_t size() const { return static_cast<uint32_t>(str_off.empty() ? 0 : str_off.size() - 1); }
    std::string_view ref_seq(uint32_t id) const {
        return std::string_view(str_data.data() + str_off[id], str_off[id + 1] - str_off[id]);
    }

    // Builds the frozen trie from the given references (moved in). Throws
    // std::invalid_argument on a character outside the alphabet.
    static Trie build(std::vector<std::string> refs, Alphabet);
};

}  // namespace seqtree
