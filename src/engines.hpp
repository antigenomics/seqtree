#pragma once
#include "seqtree/types.hpp"
#include "trie.hpp"

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace seqtree {

// Reusable per-thread scratch shared by both engines.
struct Scratch {
    std::vector<int32_t>  dp;    // seqtrie: (max_depth+2) * (qlen+1) row buffer
    std::unordered_map<uint32_t, uint32_t> seen;  // seqtm: ref_id -> index in out
    uint64_t collisions = 0;     // seqtm: times a ref was re-reached via a different edit path
};

// Resolved, engine-independent limits derived from SearchParams + matrix.
struct Limits {
    int32_t budget;        // max accumulated penalty
    int     max_sub, max_ins, max_del, max_tot;
    int32_t gap;           // linear gap cost per indel
    bool    unit;          // true => unit cost, false => use matrix
    const SubstitutionMatrix* mat;
    const PositionalMatrix* posmat;  // optional per-position penalties (Hamming path only)
    uint32_t max_hits;     // 0 => unlimited
};

Limits resolve_limits(const SearchParams&);

// Both append to out (caller clears). qcodes has length qlen, already encoded.
void search_seqtm(const Trie&, const uint8_t* qcodes, int qlen, const Limits&,
                  Scratch&, std::vector<Hit>& out);
void search_seqtrie(const Trie&, const uint8_t* qcodes, int qlen, const Limits&,
                    Scratch&, std::vector<Hit>& out);

}  // namespace seqtree
