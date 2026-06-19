#include "doctest.h"
#include "seqtree/seqtree.hpp"

#include <map>
#include <vector>

using namespace seqtree;

TEST_CASE("positional matrix: weights, masking, penalty") {
    auto base = SubstitutionMatrix::unit(24);
    auto pm = PositionalMatrix::from_weights(base, {0, 1, 2, 1, 1});  // pos0 masked, pos2 upweighted
    CHECK(pm.width() == 5);
    CHECK(pm.size() == 24);
    CHECK(pm.masked(0));
    CHECK(!pm.masked(1));
    CHECK(pm.penalty(0, 0, 1) == 0);  // masked position: free even on mismatch
    CHECK(pm.penalty(1, 0, 0) == 0);  // match
    CHECK(pm.penalty(1, 0, 1) == 1);  // weight 1 * unit mismatch
    CHECK(pm.penalty(2, 0, 1) == 2);  // weight 2 * unit mismatch
}

TEST_CASE("search with positional matrix: masked anchor is free, others count") {
    // frame width 5; position 0 is a masked "anchor".
    auto idx = Index::build({"AAAAA", "ABAAA", "AABAA"}, Alphabet::AminoAcid);
    auto base = SubstitutionMatrix::unit(24);
    auto pm = PositionalMatrix::from_weights(base, {0, 1, 1, 1, 1});

    SearchParams p;
    p.engine = Engine::SeqTm;
    p.pos_matrix = &pm;
    p.max_substitutions = 0;  // zero TCR-facing (non-masked) mismatches allowed

    Searcher s(*idx);
    std::map<uint32_t, Hit> m;
    for (const auto& h : s.search("CAAAA", p)) m[h.ref_id] = h;  // differs at masked pos 0 only

    CHECK(m.count(0));            // AAAAA: pos0 mismatch is masked -> free
    CHECK(m[0].score == 0);
    CHECK(m[0].n_subs == 0);      // masked mismatch does not count
    CHECK(!m.count(1));           // ABAAA differs at non-masked pos1 -> excluded at max_subs=0
    CHECK(!m.count(2));           // AABAA differs at non-masked pos2 -> excluded
}
