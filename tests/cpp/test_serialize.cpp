#include "doctest.h"
#include "seqtree/seqtree.hpp"

#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace seqtree;

TEST_CASE("serialize round-trip preserves search results") {
    std::vector<std::string> refs = {"CASSLAPGATNEKLFF", "CASSLELGATNEKLFF",
                                     "CASSPQGATNEKLFF", "CAT", "CAT"};  // includes a duplicate
    auto idx = Index::build(refs, Alphabet::AminoAcid);
    const char* path = "/tmp/seqtree_serialize_roundtrip.sqtree";
    idx->save(path);
    auto idx2 = Index::load(path);

    CHECK(idx2->size() == idx->size());
    for (uint32_t i = 0; i < idx->size(); ++i) CHECK(idx2->ref_seq(i) == idx->ref_seq(i));

    Searcher s1(*idx), s2(*idx2);
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 2;
    p.max_insertions = 1;
    p.max_deletions = 1;
    for (const auto& q : refs) {
        std::map<uint32_t, int32_t> a, b;
        for (const auto& h : s1.search(q, p)) a[h.ref_id] = h.score;
        for (const auto& h : s2.search(q, p)) b[h.ref_id] = h.score;
        CHECK(a == b);
    }
    std::remove(path);
}

TEST_CASE("load rejects a non-index file") {
    const char* path = "/tmp/seqtree_serialize_bad.bin";
    {
        std::ofstream os(path, std::ios::binary);
        os << "this is definitely not a seqtree index";
    }
    CHECK_THROWS_AS(Index::load(path), std::runtime_error);
    std::remove(path);
}

TEST_CASE("load of a missing file throws") {
    CHECK_THROWS_AS(Index::load("/tmp/seqtree_does_not_exist_12345.sqtree"), std::runtime_error);
}
