#include "doctest.h"
#include "seqtree/seqtree.hpp"

#include <stdexcept>

using namespace seqtree;

TEST_CASE("build and ref_seq") {
    auto idx = Index::build({"CAT", "CAR", "CST"}, Alphabet::AminoAcid);
    CHECK(idx->size() == 3);
    CHECK(idx->ref_seq(0) == "CAT");
    CHECK(idx->ref_seq(2) == "CST");
}

TEST_CASE("invalid symbol throws") {
    CHECK_THROWS_AS(Index::build({"CA1"}, Alphabet::AminoAcid), std::invalid_argument);
}

TEST_CASE("duplicate refs get distinct ids") {
    auto idx = Index::build({"CAT", "CAT"}, Alphabet::AminoAcid);
    CHECK(idx->size() == 2);
    Searcher s(*idx);
    SearchParams p;  // all caps 0 -> exact match only
    auto hits = s.search("CAT", p);
    CHECK(hits.size() == 2);  // both ids reported
}

TEST_CASE("empty index") {
    auto idx = Index::build({}, Alphabet::AminoAcid);
    CHECK(idx->size() == 0);
    Searcher s(*idx);
    SearchParams p;
    CHECK(s.search("CAT", p).empty());
}
