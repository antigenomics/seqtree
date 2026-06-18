#include "doctest.h"
#include "seqtree/seqtree.hpp"

using namespace seqtree;

TEST_CASE("unit matrix") {
    auto m = SubstitutionMatrix::unit(4);
    CHECK(m.penalty(0, 0) == 0);
    CHECK(m.penalty(0, 1) == 1);
    CHECK(m.penalty(3, 3) == 0);
}

TEST_CASE("blosum62 penalties non-negative, identity zero") {
    auto m = SubstitutionMatrix::blosum62();
    CHECK(m.size() == 24);
    for (uint8_t i = 0; i < 24; ++i) CHECK(m.penalty(i, i) == 0);
    // squared-distance: pen(a,b) = s_aa + s_bb - 2 s_ab
    // A=0, R=1: 4 + 5 - 2*(-1) = 11
    CHECK(m.penalty(0, 1) == 11);
    // W=17, Y=18: 11 + 7 - 2*2 = 14
    CHECK(m.penalty(17, 18) == 14);
    CHECK(m.penalty(0, 1) == m.penalty(1, 0));  // symmetric
}
