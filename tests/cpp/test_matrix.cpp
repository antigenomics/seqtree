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
    // A=0, R=1: sim(A,R)=-1, max(diag A=4, R=5)=5 -> penalty 6
    CHECK(m.penalty(0, 1) == 6);
    // W=17, Y=18: sim(W,Y)=2, max(diag W=11, Y=7)=11 -> penalty 9
    CHECK(m.penalty(17, 18) == 9);
    CHECK(m.penalty(0, 1) == m.penalty(1, 0));  // symmetric
}
