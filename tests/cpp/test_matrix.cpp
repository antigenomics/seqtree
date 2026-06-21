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

TEST_CASE("built-in amino-acid matrices: size, zero diagonal, symmetry") {
    for (auto m : {SubstitutionMatrix::pam250(), SubstitutionMatrix::pam100(),
                   SubstitutionMatrix::structural()}) {
        CHECK(m.size() == 24);
        for (uint8_t i = 0; i < 24; ++i) {
            CHECK(m.penalty(i, i) == 0);
            for (uint8_t j = 0; j < 24; ++j) CHECK(m.penalty(i, j) == m.penalty(j, i));
        }
    }
    // A=0, E=6 (squared distance). PAM250: 2+4-2*0=6; PAM100: 4+5-2*0=9.
    CHECK(SubstitutionMatrix::pam250().penalty(0, 6) == 6);
    CHECK(SubstitutionMatrix::pam100().penalty(0, 6) == 9);
    // structural is on a 0..10 similarity scale: A-A=10, E-E=10, A-E=6 -> 10+10-12=8.
    CHECK(SubstitutionMatrix::structural().penalty(0, 6) == 8);
}
