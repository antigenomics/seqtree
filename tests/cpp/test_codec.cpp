#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "seqtree/seqtree.hpp"

using namespace seqtree;

TEST_CASE("amino-acid codec") {
    Codec c(Alphabet::AminoAcid);
    CHECK(c.size() == 24);
    CHECK(c.encode('A') == 0);
    CHECK(c.encode('a') == 0);  // case-insensitive
    CHECK(c.decode(0) == 'A');
    CHECK(c.encode('R') == 1);
    CHECK(c.encode('*') == 23);
    CHECK(c.encode('1') == Codec::kInvalid);
    CHECK(c.encode('J') == Codec::kInvalid);  // not in BLOSUM62 order
}

TEST_CASE("nucleotide codecs") {
    Codec c(Alphabet::Nucleotide);
    CHECK(c.size() == 4);
    CHECK(c.encode('A') == 0);
    CHECK(c.encode('T') == 3);
    CHECK(c.encode('N') == Codec::kInvalid);

    Codec iu(Alphabet::NucleotideIUPAC);
    CHECK(iu.encode('N') != Codec::kInvalid);
}
