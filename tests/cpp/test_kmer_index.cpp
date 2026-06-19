#include "doctest.h"
#include "seqtree/kmer_index.hpp"

#include <cstdio>
#include <map>
#include <string>
#include <vector>

using namespace seqtree;

// Three "peptides" via their k-mer sets; peptides 0 and 1 share two k-mers.
static std::vector<std::vector<std::string>> kmers() {
    return {
        {"CASS", "ASSL", "SSLA"},  // peptide 0
        {"CASS", "ASSL", "SSLE"},  // peptide 1 (shares CASS, ASSL with 0)
        {"WXYZ", "XYZW"},          // peptide 2 (disjoint)  (X is a valid AA symbol)
    };
}

TEST_CASE("kmer index: build, seed_and_gather, shared counts") {
    auto ki = KmerIndex::build(kmers(), Alphabet::AminoAcid, {10, 20, 10});
    CHECK(ki->num_peptides() == 3);
    CHECK(ki->num_kmers() == 6);  // distinct: CASS ASSL SSLA SSLE WXYZ XYZW

    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 0;  // exact k-mer seeding
    auto res = ki->seed_and_gather({{"CASS", "ASSL"}}, p, /*min_shared=*/1, /*allele_filter=*/-1, 1);
    REQUIRE(res.size() == 1);
    std::map<uint32_t, uint32_t> shared;
    for (const auto& c : res[0]) shared[c.peptide_id] = c.shared_kmers;
    CHECK(shared[0] == 2);   // peptide 0 shares CASS + ASSL
    CHECK(shared[1] == 2);   // peptide 1 shares CASS + ASSL
    CHECK(shared.count(2) == 0);
    CHECK(res[0].front().shared_kmers == 2);  // ranked: top has the most shared
}

TEST_CASE("kmer index: allele filter") {
    auto ki = KmerIndex::build(kmers(), Alphabet::AminoAcid, {10, 20, 10});
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 0;
    // restrict to allele 20 -> only peptide 1 is eligible
    auto res = ki->seed_and_gather({{"CASS", "ASSL"}}, p, 1, /*allele_filter=*/20, 1);
    REQUIRE(res.size() == 1);
    CHECK(res[0].size() == 1);
    CHECK(res[0][0].peptide_id == 1);
}

TEST_CASE("kmer index: fuzzy seeding with one substitution") {
    auto ki = KmerIndex::build(kmers(), Alphabet::AminoAcid, {});
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 1;
    // "SSLD" is 1 sub from SSLA (pep0) and SSLE (pep1)
    auto res = ki->seed_and_gather({{"SSLD"}}, p, 1, -1, 1);
    std::map<uint32_t, uint32_t> shared;
    for (const auto& c : res[0]) shared[c.peptide_id] = c.shared_kmers;
    CHECK(shared.count(0));
    CHECK(shared.count(1));
}

TEST_CASE("kmer index: save/load round-trip") {
    auto ki = KmerIndex::build(kmers(), Alphabet::AminoAcid, {10, 20, 10});
    const char* path = "/tmp/seqtree_kmer_roundtrip.sqkm";
    ki->save(path);
    auto ki2 = KmerIndex::load(path);
    CHECK(ki2->num_peptides() == ki->num_peptides());
    CHECK(ki2->num_kmers() == ki->num_kmers());
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 0;
    auto a = ki->seed_and_gather({{"CASS"}}, p, 1, -1, 1);
    auto b = ki2->seed_and_gather({{"CASS"}}, p, 1, -1, 1);
    CHECK(a[0].size() == b[0].size());
    std::remove(path);
}
