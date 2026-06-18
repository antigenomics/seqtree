#include "doctest.h"
#include "seqtree/seqtree.hpp"

#include <map>
#include <string>
#include <vector>

using namespace seqtree;

static std::map<uint32_t, Hit> run(Index& idx, const std::string& q, SearchParams p) {
    Searcher s(idx);
    std::map<uint32_t, Hit> m;
    for (const auto& h : s.search(q, p)) m[h.ref_id] = h;
    return m;
}

TEST_CASE("very long sequences: exact, one substitution, both engines") {
    std::string r;
    for (int i = 0; i < 2000; ++i) r += "ACDEFGHIKLMNPQRSTVWY"[i % 20];
    std::string homo(2000, 'A');
    auto idx = Index::build({r, homo}, Alphabet::AminoAcid);
    CHECK(idx->size() == 2);

    CHECK(run(*idx, r, SearchParams{}).count(0));  // exact

    std::string rm = r;
    rm[1000] = (rm[1000] == 'A') ? 'C' : 'A';  // one substitution mid-sequence
    SearchParams tm;
    tm.engine = Engine::SeqTm;
    tm.max_substitutions = 1;
    auto m = run(*idx, rm, tm);
    CHECK(m.count(0));
    CHECK(m[0].n_subs == 1);

    SearchParams tr;
    tr.engine = Engine::SeqTrie;
    tr.max_total_edits = 1;
    CHECK(run(*idx, rm, tr).count(0));
}

TEST_CASE("empty reference and empty query") {
    auto idx = Index::build({"", "A", "AA"}, Alphabet::AminoAcid);
    CHECK(idx->size() == 3);
    CHECK(idx->ref_seq(0).empty());

    CHECK(run(*idx, "", SearchParams{}).count(0));  // empty query == empty ref

    SearchParams ins;
    ins.engine = Engine::SeqTm;
    ins.max_insertions = 2;
    auto mi = run(*idx, "", ins);  // reach A, AA by insertions
    CHECK(mi.count(0));
    CHECK(mi.count(1));
    CHECK(mi.count(2));

    CHECK(run(*idx, "A", SearchParams{}).count(1));  // exact against an index holding ""
}

TEST_CASE("query much longer than references (deletions)") {
    auto idx = Index::build({"CAT"}, Alphabet::AminoAcid);
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_deletions = 3;
    auto m = run(*idx, "CATWYF", p);  // three extra query chars -> three deletions
    CHECK(m.count(0));
    CHECK(m[0].n_dels == 3);
}

TEST_CASE("homopolymer indels reach neighbors of different length") {
    auto idx = Index::build({"AAAAAAAA", "AAAAAAA", "AAAAAAAAA"}, Alphabet::AminoAcid);  // 8,7,9
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_insertions = 1;
    p.max_deletions = 1;
    auto m = run(*idx, "AAAAAAAA", p);
    CHECK(m.count(0));  // exact
    CHECK(m.count(1));  // shorter, via deletion
    CHECK(m.count(2));  // longer, via insertion
}

TEST_CASE("single-character references") {
    auto idx = Index::build({"A", "C", "W"}, Alphabet::AminoAcid);
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 1;
    auto m = run(*idx, "A", p);
    CHECK(m.size() == 3);  // every length-1 ref is within one substitution
}
