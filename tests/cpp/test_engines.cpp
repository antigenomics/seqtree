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

TEST_CASE("exact match only by default") {
    auto idx = Index::build({"CAT", "CAR", "CST", "WAT"}, Alphabet::AminoAcid);
    auto m = run(*idx, "CAT", SearchParams{});
    CHECK(m.size() == 1);
    CHECK(m.count(0));
    CHECK(m[0].score == 0);
}

TEST_CASE("substitutions (seqtm)") {
    auto idx = Index::build({"CAT", "CAR", "CST", "WST"}, Alphabet::AminoAcid);
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 1;
    auto m = run(*idx, "CAT", p);
    // within 1 sub of CAT: CAT(0), CAR(T->R), CST(A->S). WST is 2 subs (C->W, A->S).
    CHECK(m.count(0));
    CHECK(m.count(1));
    CHECK(m.count(2));
    CHECK(!m.count(3));
    CHECK(m[1].n_subs == 1);
    CHECK(m[1].score == 1);
}

TEST_CASE("per-type deletion vs insertion caps (seqtm)") {
    // ids: CASS=0, CAS=1, CASSS=2. Query "CASS".
    auto idx = Index::build({"CASS", "CAS", "CASSS"}, Alphabet::AminoAcid);

    SearchParams del;
    del.engine = Engine::SeqTm;
    del.max_deletions = 1;
    auto md = run(*idx, "CASS", del);
    CHECK(md.count(0));         // exact
    CHECK(md.count(1));         // CAS: drop one query char
    CHECK(!md.count(2));        // CASSS needs an insertion
    CHECK(md[1].n_dels == 1);

    SearchParams ins;
    ins.engine = Engine::SeqTm;
    ins.max_insertions = 1;
    auto mi = run(*idx, "CASS", ins);
    CHECK(mi.count(0));
    CHECK(mi.count(2));         // CASSS: add one char
    CHECK(!mi.count(1));        // CAS needs a deletion
    CHECK(mi[2].n_ins == 1);
}

TEST_CASE("matrix-weighted budget (seqtrie)") {
    auto blosum = SubstitutionMatrix::blosum62();
    auto idx = Index::build({"CAT", "CAR", "WAT"}, Alphabet::AminoAcid);
    SearchParams p;
    p.engine = Engine::SeqTrie;
    p.matrix = &blosum;
    p.gap_open = 100;        // make indels irrelevant
    p.max_score_penalty = 6;
    auto m = run(*idx, "CAT", p);
    // CAT=0; CAR: T->R penalty 6; WAT: C->W penalty 13 (excluded)
    CHECK(m.count(0));
    CHECK(m[0].score == 0);
    CHECK(m.count(1));
    CHECK(m[1].score == 6);
    CHECK(!m.count(2));
}

TEST_CASE("top hit") {
    auto idx = Index::build({"CAT", "CAR", "WAT"}, Alphabet::AminoAcid);
    Searcher s(*idx);
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 2;
    Hit top;
    CHECK(s.search_top("CAT", p, top));
    CHECK(top.ref_id == 0);
    CHECK(top.score == 0);
}

TEST_CASE("seqtrie driven by max_total_edits alone") {
    // Per-type caps default to 0; seqtrie ignores them and uses the total budget.
    auto idx = Index::build({"CASSLAPGATNEKLFF"}, Alphabet::AminoAcid);
    SearchParams p;
    p.engine = Engine::SeqTrie;
    p.max_total_edits = 1;
    auto m = run(*idx, "CASSLWPGATNEKLFF", p);  // one substitution
    CHECK(m.count(0));
    CHECK(m[0].score == 1);
}

TEST_CASE("engine agreement within total-edit budget") {
    std::vector<std::string> refs = {
        "CASSLAPGATNEKLFF", "CASSLELGATNEKLFF", "CASSPQGATNEKLFF",
        "CASSLAPGATNEKLF",  "CSSSLAPGATNEKLFF", "WASSLAPGATNEKLFF"};
    auto idx = Index::build(refs, Alphabet::AminoAcid);

    SearchParams pa;
    pa.engine = Engine::SeqTm;
    pa.max_substitutions = 2;
    pa.max_insertions = 2;
    pa.max_deletions = 2;
    pa.max_total_edits = 2;
    SearchParams pb = pa;
    pb.engine = Engine::SeqTrie;

    auto ma = run(*idx, "CASSLAPGATNEKLFF", pa);
    auto mb = run(*idx, "CASSLAPGATNEKLFF", pb);
    CHECK(ma.size() == mb.size());
    for (const auto& [id, h] : ma) {
        CHECK(mb.count(id));
        CHECK(mb[id].score == h.score);
    }
}

TEST_CASE("alignment ops match reported edit counts") {
    auto idx = Index::build({"CASSLAPGATNEKLFF"}, Alphabet::AminoAcid);
    Searcher s(*idx);
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 2;
    std::string q = "CASSLELGATNEKLFF";  // 2 substitutions vs ref
    auto hits = s.search(q, p);
    REQUIRE(hits.size() == 1);
    CHECK(hits[0].n_subs == 2);

    auto al = idx->align(q, 0, p);
    int S = 0, I = 0, D = 0;
    for (char c : al.ops) {
        if (c == 'S') ++S;
        else if (c == 'I') ++I;
        else if (c == 'D') ++D;
    }
    CHECK(S == hits[0].n_subs);
    CHECK(I == hits[0].n_ins);
    CHECK(D == hits[0].n_dels);
    CHECK(al.score == hits[0].score);
}
