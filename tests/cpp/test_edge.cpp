#include "doctest.h"
#include "count_unit.hpp"
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

TEST_CASE("grouped gapblock population buffers and full-key puncture") {
    seqtree::GroupedGeometry geometry;
    seqtree::GroupedCountBatch batch({"AAA","AAC"},{"AAA","ACA"},{0,2},
        {0,2,3},{0,1,1},{0,1},{-1,1},{0,1},{1,0},1,
        seqtree::Alphabet::AminoAcid,nullptr,1,1,{},0,1,true,geometry);
    auto first=batch.next();
    REQUIRE(first.has_value());
    CHECK(first->start==0);
    CHECK(first->stop==1);
    CHECK(first->counts.shape==std::vector<size_t>{1,7,2});
    CHECK(first->counts.data[0]==0);
    CHECK(first->counts.data[1]==1);
    CHECK(first->sizes.data[0]==1);
    CHECK(first->exact_counts.data==std::vector<uint32_t>{1,1});
    REQUIRE(batch.next().has_value());
    CHECK_FALSE(batch.next().has_value());
    CHECK_THROWS_AS(seqtree::GroupedCountBatch({"AAA"},{"AAA"},{0,0},
        {0,1},{0},{0},{0},{0},{0},1,seqtree::Alphabet::AminoAcid,nullptr,1,1,{},0,1,false,geometry),std::invalid_argument);
}

// The bounded Myers rejection must preserve every within-cap distance.
TEST_CASE("bounded unit bitvector distance preserves all short binary distances") {
    Codec codec(Alphabet::AminoAcid);
    std::vector<std::string> strings{ "" };
    for (unsigned n=1;n<=6;++n) for(unsigned mask=0;mask<(1u<<n);++mask) {
        std::string s(n,'A');
        for(unsigned i=0;i<n;++i) if(mask&(1u<<i)) s[i]='C';
        strings.push_back(s);
    }
    for (const auto& q:strings) {
        if(q.empty()) continue;
        count_unit::Pattern pattern(q,codec);
        for(const auto& r:strings) {
            std::vector<int> row(q.size()+1);
            for(size_t i=0;i<row.size();++i) row[i]=int(i);
            for(size_t j=0;j<r.size();++j) {
                int diagonal=row[0]; row[0]=int(j+1);
                for(size_t i=1;i<row.size();++i) {
                    int old=row[i];
                    row[i]=std::min({old+1,row[i-1]+1,diagonal+int(q[i-1]!=r[j])});
                    diagonal=old;
                }
            }
            for(int cap=1;cap<=5;++cap)
                CHECK(pattern.distance(r,cap)==std::min(row.back(),cap+1));
        }
    }
    count_unit::Pattern long_pattern(std::string(63,'A'),codec);
    CHECK(long_pattern.distance(std::string(63,'A'),3)==0);
    CHECK(long_pattern.distance(std::string(63,'C'),3)==4);
    CHECK(long_pattern.distance(std::string(65,'A'),3)==2);
}
