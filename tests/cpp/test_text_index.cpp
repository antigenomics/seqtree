// Completeness is a proof obligation here, not a benchmark. The caller this was built for
// refuses a peptide whose nearest parents name more than one gene, so a MISSED hit does not
// degrade its answer -- it changes it, from ambiguous to confidently wrong. These tests assert
// set equality against an independent brute-force scan, over a grid of query length, max_subs
// and k, rather than merely that hits were found.

#include "doctest.h"
#include "seqtree/text_index.hpp"
#include "tmp_path.hpp"

#include <cstdio>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

using namespace seqtree;

namespace {

std::vector<std::string> random_text(unsigned seed, size_t n_refs, size_t min_len, size_t max_len) {
    std::mt19937 rng(seed);
    // A deliberately small alphabet: it makes near-matches common, so the equivalence tests
    // exercise real candidate lists rather than empty ones.
    const std::string alpha = "ACDEFG";
    std::vector<std::string> refs;
    for (size_t i = 0; i < n_refs; ++i) {
        size_t len = min_len + rng() % (max_len - min_len + 1);
        std::string s;
        for (size_t j = 0; j < len; ++j) s += alpha[rng() % alpha.size()];
        refs.push_back(std::move(s));
    }
    return refs;
}

using Triple = std::tuple<uint32_t, uint32_t, uint16_t>;  // ref_id, offset, n_subs

// Every position in every record, scanned directly. No index, no seeds, no shared code.
std::set<Triple> brute_force(const std::vector<std::string>& refs, const std::string& q,
                             uint16_t max_subs) {
    std::set<Triple> out;
    for (uint32_t r = 0; r < refs.size(); ++r) {
        const std::string& ref = refs[r];
        if (ref.size() < q.size()) continue;
        for (size_t off = 0; off + q.size() <= ref.size(); ++off) {
            uint16_t n = 0;
            for (size_t i = 0; i < q.size() && n <= max_subs; ++i)
                if (ref[off + i] != q[i]) ++n;
            if (n <= max_subs) out.insert({r, uint32_t(off), n});
        }
    }
    return out;
}

std::set<Triple> from_result(const TextResult& res, size_t query) {
    std::set<Triple> out;
    for (uint32_t h = res.query_begin[query]; h < res.query_begin[query + 1]; ++h)
        out.insert({res.ref_id[h], res.offset[h], res.n_subs[h]});
    return out;
}

}  // namespace

TEST_CASE("text index: matches a brute-force scan exactly across L, max_subs and k") {
    auto refs = random_text(20260906, /*n_refs=*/12, /*min_len=*/30, /*max_len=*/90);

    // Queries of three kinds: exact windows of the text, mutated copies, and random strings.
    std::mt19937 rng(7);
    const std::string alpha = "ACDEFG";

    for (uint8_t k : {3, 4, 5}) {
        auto ix = TextIndex::build(refs, Alphabet::AminoAcid, k);
        CHECK(ix->k() == k);
        CHECK(ix->num_refs() == refs.size());

        for (size_t L = 6; L <= 30; ++L) {
            std::vector<std::string> queries;
            for (int rep = 0; rep < 3; ++rep) {  // windows of the text
                const std::string& ref = refs[rng() % refs.size()];
                if (ref.size() < L) continue;
                queries.push_back(ref.substr(rng() % (ref.size() - L + 1), L));
            }
            if (!queries.empty()) {  // a mutated copy of one of them
                std::string mut = queries[0];
                mut[rng() % mut.size()] = alpha[rng() % alpha.size()];
                queries.push_back(mut);
            }
            std::string rand_q;  // pure noise, usually with no hits at all
            for (size_t i = 0; i < L; ++i) rand_q += alpha[rng() % alpha.size()];
            queries.push_back(rand_q);
            if (queries.empty()) continue;

            for (uint16_t m = 0; m <= 3; ++m) {
                TextQueryOpts opts;
                opts.max_subs = m;
                auto res = ix->search_batch(queries, opts, 1);
                for (size_t qi = 0; qi < queries.size(); ++qi) {
                    INFO("k=", int(k), " L=", L, " m=", m, " q=", queries[qi]);
                    CHECK(from_result(res, qi) == brute_force(refs, queries[qi], m));
                }
            }
        }
    }
}

TEST_CASE("text index: both query paths are exercised and agree") {
    // The dispatch is on L / (max_subs + 1) vs k, so at k = 4 an L = 9 query is on the ball
    // path for m >= 2 and on the seed path for m <= 1. Both must give the same answer.
    auto refs = random_text(11, 8, 40, 60);
    auto ix = TextIndex::build(refs, Alphabet::AminoAcid, 4);
    const std::string q = refs[0].substr(5, 9);

    for (uint16_t m = 0; m <= 3; ++m) {
        TextQueryOpts opts;
        opts.max_subs = m;
        auto res = ix->search_batch({q}, opts, 1);
        CHECK(from_result(res, 0) == brute_force(refs, q, m));
    }
}

TEST_CASE("text index: a hit can never span a record boundary") {
    // "AAAA" + "AAAA" as two records. A length-8 all-A query matches inside NEITHER, and must
    // not be allowed to buy the sentinel with one of its two permitted substitutions.
    std::vector<std::string> refs{"AAAA", "AAAA"};
    auto ix = TextIndex::build(refs, Alphabet::AminoAcid, 3);
    TextQueryOpts opts;
    opts.max_subs = 2;
    auto res = ix->search_batch({"AAAAAAAA"}, opts, 1);
    CHECK(res.num_hits() == 0);
}

TEST_CASE("text index: a residue outside the alphabet is a hole, not a wildcard") {
    // U (selenocysteine) is absent from the amino-acid codec -- 36 residues of the human
    // proteome are U -- so it encodes to kInvalid and no hit may cross it.
    // In the TEXT it is kept as a hole and counted -- refusing the build would make the class
    // unusable on the very proteomes it exists for.
    auto ix = TextIndex::build({"AAAUAAA"}, Alphabet::AminoAcid, 3);
    CHECK(ix->num_unknown() == 1);
    TextQueryOpts opts;
    opts.max_subs = 1;
    // No window of AAAUAAA can be reached, because none of them avoids the hole.
    CHECK(ix->search_batch({"AAAAAAA"}, opts, 1).num_hits() == 0);
    CHECK(ix->search_batch({"AAAA"}, opts, 1).num_hits() == 0);
    // ... but a window entirely on one side of it is reachable.
    CHECK(ix->search_batch({"AAA"}, opts, 1).num_hits() == 2);

    // In a QUERY it is refused by name: there is no sensible answer to "match this hole".
    CHECK_THROWS_AS(ix->search_batch({"AAAUAAA"}, opts, 1), std::invalid_argument);
}

TEST_CASE("text index: ambiguity codes match only themselves") {
    // X is a real code in this codec, not a wildcard: XAA is one substitution from AAA.
    // A single length-3 record, so the record has exactly one window to match.
    auto ix = TextIndex::build({"XAA"}, Alphabet::AminoAcid, 3);
    TextQueryOpts exact;
    CHECK(ix->search_batch({"XAA"}, exact, 1).num_hits() == 1);
    CHECK(ix->search_batch({"AAA"}, exact, 1).num_hits() == 0);
    TextQueryOpts one;
    one.max_subs = 1;
    CHECK(ix->search_batch({"AAA"}, one, 1).num_hits() == 1);
}

TEST_CASE("text index: mismatch detail names the pair, not just the position") {
    auto ix = TextIndex::build({"CASSLGQYF"}, Alphabet::AminoAcid, 4);
    TextQueryOpts opts;
    opts.max_subs = 1;
    auto res = ix->search_batch({"CASSPGQYF"}, opts, 1);
    REQUIRE(res.num_hits() == 1);
    REQUIRE(res.n_subs[0] == 1);
    REQUIRE(res.mm_begin[1] - res.mm_begin[0] == 1);
    CHECK(res.mm_pos[0] == 4);
    CHECK(res.mm_query_aa[0] == 'P');
    CHECK(res.mm_text_aa[0] == 'L');
}

TEST_CASE("text index: best_only returns the whole nearest shell") {
    // Three records at distance 0, 1 and 2 from the query. best_only must stop at the nearest
    // shell that has anything in it, and exclude_exact must move it out to the next one.
    std::vector<std::string> refs{"CASSLGQYF", "CASSLGQYW", "CASSLGQWW"};
    auto ix = TextIndex::build(refs, Alphabet::AminoAcid, 4);
    TextQueryOpts opts;
    opts.max_subs = 3;
    opts.best_only = true;
    auto res = ix->search_batch({"CASSLGQYF"}, opts, 1);
    CHECK(res.num_hits() == 1);
    CHECK(res.n_subs[0] == 0);

    opts.exclude_exact = true;
    res = ix->search_batch({"CASSLGQYF"}, opts, 1);
    REQUIRE(res.num_hits() == 1);
    CHECK(res.n_subs[0] == 1);
    CHECK(res.ref_id[0] == 1);
}

TEST_CASE("text index: group fold collapses hits and makes a tie visible") {
    std::vector<std::string> refs{"CASSLGQYF", "CASSLGQYW", "CASSLGQYA"};
    auto ix = TextIndex::build(refs, Alphabet::AminoAcid, 4, {7, 7, 9});
    TextQueryOpts opts;
    opts.max_subs = 1;
    opts.group_by = true;
    auto res = ix->search_batch({"CASSLGQYF"}, opts, 1);
    REQUIRE(res.group_begin.size() == 2);
    REQUIRE(res.group_begin[1] == 2);  // two distinct groups reached
    CHECK(res.group_id[0] == 7);
    CHECK(res.group_min_subs[0] == 0);
    CHECK(res.group_n_hits[0] == 2);
    CHECK(res.group_id[1] == 9);
    CHECK(res.group_min_subs[1] == 1);
}

TEST_CASE("text index: a cap keeps the best hits and always reports itself") {
    std::vector<std::string> refs{"CASSLGQYF", "CASSLGQYW", "CASSLGQYA"};
    auto ix = TextIndex::build(refs, Alphabet::AminoAcid, 4);
    TextQueryOpts opts;
    opts.max_subs = 1;
    auto res = ix->search_batch({"CASSLGQYF"}, opts, 1);
    CHECK(res.num_hits() == 3);
    CHECK(res.truncated[0] == 0);

    opts.max_hits = 2;
    res = ix->search_batch({"CASSLGQYF"}, opts, 1);
    REQUIRE(res.num_hits() == 2);
    CHECK(res.truncated[0] == 1);
    CHECK(res.n_subs[0] == 0);  // sorted before the cap, so the exact hit survives it
}

TEST_CASE("text index: a query shorter than k is refused, not answered incompletely") {
    auto ix = TextIndex::build({"CASSLGQYF"}, Alphabet::AminoAcid, 4);
    TextQueryOpts opts;
    CHECK_THROWS_AS(ix->search_batch({"CAS"}, opts, 1), std::invalid_argument);
}

TEST_CASE("text index: k is capped so a direct-addressed table cannot ask for gigabytes") {
    CHECK_THROWS_AS(TextIndex::build({"CASSLGQYF"}, Alphabet::AminoAcid, 8), std::invalid_argument);
    CHECK_NOTHROW(TextIndex::build({"CASSLGQYF"}, Alphabet::AminoAcid, 5));
}

TEST_CASE("text index: hit order is deterministic across thread counts") {
    auto refs = random_text(99, 40, 40, 80);
    auto ix = TextIndex::build(refs, Alphabet::AminoAcid, 4);
    std::vector<std::string> queries;
    for (const auto& r : refs) queries.push_back(r.substr(0, 12));

    TextQueryOpts opts;
    opts.max_subs = 2;
    auto one = ix->search_batch(queries, opts, 1);
    for (int t : {2, 4, 8}) {
        auto many = ix->search_batch(queries, opts, t);
        CHECK(many.query_begin == one.query_begin);
        CHECK(many.ref_id == one.ref_id);
        CHECK(many.offset == one.offset);
        CHECK(many.n_subs == one.n_subs);
        CHECK(many.mm_pos == one.mm_pos);
    }
}

TEST_CASE("text index: save/load round-trips, both mapped and read") {
    auto refs = random_text(5, 10, 30, 50);
    auto ix = TextIndex::build(refs, Alphabet::AminoAcid, 4, {1, 1, 2, 2, 3, 3, 4, 4, 5, 5});
    const std::string path = tmp_path("seqtree_text_index_test.sti");
    ix->save(path);

    TextQueryOpts opts;
    opts.max_subs = 2;
    opts.group_by = true;
    std::vector<std::string> queries{refs[0].substr(0, 12), refs[3].substr(2, 15)};
    auto want = ix->search_batch(queries, opts, 1);

    for (bool mmap : {false, true}) {
        auto back = TextIndex::load(path, mmap);
        CHECK(back->num_refs() == ix->num_refs());
        CHECK(back->num_residues() == ix->num_residues());
        CHECK(back->k() == ix->k());
        CHECK(back->ref_seq(0) == refs[0]);
        auto got = back->search_batch(queries, opts, 1);
        CHECK(got.query_begin == want.query_begin);
        CHECK(got.ref_id == want.ref_id);
        CHECK(got.offset == want.offset);
        CHECK(got.n_subs == want.n_subs);
        CHECK(got.group_id == want.group_id);
    }
    std::remove(path.c_str());
}

TEST_CASE("text index: a substitution matrix scores a hit without changing which hits return") {
    auto ix = TextIndex::build({"CASSLGQYF"}, Alphabet::AminoAcid, 4);
    auto blosum = SubstitutionMatrix::blosum62();
    TextQueryOpts opts;
    opts.max_subs = 1;
    auto plain = ix->search_batch({"CASSIGQYF"}, opts, 1);  // L -> I, a conservative swap
    opts.matrix = &blosum;
    auto scored = ix->search_batch({"CASSIGQYF"}, opts, 1);
    CHECK(plain.ref_id == scored.ref_id);
    CHECK(plain.offset == scored.offset);
    REQUIRE(scored.num_hits() == 1);
    CHECK(plain.score[0] == 0);
    CHECK(scored.score[0] == blosum.similarity(9, 10));  // I, L in codec order
}

TEST_CASE("text index: nucleotide alphabet works and k may be larger there") {
    std::vector<std::string> refs{"ACGTACGTACGTTTGCA", "TTTTACGTACGTAAACG"};
    auto ix = TextIndex::build(refs, Alphabet::Nucleotide, 8);
    TextQueryOpts opts;
    opts.max_subs = 1;
    auto res = ix->search_batch({"ACGTACGTACGT"}, opts, 1);
    CHECK(from_result(res, 0) == brute_force(refs, "ACGTACGTACGT", 1));
}
