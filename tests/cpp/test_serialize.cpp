#include "doctest.h"
#include "seqtree/seqtree.hpp"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace seqtree;

TEST_CASE("serialize round-trip preserves search results") {
    std::vector<std::string> refs = {"CASSLAPGATNEKLFF", "CASSLELGATNEKLFF",
                                     "CASSPQGATNEKLFF", "CAT", "CAT"};  // includes a duplicate
    auto idx = Index::build(refs, Alphabet::AminoAcid);
    const char* path = "/tmp/seqtree_serialize_roundtrip.sqtree";
    idx->save(path);
    auto idx2 = Index::load(path);

    CHECK(idx2->size() == idx->size());
    for (uint32_t i = 0; i < idx->size(); ++i) CHECK(idx2->ref_seq(i) == idx->ref_seq(i));

    Searcher s1(*idx), s2(*idx2);
    SearchParams p;
    p.engine = Engine::SeqTm;
    p.max_substitutions = 2;
    p.max_insertions = 1;
    p.max_deletions = 1;
    for (const auto& q : refs) {
        std::map<uint32_t, int32_t> a, b;
        for (const auto& h : s1.search(q, p)) a[h.ref_id] = h.score;
        for (const auto& h : s2.search(q, p)) b[h.ref_id] = h.score;
        CHECK(a == b);
    }
    std::remove(path);
}

TEST_CASE("load rejects a non-index file") {
    const char* path = "/tmp/seqtree_serialize_bad.bin";
    {
        std::ofstream os(path, std::ios::binary);
        os << "this is definitely not a seqtree index";
    }
    CHECK_THROWS_AS(Index::load(path), std::runtime_error);
    std::remove(path);
}

TEST_CASE("load of a missing file throws") {
    CHECK_THROWS_AS(Index::load("/tmp/seqtree_does_not_exist_12345.sqtree"), std::runtime_error);
}

// Concurrency: save writes through a temporary and renames it into place, so a reader never sees
// a partial index. These guard the invariants that make that true.

TEST_CASE("save leaves no temporary behind, and overwrites cleanly") {
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "seqtree_atomic_test";
    fs::remove_all(dir);
    fs::create_directories(dir);
    const std::string path = (dir / "idx.sqtree").string();

    auto idx = Index::build({"CASSLGQAYEQYF", "CASSPGQAYEQYF"}, Alphabet::AminoAcid);
    idx->save(path);
    idx->save(path);  // overwrite an existing index

    size_t files = 0, strays = 0;
    for (const auto& e : fs::directory_iterator(dir)) {
        ++files;
        if (e.path().string().find(".tmp.") != std::string::npos) ++strays;
    }
    CHECK(files == 1);
    CHECK(strays == 0);
    CHECK(Index::load(path)->size() == 2);
    fs::remove_all(dir);
}

TEST_CASE("a failed save cleans up its temporary and does not clobber the old index") {
    namespace fs = std::filesystem;
    const fs::path dir = fs::temp_directory_path() / "seqtree_atomic_fail";
    fs::remove_all(dir);
    fs::create_directories(dir);

    // A good index is already in place.
    const std::string good = (dir / "idx.sqtree").string();
    Index::build({"CASSLGQAYEQYF"}, Alphabet::AminoAcid)->save(good);

    // Saving into a directory that does not exist must throw, and litter nothing.
    const std::string bad = (dir / "missing" / "idx.sqtree").string();
    CHECK_THROWS(Index::build({"CASSLGQAYEQYF"}, Alphabet::AminoAcid)->save(bad));

    size_t strays = 0;
    for (const auto& e : fs::directory_iterator(dir))
        if (e.path().string().find(".tmp.") != std::string::npos) ++strays;
    CHECK(strays == 0);
    CHECK(Index::load(good)->size() == 1);  // the pre-existing index survived
    fs::remove_all(dir);
}
