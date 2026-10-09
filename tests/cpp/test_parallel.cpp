#include "doctest.h"
#include "seqtree/parallel.hpp"
#include <memory>
#include <stdexcept>

TEST_CASE("serial batches run inline and create one scratch") {
    const auto caller = std::this_thread::get_id();
    size_t factories = 0, seen = 0;
    seqtree::parallel_for(17, 1, [&] { ++factories; return size_t{0}; },
                         [&](size_t i, size_t& scratch) {
        CHECK(std::this_thread::get_id() == caller);
        CHECK(i == seen++);
        CHECK(scratch++ == i);
    });
    CHECK(factories == 1);
    CHECK(seen == 17);
    seqtree::parallel_for(0, 8, [&] { ++factories; return 0; }, [](size_t, int&) {});
    CHECK(factories == 1);
    seqtree::parallel_for(1, 8, [] { return 0; }, [&](size_t, int&) {
        CHECK(std::this_thread::get_id() == caller);
    });
}

TEST_CASE("parallel batches visit each item once and join every scratch") {
    struct Local {
        std::atomic<int>& alive;
        explicit Local(std::atomic<int>& a) : alive(a) { ++alive; }
        ~Local() { --alive; }
    };
    std::atomic<int> alive{0}, factories{0};
    std::vector<std::atomic<int>> visits(1031);
    const auto caller = std::this_thread::get_id();
    seqtree::parallel_for(visits.size(), 4, [&] {
        ++factories;
        return std::make_unique<Local>(alive);
    }, [&](size_t i, const auto&) {
        CHECK(std::this_thread::get_id() != caller);
        ++visits[i];
    });
    CHECK(factories == 4);
    CHECK(alive == 0);
    for (const auto& n : visits) CHECK(n == 1);
}

TEST_CASE("factory and body exceptions propagate after workers join") {
    for (int threads : {1, 4}) {
        CHECK_THROWS_WITH_AS(seqtree::parallel_for(100, threads, []() -> int {
            throw std::runtime_error("scratch failed");
        }, [](size_t, int&) {}), "scratch failed", std::runtime_error);
        std::atomic<int> alive{0};
        auto make = [&] {
            ++alive;
            return std::shared_ptr<int>(new int, [&](int* p) { delete p; --alive; });
        };
        CHECK_THROWS_WITH_AS(seqtree::parallel_for(100, threads, make,
                             [](size_t, const auto&) { throw std::runtime_error("body failed"); }),
                             "body failed", std::runtime_error);
        CHECK(alive == 0);
    }
}
