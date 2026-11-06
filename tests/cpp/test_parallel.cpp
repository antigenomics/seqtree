#include "doctest.h"
#include "seqtree/parallel.hpp"
#include <memory>
#include <chrono>
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
    seqtree::parallel_for(0, 1, [&] { ++factories; return 0; }, [](size_t, int&) {});
    CHECK(factories == 1);
    seqtree::parallel_for(1, 1, [] { return 0; }, [&](size_t, int&) {
        CHECK(std::this_thread::get_id() == caller);
    });
}

TEST_CASE("parallel batches visit each item once and join every scratch") {
    struct Local {
        std::atomic<int>& alive;
        explicit Local(std::atomic<int>& a) : alive(a) { ++alive; }
        ~Local() { --alive; }
    };
    std::atomic<int> alive{0}, factories{0}, caller_factories{0};
    std::vector<std::atomic<int>> visits(1031);
    const auto caller = std::this_thread::get_id();
    seqtree::parallel_for(visits.size(), 2, [&] {
        ++factories;
        if (std::this_thread::get_id() == caller) ++caller_factories;
        return std::make_unique<Local>(alive);
    }, [&](size_t i, const auto&) {
        ++visits[i];
    });
    CHECK(factories == 2);
    CHECK(caller_factories == 1);
    CHECK(alive == 0);
    for (const auto& n : visits) CHECK(n == 1);
}

TEST_CASE("factory and body exceptions propagate after workers join") {
    for (int threads : {1, 2}) {
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

TEST_CASE("parallel batches cap an excessive live worker budget") {
    const unsigned available = seqtree::available_cpus();
    REQUIRE(available >= 1);
    std::atomic<unsigned> factories{0};
    seqtree::parallel_for(size_t(available) + 1, int(available) + 1,
                         [&] { ++factories; return 0; }, [](size_t, int&) {});
    CHECK(factories == std::min(available, 999u));
}


TEST_CASE("automatic batches use physical cores without capping explicit allocation") {
    for (int requested : {0, int(seqtree::available_cpus())}) {
        std::atomic<unsigned> factories{0};
        seqtree::parallel_for(2048, requested, [&] { ++factories; return 0; },
                             [](size_t, int&) {});
        CHECK(factories == (requested ? unsigned(requested) :
                           seqtree::physical_cpus()));
    }
}

TEST_CASE("SMT topology counts only cores represented by affinity") {
    const auto root = std::filesystem::temp_directory_path() /
        ("seqtree-topology-" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
    std::vector<unsigned> cpus;
    for (unsigned cpu = 0; cpu < 40; ++cpu) {
        const auto path = root / ("cpu" + std::to_string(cpu)) / "topology";
        std::filesystem::create_directories(path);
        std::ofstream(path / "physical_package_id") << 0;
        std::ofstream(path / "core_id") << cpu % 20;
        cpus.push_back(cpu);
    }
    CHECK(seqtree::physical_cores_for_cpus(cpus, 40, root.string() + "/") == 20);
    CHECK(seqtree::physical_cores_for_cpus({0, 1, 20, 21}, 4, root.string() + "/") == 2);
    CHECK(seqtree::physical_cores_for_cpus(cpus, 3, root.string() + "/") == 3);
    std::filesystem::remove_all(root);
}

TEST_CASE("Windows hard CPU rates respect affinity and floor fractional CPUs") {
    CHECK(seqtree::cpu_rate_budget(40, 40, 5000) == 20);
    CHECK(seqtree::cpu_rate_budget(4, 40, 5000) == 4);
    CHECK(seqtree::cpu_rate_budget(40, 40, 100) == 1);
    CHECK_THROWS_AS(seqtree::cpu_rate_budget(40, 40, 0), std::invalid_argument);
#ifdef _WIN32
    const auto [logical, physical] = seqtree::windows_cpu_budget();
    CHECK(logical >= physical);
    CHECK(physical >= 1);
#endif
}
