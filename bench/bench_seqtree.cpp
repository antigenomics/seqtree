// Dependency-free chrono benchmark for the seqtree core.
//   ./seqtree_bench [size ...]      default sizes: 1000 10000
// Reports (TSV): build time + peak RSS, single-query latency, batch throughput,
// and thread scaling, on seeded synthetic amino-acid data.
#include "seqtree/seqtree.hpp"

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#ifdef __APPLE__
#include <mach/mach.h>
static size_t peak_rss_mb() {
    mach_task_basic_info info;
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO, (task_info_t)&info, &count) == KERN_SUCCESS)
        return info.resident_size_max / (1024 * 1024);
    return 0;
}
#else
static size_t peak_rss_mb() { return 0; }
#endif

using clk = std::chrono::steady_clock;
static double ms(clk::duration d) { return std::chrono::duration<double, std::milli>(d).count(); }

static const char* kAA = "ACDEFGHIKLMNPQRSTVWY";

static std::vector<std::string> random_db(size_t n, std::mt19937_64& rng, int lo = 12, int hi = 18) {
    std::uniform_int_distribution<int> len(lo, hi), sym(0, 19);
    std::vector<std::string> db(n);
    for (size_t i = 0; i < n; ++i) {
        int L = len(rng);
        db[i].resize(L);
        for (int j = 0; j < L; ++j) db[i][j] = kAA[sym(rng)];
    }
    return db;
}

// Mutate db members with `nsubs` substitutions to make queries with known answers.
static std::vector<std::string> mutated_queries(const std::vector<std::string>& db, size_t n,
                                                 int nsubs, std::mt19937_64& rng) {
    std::uniform_int_distribution<size_t> pick(0, db.size() - 1);
    std::uniform_int_distribution<int> sym(0, 19);
    std::vector<std::string> q(n);
    for (size_t i = 0; i < n; ++i) {
        std::string s = db[pick(rng)];
        for (int k = 0; k < nsubs && !s.empty(); ++k) {
            std::uniform_int_distribution<size_t> pos(0, s.size() - 1);
            s[pos(rng)] = kAA[sym(rng)];
        }
        q[i] = s;
    }
    return q;
}

int main(int argc, char** argv) {
    std::vector<size_t> sizes;
    for (int i = 1; i < argc; ++i) sizes.push_back(std::stoul(argv[i]));
    if (sizes.empty()) sizes = {1000, 10000};

    std::mt19937_64 rng(42);
    std::printf("size\tbuild_ms\trss_mb\tq_med_us\tq_p99_us\tbatch_qps\tt1_qps\tt4_qps\tt8_qps\tt16_qps\n");

    std::string methods;  // second table: two engines + alignment cost, printed after

    for (size_t N : sizes) {
        auto db = random_db(N, rng);
        size_t nq = std::min<size_t>(N, 20000);
        auto queries = mutated_queries(db, nq, 1, rng);

        auto t0 = clk::now();
        auto idx = seqtree::Index::build(db, seqtree::Alphabet::AminoAcid);
        double build_ms = ms(clk::now() - t0);

        seqtree::SearchParams p;
        p.engine = seqtree::Engine::SeqTm;
        p.max_substitutions = 1;

        // single-query latency
        std::vector<double> lat;
        lat.reserve(std::min<size_t>(nq, 5000));
        {
            seqtree::Searcher s(*idx);
            std::vector<seqtree::Hit> out;
            for (size_t i = 0; i < lat.capacity(); ++i) {
                auto a = clk::now();
                s.search_into(queries[i % queries.size()], p, out);
                lat.push_back(ms(clk::now() - a) * 1000.0);  // us
            }
        }
        std::sort(lat.begin(), lat.end());
        double med = lat[lat.size() / 2];
        double p99 = lat[std::min(lat.size() - 1, size_t(lat.size() * 0.99))];

        auto qps_at = [&](int threads) {
            auto a = clk::now();
            auto res = idx->search_batch(queries, p, threads);
            double sec = ms(clk::now() - a) / 1000.0;
            (void)res;
            return queries.size() / sec;
        };
        double batch_qps = qps_at(0);
        double t1 = qps_at(1), t4 = qps_at(4), t8 = qps_at(8), t16 = qps_at(16);

        std::printf("%zu\t%.1f\t%zu\t%.2f\t%.2f\t%.0f\t%.0f\t%.0f\t%.0f\t%.0f\n",
                    N, build_ms, peak_rss_mb(), med, p99, batch_qps, t1, t4, t8, t16);
        std::fflush(stdout);

        // --- two-engine comparison across scopes + alignment cost ---
        char buf[256];
        for (int scope = 1; scope <= 2; ++scope) {
            auto sq = mutated_queries(db, nq, scope, rng);
            for (auto eng : {seqtree::Engine::SeqTm, seqtree::Engine::SeqTrie}) {
                seqtree::SearchParams mp;
                mp.engine = eng;
                mp.max_substitutions = scope;
                mp.max_total_edits = scope;
                auto a = clk::now();
                auto res = idx->search_batch(sq, mp, 0);
                double q = sq.size() / (ms(clk::now() - a) / 1000.0);
                size_t hits = 0;
                for (auto& r : res) hits += r.size();
                std::snprintf(buf, sizeof buf, "%zu\t%s\tmax_edits=%d\t%.0f\t%.2f\n", N,
                              eng == seqtree::Engine::SeqTm ? "seqtm" : "seqtrie", scope, q,
                              double(hits) / sq.size());
                methods += buf;
            }
        }
        // alignment fetch cost: align each query against its top hit
        {
            seqtree::SearchParams ap;
            ap.engine = seqtree::Engine::SeqTm;
            ap.max_substitutions = 1;
            auto aq = mutated_queries(db, std::min<size_t>(nq, 5000), 1, rng);
            seqtree::Searcher s(*idx);
            std::vector<std::pair<uint32_t, std::string>> pairs;
            for (auto& q : aq) {
                seqtree::Hit h;
                if (s.search_top(q, ap, h)) pairs.emplace_back(h.ref_id, q);
            }
            auto a = clk::now();
            for (auto& [rid, q] : pairs) idx->align(q, rid, ap);
            double per = pairs.empty() ? 0 : ms(clk::now() - a) * 1000.0 / pairs.size();
            std::snprintf(buf, sizeof buf, "%zu\talign\tper_call_us\t%.2f\t%zu\n", N, per,
                          pairs.size());
            methods += buf;
        }
    }

    std::printf("\nsize\tengine\tknob\tqps\thits_per_q\n%s", methods.c_str());
    return 0;
}
