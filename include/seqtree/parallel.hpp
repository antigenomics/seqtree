#pragma once
// One parallel-for scaffold for the whole tree. Every batch entry point in seqtree used to
// carry its own copy of this: an atomic cursor, an adaptive chunk, an exception_ptr, a vector
// of threads, and a rethrow after the join. Seven copies drifted -- one of them had no
// exception plumbing at all.
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <mutex>
#include <thread>
#include <vector>

namespace seqtree {

// Runs body(i, local) for every i in [0, n), across nt workers pulling contiguous chunks off
// one atomic cursor. `threads <= 0` means hardware concurrency.
//
// `make_local()` is called once per worker and its result handed to every body() that worker
// runs -- that is how a Searcher or a scratch buffer gets allocated once per thread rather
// than once per item. Return a dummy for kernels that need no state.
//
// The chunk is adaptive: large batches keep ~1024 items per grab so the atomic is not the
// bottleneck, small ones split ~8 chunks per thread so no worker goes idle. Raise `chunk_min`
// when the per-item output struct is small enough that neighbouring workers would share a
// cache line (TextIndex's QueryOut is ~80 bytes, so it passes 8).
//
// Workers catch, the first exception is kept, and it is rethrown after the join. This is not
// optional: an exception escaping a std::thread entry function calls std::terminate, which
// from Python is an uncatchable SIGABRT rather than a raise.
template <class MakeLocal, class Body>
void parallel_for(size_t n, int threads, MakeLocal make_local, Body body,
                  size_t chunk_min = 1, size_t chunk_max = 1024) {
    if (n == 0) return;

    unsigned nt = threads > 0 ? unsigned(threads)
                              : std::max(1u, std::thread::hardware_concurrency());
    nt = unsigned(std::min<size_t>(nt, n));  // n >= 1 here, so nt >= 1

    std::atomic<size_t> next{0};
    const size_t chunk = std::clamp<size_t>(n / (size_t(nt) * 8), chunk_min, chunk_max);
    std::exception_ptr err;
    std::mutex emu;

    auto worker = [&] {
        auto local = make_local();
        for (;;) {
            const size_t start = next.fetch_add(chunk);
            if (start >= n) break;
            const size_t end = std::min(n, start + chunk);
            for (size_t i = start; i < end; ++i) {
                try {
                    body(i, local);
                } catch (...) {
                    std::lock_guard<std::mutex> lk(emu);
                    if (!err) err = std::current_exception();
                    return;
                }
            }
        }
    };

    std::vector<std::thread> pool;
    pool.reserve(nt);
    for (unsigned t = 0; t < nt; ++t) pool.emplace_back(worker);
    for (auto& th : pool) th.join();
    if (err) std::rethrow_exception(err);
}

}  // namespace seqtree
