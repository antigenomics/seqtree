#pragma once
// One parallel-for scaffold for the whole tree. Every batch entry point in seqtree used to
// carry its own copy of this: an atomic cursor, an adaptive chunk, an exception_ptr, a vector
// of threads, and a rethrow after the join. Seven copies drifted -- one of them had no
// exception plumbing at all.
#include <algorithm>
#include <atomic>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

#ifdef __linux__
#include <sched.h>
#elif defined(__APPLE__)
#include <sys/sysctl.h>
#elif defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#include <windows.h>
#endif

namespace seqtree {

// Resolve once per batch; affinity/quota limits can change between calls.
inline unsigned available_cpus() {
    unsigned count = std::max(1u, std::thread::hardware_concurrency());
#ifdef __linux__
    cpu_set_t affinity;
    if (sched_getaffinity(0, sizeof(affinity), &affinity) == 0)
        count = std::min(count, unsigned(std::max(1, CPU_COUNT(&affinity))));
    std::ifstream limits("/sys/fs/cgroup/cpu.max");
    uint64_t quota, period;
    if (limits >> quota >> period && period)
        count = unsigned(std::min<uint64_t>(count, std::max<uint64_t>(1, quota / period)));
    for (const char* name : {"SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"}) {
        if (const char* value = std::getenv(name)) {
            char* end = nullptr;
            const auto allocated = std::strtoul(value, &end, 10);
            if (end != value && !*end && allocated)
                count = unsigned(std::min<unsigned long>(count, allocated));
        }
    }
#elif defined(__APPLE__)
    unsigned physical = count;
    size_t size = sizeof(physical);
    if (sysctlbyname("hw.physicalcpu", &physical, &size, nullptr, 0) == 0 && physical)
        count = std::min(count, physical);
#elif defined(_WIN32)
    DWORD_PTR process_mask, system_mask;
    if (GetProcessAffinityMask(GetCurrentProcess(), &process_mask, &system_mask)) {
        unsigned allowed = 0;
        for (; process_mask; process_mask &= process_mask - 1) ++allowed;
        if (allowed) count = std::min(count, allowed);
    }
#endif
    return count;
}

// Runs body(i, local) for every i in [0, n), across nt workers pulling contiguous chunks off
// one atomic cursor. `threads <= 0` means the detected available CPU budget.
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
    const unsigned available = available_cpus();
    unsigned nt = threads > 0 ? unsigned(threads) : available;
    if (nt > available || nt >= 1000)
        throw std::invalid_argument("worker count exceeds the available CPU budget; request fewer threads");
    if (n == 0) return;
    nt = unsigned(std::min<size_t>(nt, n));  // Do not start idle workers for absent work.

    // A serial budget runs on the calling thread; spawning one worker only adds overhead.
    if (nt == 1) {
        auto local = make_local();
        for (size_t i = 0; i < n; ++i) body(i, local);
        return;
    }

    std::atomic<size_t> next{0};
    const size_t chunk = std::clamp<size_t>(n / (size_t(nt) * 8), chunk_min, chunk_max);
    std::exception_ptr err;
    std::mutex emu;

    auto worker = [&] {
        // Scratch construction can allocate and throw too. Nothing may escape the entry point.
        try {
            auto local = make_local();
            for (;;) {
                const size_t start = next.fetch_add(chunk);
                if (start >= n) break;
                const size_t end = std::min(n, start + chunk);
                for (size_t i = start; i < end; ++i) {
                    body(i, local);
                }
            }
        } catch (...) {
            std::lock_guard<std::mutex> lk(emu);
            if (!err) err = std::current_exception();
        }
    };

    // A failed thread launch must also join workers already started before unwinding.
    std::vector<std::jthread> pool;
    pool.reserve(nt);
    // The caller is one compute thread: a budget B creates only B-1 OS workers.
    for (unsigned t = 1; t < nt; ++t) pool.emplace_back(worker);
    worker();
    for (auto& th : pool) th.join();
    if (err) std::rethrow_exception(err);
}

}  // namespace seqtree
