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
#include <filesystem>
#include <sstream>
#include <mutex>
#include <set>
#include <string>
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

inline unsigned cpu_rate_budget(unsigned available, unsigned system_cpus, unsigned rate) {
    if (!rate || rate > 10000) throw std::invalid_argument("invalid Windows CPU rate limit");
    return std::min(available, std::max(1u, unsigned(uint64_t(system_cpus) * rate / 10000)));
}

#ifdef _WIN32
inline std::pair<unsigned, unsigned> windows_cpu_budget() {
    const HANDLE process = GetCurrentProcess();
    USHORT group_count = 1, group = 0;
    if (!GetProcessGroupAffinity(process, &group_count, &group) || group_count != 1)
        throw std::runtime_error("CPU budgeting requires a single Windows processor group; bind the process to one group");
    DWORD_PTR mask = 0, system_mask = 0;
    if (!GetProcessAffinityMask(process, &mask, &system_mask) || !mask)
        throw std::runtime_error("cannot read Windows process CPU affinity");
    unsigned logical = 0;
    for (auto bits = mask; bits; bits &= bits - 1) ++logical;
    DWORD bytes = 0;
    GetLogicalProcessorInformationEx(RelationProcessorCore, nullptr, &bytes);
    if (!bytes) throw std::runtime_error("cannot read Windows physical CPU topology");
    std::vector<unsigned char> buffer(bytes);
    if (!GetLogicalProcessorInformationEx(RelationProcessorCore,
            reinterpret_cast<PSYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX>(buffer.data()), &bytes))
        throw std::runtime_error("cannot read Windows physical CPU topology");
    unsigned physical = 0;
    for (size_t offset = 0; offset < bytes;) {
        const auto* info = reinterpret_cast<const SYSTEM_LOGICAL_PROCESSOR_INFORMATION_EX*>(buffer.data() + offset);
        if (!info->Size || offset + info->Size > bytes)
            throw std::runtime_error("invalid Windows physical CPU topology");
        bool allowed = false;
        for (WORD i = 0; i < info->Processor.GroupCount; ++i)
            allowed |= info->Processor.GroupMask[i].Group == group &&
                       (info->Processor.GroupMask[i].Mask & mask) != 0;
        physical += allowed;
        offset += info->Size;
    }
    if (!physical) throw std::runtime_error("Windows affinity contains no physical cores");
    BOOL in_job = FALSE;
    if (!IsProcessInJob(process, nullptr, &in_job))
        throw std::runtime_error("cannot inspect Windows CPU job membership");
    if (in_job) {
        JOBOBJECT_CPU_RATE_CONTROL_INFORMATION rate{};
        if (!QueryInformationJobObject(nullptr, JobObjectCpuRateControlInformation, &rate, sizeof(rate), nullptr))
            throw std::runtime_error("cannot inspect Windows CPU job limit");
        if (rate.ControlFlags & JOB_OBJECT_CPU_RATE_CONTROL_ENABLE) {
            unsigned fraction = 0;
            if (rate.ControlFlags & JOB_OBJECT_CPU_RATE_CONTROL_HARD_CAP) fraction = rate.CpuRate;
            else if (rate.ControlFlags & JOB_OBJECT_CPU_RATE_CONTROL_MIN_MAX_RATE) fraction = rate.MaxRate;
            if (fraction) {
                const auto total = GetActiveProcessorCount(ALL_PROCESSOR_GROUPS);
                logical = cpu_rate_budget(logical, total, fraction);
            }
        }
    }
    return {logical, std::min(logical, physical)};
}
#endif

// Resolve once per batch; affinity/quota limits can change between calls.
inline unsigned available_cpus() {
    unsigned count = std::max(1u, std::thread::hardware_concurrency());
#ifdef __linux__
    cpu_set_t affinity;
    if (sched_getaffinity(0, sizeof(affinity), &affinity) == 0)
        count = std::min(count, unsigned(std::max(1, CPU_COUNT(&affinity))));
    namespace fs = std::filesystem;
    const fs::path root("/sys/fs/cgroup");
    std::vector<fs::path> bases{root, root / "cpu", root / "cpu,cpuacct", root / "cpuacct,cpu"};
    std::vector<fs::path> paths = bases;
    std::ifstream membership("/proc/self/cgroup");
    std::string line;
    while (std::getline(membership, line)) {
        const auto first = line.find(':');
        const auto second = line.find(':', first + 1);
        if (first == std::string::npos || second == std::string::npos) continue;
        const auto controllers = line.substr(first + 1, second - first - 1);
        if (!controllers.empty() && ("," + controllers + ",").find(",cpu,") == std::string::npos) continue;
        const fs::path relative = fs::path(line.substr(second + 1)).relative_path();
        bool safe = true;
        for (const auto& component : relative) if (component == "..") safe = false;
        if (!safe) continue;
        for (const auto& base : bases) {
            if (controllers.empty() && base != root) continue;
            for (auto path = base / relative; path != base; path = path.parent_path())
                paths.push_back(path);
        }
    }
    for (const auto& path : paths) {
        std::string quota_text;
        uint64_t period = 0;
        std::ifstream v2(path / "cpu.max");
        if (!(v2 >> quota_text >> period)) {
            std::ifstream v1quota(path / "cpu.cfs_quota_us"), v1period(path / "cpu.cfs_period_us");
            if (!(v1quota >> quota_text) || !(v1period >> period)) continue;
        }
        if (quota_text == "max" || quota_text == "-1" || !period) continue;
        uint64_t quota = 0;
        std::istringstream parsed(quota_text);
        if (parsed >> quota && quota)
            count = unsigned(std::min<uint64_t>(count, std::max<uint64_t>(1, quota / period)));
    }
    for (const char* name : {"SLURM_CPUS_PER_TASK", "SLURM_CPUS_ON_NODE"}) {
        if (const char* value = std::getenv(name)) {
            char* end = nullptr;
            const auto allocated = std::strtoul(value, &end, 10);
            if (end != value && !*end && allocated)
                count = unsigned(std::min<unsigned long>(count, allocated));
        }
    }
#elif defined(_WIN32)
    count = std::min(count, windows_cpu_budget().first);
#endif
    return count;
}

inline unsigned physical_cores_for_cpus(const std::vector<unsigned>& cpus,
                                       unsigned available,
                                       const std::string& root = "/sys/devices/system/cpu/") {
    std::set<std::pair<int, int>> cores;
    for (unsigned cpu : cpus) {
        const auto path = root + "cpu" + std::to_string(cpu) + "/topology/";
        int package, core;
        if (!(std::ifstream(path + "physical_package_id") >> package) ||
            !(std::ifstream(path + "core_id") >> core) || package < 0 || core < 0) return available;
        cores.emplace(package, core);
    }
    return cores.empty() ? available : std::min(available, unsigned(cores.size()));
}

// Topology is separate from the logical allocation: explicit budgets may use SMT.
inline unsigned physical_cpus() {
    const unsigned available = available_cpus();
#ifdef __linux__
    cpu_set_t affinity;
    if (sched_getaffinity(0, sizeof(affinity), &affinity) == 0) {
        std::vector<unsigned> cpus;
        for (int cpu = 0; cpu < CPU_SETSIZE; ++cpu)
            if (CPU_ISSET(cpu, &affinity)) cpus.push_back(unsigned(cpu));
        return physical_cores_for_cpus(cpus, available);
    }
#elif defined(__APPLE__)
    unsigned physical = available;
    size_t size = sizeof(physical);
    if (sysctlbyname("hw.physicalcpu", &physical, &size, nullptr, 0) == 0 && physical)
        return std::min(available, physical);
#elif defined(_WIN32)
    return windows_cpu_budget().second;
#endif
    return available;
}

// Runs body(i, local) for every i in [0, n), across nt workers pulling contiguous chunks off
// one atomic cursor. `threads <= 0` selects available physical cores; positive counts may use logical CPUs.
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
    unsigned nt = threads > 0 ? unsigned(threads) : physical_cpus();
    nt = std::min({nt, available, 999u});
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
