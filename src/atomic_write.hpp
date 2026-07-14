#pragma once
// Crash- and concurrency-safe file replacement.
//
// A serialized index is only useful if a reader never observes it half-written. Writing straight
// into the destination gives a reader a window -- the whole duration of the write -- in which the
// file exists but is truncated, which is exactly what a cold-cache fan-out hits: one worker is
// still saving while the next already sees the path and calls load(). Measured on a 45 MB control
// index that window is ~55 ms and it corrupted 10 out of 10 attempts.
//
// So write a uniquely-named temporary beside the destination and rename it into place. Rename is
// atomic on the same filesystem, on POSIX and on Windows alike, so a reader sees either the
// previous complete file or the new complete file -- never a partial one. On POSIX a reader that
// already holds the old file open keeps reading the old inode, which is fine.

#include <atomic>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <thread>

#ifdef _WIN32
#include <process.h>
#else
#include <unistd.h>
#endif

namespace seqtree::detail {

// Unique per process *and* per concurrent save within a process, so two savers never share a
// temporary. The pid alone is not enough: one process may save the same path from two threads.
inline std::string temp_path(const std::string& path) {
    static std::atomic<uint64_t> seq{0};
#ifdef _WIN32
    const long long pid = _getpid();
#else
    const long long pid = ::getpid();
#endif
    return path + ".tmp." + std::to_string(pid) + "." + std::to_string(seq.fetch_add(1));
}

// Serialize with `write` into a temporary, then atomically rename it over `path`. Throws on
// failure, and never leaves a stray temporary behind.
template <typename Writer>
void atomic_write(const std::string& path, Writer&& write) {
    const std::string tmp = temp_path(path);
    try {
        {
            std::ofstream os(tmp, std::ios::binary);
            if (!os) throw std::runtime_error("seqtree: cannot open '" + tmp + "' for writing");
            write(os);
            os.flush();
            if (!os) throw std::runtime_error("seqtree: write failed for '" + path + "'");
        }  // must close before renaming: Windows will not rename a file that is still open

        // Windows fails the rename if a reader currently holds the destination open, unlike
        // POSIX. That window is short (a reader opens, slurps, closes), so retry briefly rather
        // than failing a save that would have succeeded a millisecond later.
        std::error_code ec;
        for (int attempt = 0; attempt < 20; ++attempt) {
            std::filesystem::rename(tmp, path, ec);
            if (!ec) return;
            std::this_thread::sleep_for(std::chrono::milliseconds(5));
        }
        throw std::runtime_error("seqtree: cannot replace '" + path + "': " + ec.message());
    } catch (...) {
        std::error_code ec;
        std::filesystem::remove(tmp, ec);
        throw;
    }
}

}  // namespace seqtree::detail
