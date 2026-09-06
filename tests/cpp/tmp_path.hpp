#pragma once
// Tests write real files. "/tmp/..." is not a path on Windows, which the windows-latest CI job
// found the moment it was added, so ask the standard library where scratch files go.

#include <filesystem>
#include <string>

inline std::string tmp_path(const std::string& name) {
    return (std::filesystem::temp_directory_path() / name).string();
}
