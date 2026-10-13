#pragma once
#include "seqtree/seqtree.hpp"
#include <array>
namespace seqtree::detail {
using PairwiseScratch = std::array<std::vector<int64_t>, 6>;
int64_t align_score_encoded(const uint8_t*, size_t, const uint8_t*, size_t,
    const SubstitutionMatrix&, AlignMode, int32_t, int32_t, PairwiseScratch&);
}
