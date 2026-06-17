#include "seqtree/seqtree.hpp"

#include <algorithm>

namespace seqtree {

#include "blosum62.inc"
#include "pam50.inc"

SubstitutionMatrix SubstitutionMatrix::unit(uint8_t size) {
    SubstitutionMatrix m;
    m.size_ = size;
    m.pen_.assign(size_t(size) * size, 1);
    for (uint8_t i = 0; i < size; ++i) m.pen_[i * size + i] = 0;
    return m;
}

SubstitutionMatrix SubstitutionMatrix::from_similarity(uint8_t size, const int32_t* sim) {
    SubstitutionMatrix m;
    m.size_ = size;
    m.pen_.resize(size_t(size) * size);
    for (uint8_t a = 0; a < size; ++a) {
        for (uint8_t b = 0; b < size; ++b) {
            int32_t best = std::max(sim[a * size + a], sim[b * size + b]);
            int32_t p = best - sim[a * size + b];  // >= 0; identity => 0
            m.pen_[a * size + b] = p < 0 ? 0 : p;
        }
    }
    return m;
}

SubstitutionMatrix SubstitutionMatrix::blosum62() {
    return from_similarity(static_cast<uint8_t>(kBlosum62Size), kBlosum62);
}

SubstitutionMatrix SubstitutionMatrix::pam50() {
    return from_similarity(static_cast<uint8_t>(kPam50Size), kPam50);
}

}  // namespace seqtree
