#include "seqtree/seqtree.hpp"

namespace seqtree {

#include "blosum62.inc"
#include "pam250.inc"
#include "pam100.inc"
#include "structural.inc"

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
    // Gram -> squared-distance transform: treating sim as an inner product
    // s(a,b)=<phi(a),phi(b)>, this is ||phi(a)-phi(b)||^2 = s_aa + s_bb - 2 s_ab.
    // Symmetric, 0 on the identity, and >= 0 for BLOSUM/PAM (the diagonal is each
    // row's max, so s_ab <= min(s_aa,s_bb)). The clamp guards pathological pairs.
    for (uint8_t a = 0; a < size; ++a) {
        for (uint8_t b = 0; b < size; ++b) {
            int32_t p = sim[a * size + a] + sim[b * size + b] - 2 * sim[a * size + b];
            m.pen_[a * size + b] = p < 0 ? 0 : p;
        }
    }
    return m;
}

SubstitutionMatrix SubstitutionMatrix::blosum62() {
    return from_similarity(static_cast<uint8_t>(kBlosum62Size), kBlosum62);
}

SubstitutionMatrix SubstitutionMatrix::pam250() {
    return from_similarity(static_cast<uint8_t>(kPam250Size), kPam250);
}

SubstitutionMatrix SubstitutionMatrix::pam100() {
    return from_similarity(static_cast<uint8_t>(kPam100Size), kPam100);
}

SubstitutionMatrix SubstitutionMatrix::structural() {
    return from_similarity(static_cast<uint8_t>(kStructuralSize), kStructural);
}

}  // namespace seqtree
