#include "seqtree/seqtree.hpp"

#include <algorithm>
#include <vector>

namespace seqtree {

#include "blosum62.inc"
#include "blosum45.inc"
#include "blosum80.inc"
#include "pam250.inc"
#include "pam100.inc"
#include "structural.inc"

SubstitutionMatrix SubstitutionMatrix::unit(uint8_t size) {
    SubstitutionMatrix m;
    m.size_ = size;
    m.pen_.assign(size_t(size) * size, 1);
    for (uint8_t i = 0; i < size; ++i) m.pen_[i * size + i] = 0;
    // The unit matrix has no underlying log-odds, so its similarity is the plain +1/0
    // match score. Note this is the one matrix where sim and pen are NOT related by the
    // Gram transform (that would give an off-diagonal penalty of 2, not 1); they are simply
    // the two conventions for "identical or not".
    m.sim_.assign(size_t(size) * size, 0);
    for (uint8_t i = 0; i < size; ++i) m.sim_[i * size + i] = 1;
    return m;
}

SubstitutionMatrix SubstitutionMatrix::from_similarity(uint8_t size, const int32_t* sim) {
    SubstitutionMatrix m;
    m.size_ = size;
    m.pen_.resize(size_t(size) * size);
    // The raw log-odds are kept as well as the penalty. The Gram transform below is lossy --
    // it forces the diagonal to zero, so s(a,a) is destroyed and cannot be recovered from
    // pen alone. Search and the E-value ball need the penalty (>= 0, zero on the identity);
    // a Needleman-Wunsch / Smith-Waterman aligner needs the raw similarity, which is signed
    // and maximised. Both are views of the same matrix.
    m.sim_.assign(sim, sim + size_t(size) * size);
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

SubstitutionMatrix SubstitutionMatrix::blosum45() {
    return from_similarity(static_cast<uint8_t>(kBlosum45Size), kBlosum45);
}

SubstitutionMatrix SubstitutionMatrix::blosum80() {
    return from_similarity(static_cast<uint8_t>(kBlosum80Size), kBlosum80);
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

int32_t SubstitutionMatrix::scale() const {
    if (size_ < 2) return 0;
    std::vector<int32_t> off;
    off.reserve(size_t(size_) * (size_ - 1));
    for (uint8_t a = 0; a < size_; ++a)
        for (uint8_t b = 0; b < size_; ++b)
            if (a != b) off.push_back(penalty(a, b));
    auto mid = off.begin() + off.size() / 2;
    std::nth_element(off.begin(), mid, off.end());
    return *mid;
}

}  // namespace seqtree
