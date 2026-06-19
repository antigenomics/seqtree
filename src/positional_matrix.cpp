#include "seqtree/seqtree.hpp"

#include <stdexcept>

namespace seqtree {

PositionalMatrix PositionalMatrix::from_weights(const SubstitutionMatrix& base,
                                                const std::vector<int32_t>& weights) {
    PositionalMatrix m;
    m.size_ = base.size();
    m.width_ = static_cast<uint16_t>(weights.size());
    m.pen_.resize(size_t(m.width_) * m.size_ * m.size_);
    m.masked_.assign(m.width_, 0);
    for (uint16_t p = 0; p < m.width_; ++p) {
        int32_t w = weights[p];
        m.masked_[p] = (w == 0) ? 1 : 0;
        for (uint8_t a = 0; a < m.size_; ++a)
            for (uint8_t b = 0; b < m.size_; ++b)
                m.pen_[(size_t(p) * m.size_ + a) * m.size_ + b] = w * base.penalty(a, b);
    }
    return m;
}

PositionalMatrix PositionalMatrix::from_tables(uint8_t size, uint16_t width,
                                               const std::vector<int32_t>& data,
                                               const std::vector<uint8_t>& masked) {
    if (data.size() != size_t(width) * size * size)
        throw std::invalid_argument("PositionalMatrix::from_tables: data size mismatch");
    PositionalMatrix m;
    m.size_ = size;
    m.width_ = width;
    m.pen_ = data;
    if (masked.empty())
        m.masked_.assign(width, 0);
    else if (masked.size() == width)
        m.masked_ = masked;
    else
        throw std::invalid_argument("PositionalMatrix::from_tables: masked size mismatch");
    return m;
}

}  // namespace seqtree
