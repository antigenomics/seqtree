#include "seqtree/seqtree.hpp"

#include <cstring>

namespace seqtree {

namespace {
// AminoAcid order matches the BLOSUM62 table in blosum62.inc.
constexpr const char* kAminoAcid = "ARNDCQEGHILKMFPSTWYVBZX*";
constexpr const char* kNucleotide = "ACGT";
constexpr const char* kNucleotideIUPAC = "ACGTRYSWKMBDHVN";

const char* symbols_for(Alphabet a) {
    switch (a) {
        case Alphabet::AminoAcid:       return kAminoAcid;
        case Alphabet::Nucleotide:      return kNucleotide;
        case Alphabet::NucleotideIUPAC: return kNucleotideIUPAC;
    }
    return kAminoAcid;
}
}  // namespace

Codec::Codec(Alphabet a) : alphabet_(a) {
    std::memset(enc_, kInvalid, sizeof(enc_));
    std::memset(dec_, 0, sizeof(dec_));
    const char* syms = symbols_for(a);
    size_ = 0;
    for (const char* p = syms; *p; ++p, ++size_) {
        char up = *p;
        char lo = (up >= 'A' && up <= 'Z') ? char(up - 'A' + 'a') : up;
        enc_[static_cast<uint8_t>(up)] = size_;
        enc_[static_cast<uint8_t>(lo)] = size_;  // case-insensitive
        dec_[size_] = up;
    }
}

}  // namespace seqtree
