// Private prototype: exact unit histogram acceleration; not installed or integrated.
#pragma once
#include "seqtree/seqtree.hpp"
#include <array>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cstring>

namespace seqtree::count_unit {
struct Pattern {
    std::array<uint64_t,256> eq{};
    std::array<uint8_t,63> codes{};
    size_t length;
    Pattern(std::string_view query, const Codec& codec):length(query.size()) {
        if (!length || length>63) throw std::invalid_argument("prototype pattern length must be 1..63");
        for (size_t i=0;i<length;++i) {
            codes[i]=codec.encode(query[i]);
            if (codes[i]==Codec::kInvalid) throw std::invalid_argument("invalid query symbol");
        }
        // Include all byte aliases, preserving Codec case folding and IUPAC symbols.
        for (unsigned c=0;c<256;++c)
            for (size_t i=0;i<length;++i)
                if (codec.encode(char(c))==codes[i]) eq[c]|=uint64_t(1)<<i;
    }
    int distance(std::string_view ref) const {
        uint64_t pv=~uint64_t(0),mv=0,high=uint64_t(1)<<(length-1);
        int score=int(length);
        for (unsigned char c:ref) {
            const uint64_t e=eq[c],xv=e|mv,xh=(((e&pv)+pv)^pv)|e;
            uint64_t ph=mv|~(xh|pv),mh=pv&xh;
            score+=bool(ph&high)-bool(mh&high);
            ph=(ph<<1)|1;mh<<=1;
            pv=mh|~(xv|ph);mv=ph&xv;
        }
        return score;
    }
};
// Minimal substitution count at each prefix and exact insertion/deletion counts.
// Fixed scratch: no allocation; j=i+ins-del determines the reference prefix.
inline int restricted(const Pattern& q,std::string_view r,const Codec& codec,int cap) {
    uint8_t dp[64][3][3];
    std::memset(dp,255,sizeof(dp));
    dp[0][0][0]=0;
    int best=cap+1;
    for (int i=0;i<=int(q.length);++i) for (int ins=0;ins<=2;++ins) for (int del=0;del<=2;++del) {
        int j=i+ins-del,s=dp[i][ins][del];
        if (j<0 || j>int(r.size()) || s+ins+del>cap) continue;
        if (i==int(q.length) && j==int(r.size())) best=std::min(best,s+ins+del);
        auto put=[&](int ii,int ni,int nd,int ns) {
            if (ns+ni+nd<=cap) dp[ii][ni][nd]=std::min(dp[ii][ni][nd],uint8_t(ns));
        };
        if (i<int(q.length) && j<int(r.size())) put(i+1,ins,del,s+(q.codes[i]!=codec.encode(r[j])));
        if (ins<2 && j<int(r.size())) put(i,ins+1,del,s);
        if (del<2 && i<int(q.length)) put(i+1,ins,del+1,s);
    }
    return best;
}
inline int exact_distance(const Pattern& q,std::string_view r,const Codec& codec,int cap) {
    if (r.size()>q.length+2 || q.length>r.size()+2) return cap+1;
    const int delta=int(r.size())-int(q.length);
    const int d=q.distance(r);
    if (d>cap) return cap+1;
    // For any path, ins+del<=d and |ins-del|=|delta|. Thus max(ins,del)
    // <=floor((d+|delta|)/2); both per-type caps are automatic when <=2.
    if (d+std::abs(delta)<=5) return d;
    return restricted(q,r,codec,cap);
}
// Bounded installed-API profiles favor trie traversal at total caps1/2.
inline bool eligible(const SearchParams& p) {
    return !p.matrix && !p.pos_matrix && p.gap_open==1 &&
        p.max_total_edits>=3 && p.max_total_edits<=5 &&
        p.max_substitutions>=p.max_total_edits && p.max_insertions==2 && p.max_deletions==2 &&
        (p.max_score_penalty<=0 || p.max_score_penalty>=p.max_total_edits) &&
        p.engine!=Engine::SeqTrie && p.mode==Mode::AllHits && !p.max_hits;
}
inline bool fill_histogram(const Index& idx,std::string_view query,const SearchParams& p,
                           bool exclude_exact,std::vector<uint64_t>& out) {
    if (!eligible(p) || query.empty() || query.size()>63) return false;
    Pattern pattern(query,idx.codec());
    for(uint32_t id=0;id<idx.size();++id) {
        const int d=exact_distance(pattern,idx.ref_seq(id),idx.codec(),p.max_total_edits);
        if (d<=p.max_total_edits && (!exclude_exact || d)) ++out[d];
    }
    return true;
}
} // namespace seqtree::count_unit
