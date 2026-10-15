#include "seqtree/seqtree.hpp"
#include "seqtree/parallel.hpp"
#include "count_unit.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <unordered_map>

namespace seqtree {
namespace {
struct Cell { double score=-std::numeric_limits<double>::infinity(), raw=0; uint32_t gaps=0; };
Cell best(Cell a, const Cell& b, const Cell& c) {
    if (b.score>a.score) a=b;
    if (c.score>a.score) a=c;
    return a;
}
using Rows=std::array<std::vector<Cell>,6>;
// Existing strict NW ties: M/X/Y for M and X, M/Y/X for Y; terminal M/X/Y.
Cell nw(const std::vector<uint8_t>& q,const std::vector<uint8_t>& r,
        const std::vector<std::vector<double>>& sim,double gap,Rows& rows) {
    for (auto& row:rows) row.assign(r.size()+1,Cell{});
    auto& mp=rows[0]; auto& xp=rows[1]; auto& yp=rows[2];
    auto& mc=rows[3]; auto& xc=rows[4]; auto& yc=rows[5];
    mp[0]={0,0,0};
    for (size_t j=1;j<=r.size();++j) yp[j]={-gap*double(j),0,uint32_t(j)};
    for (size_t i=1;i<=q.size();++i) {
        mc[0]=yc[0]=Cell{}; xc[0]={-gap*double(i),0,uint32_t(i)};
        for (size_t j=1;j<=r.size();++j) {
            mc[j]=best(mp[j-1],xp[j-1],yp[j-1]);
            const double s=sim[q[i-1]][r[j-1]]; mc[j].score+=s; mc[j].raw+=s;
            xc[j]=best(mp[j],xp[j],yp[j]); xc[j].score-=gap; ++xc[j].gaps;
            yc[j]=best(mc[j-1],yc[j-1],xc[j-1]); yc[j].score-=gap; ++yc[j].gaps;
        }
        mp.swap(mc); xp.swap(xc); yp.swap(yc);
    }
    return best(mp.back(),xp.back(),yp.back());
}
void require(bool ok,const char* message) { if (!ok) throw std::invalid_argument(message); }
void csr(const std::vector<uint32_t>& offsets,size_t rows,size_t entries) {
    require(offsets.size()==rows+1 && offsets.front()==0 && offsets.back()==entries &&
        std::is_sorted(offsets.begin(),offsets.end()),"invalid CSR offsets");
}
void add(uint64_t& a,uint64_t b) {
    if (b>std::numeric_limits<uint64_t>::max()-a) throw std::overflow_error("annotation count overflow");
    a+=b;
}
}

AlignmentMaxResult Index::alignment_max_batch(const std::vector<std::string>& queries,
    const SearchParams& p,const std::vector<std::vector<double>>& similarity,const std::string& symbols,int max_indels,
    const std::vector<uint32_t>& qo,const std::vector<uint32_t>& qt,
    const std::vector<uint32_t>& ro,const std::vector<uint32_t>& rt,const std::vector<int32_t>& tc,
    const std::vector<uint32_t>& co,const std::vector<int32_t>& ci,
    const std::vector<std::vector<int32_t>>& qg,const std::vector<std::vector<int32_t>>& rg,
    const std::vector<std::vector<std::vector<double>>>& offsets,
    const std::vector<uint64_t>& other,const std::vector<uint64_t>& multiplicities,
    const std::vector<double>& coefficients,int threads,bool exclude_exact) const {
    const size_t Q=queries.size(),R=size(),L=rt.size(),T=tc.size(),M=qt.size();
    require(p.engine!=Engine::SeqTrie && p.mode==Mode::AllHits && !p.max_hits && !p.matrix && !p.pos_matrix &&
        p.gap_open==1 && p.gap_extend==1 && !p.max_score_penalty,"candidate scope requires uncapped unit-cost seqtm");
    require(!exclude_exact || p.max_insertions==p.max_deletions,"query puncture requires symmetric candidate indel caps");
    require(max_indels>=0 && max_indels<=65535,"max_indels must be in [0,65535]");
    require(coefficients.size()==3 && std::all_of(coefficients.begin(),coefficients.end(),[](double x){return std::isfinite(x);}),"three finite coefficients required");
    require(!symbols.empty() && similarity.size()==symbols.size(),"similarity dimension does not match symbols");
    std::array<int,256> code; code.fill(-1); double gap=0;
    for (size_t a=0;a<symbols.size();++a) {
        const auto c=codec().encode(symbols[a]);
        require(c!=Codec::kInvalid && code[c]<0,"invalid or duplicate similarity symbol"); code[c]=int(a);
        require(similarity[a].size()==symbols.size(),"similarity matrix must be square");
        for (double x:similarity[a]) { require(std::isfinite(x),"similarity values must be finite"); gap=std::max(gap,std::abs(x)); }
    }
    require(gap>0,"NW selection gap max(abs(similarity)) must be positive");
    csr(qo,Q,M); csr(ro,R,L); csr(co,R,ci.size());
    require(other.size()==L && multiplicities.size()==L,"annotation arrays have inconsistent lengths");
    require(qg.size()==offsets.size() && rg.size()==offsets.size(),"group offset components disagree");
    for (size_t c=0;c<offsets.size();++c) {
        require(qg[c].size()==Q && rg[c].size()==L,"group ID lengths disagree");
        const size_t columns=offsets[c].empty()?0:offsets[c][0].size();
        for (const auto& row:offsets[c]) require(row.size()==columns && std::all_of(row.begin(),row.end(),[](double x){return std::isfinite(x);}),"group offsets must be finite rectangular matrices");
        for (int32_t x:qg[c]) require(x>=-1 && (x<0 || size_t(x)<offsets[c].size()),"query group ID out of range");
        for (int32_t x:rg[c]) require(x>=-1 && (x<0 || size_t(x)<columns),"reference group ID out of range");
    }
    std::vector<size_t> seen(T,std::numeric_limits<size_t>::max());
    for (size_t i=0;i<Q;++i) for (size_t k=qo[i];k<qo[i+1];++k) {
        require(qt[k]<T && seen[qt[k]]!=i,"invalid or duplicate requested target"); seen[qt[k]]=i;
    }
    for (int32_t x:tc) require(x>=0,"target class IDs must be nonnegative");
    std::vector<uint8_t> modeled(L,1);
    for (size_t r=0;r<R;++r) {
        require(co[r]<co[r+1],"each reference requires at least one class role");
        std::vector<int32_t> labels(ci.begin()+co[r],ci.begin()+co[r+1]); std::sort(labels.begin(),labels.end());
        require(labels.front()>=0 && std::adjacent_find(labels.begin(),labels.end())==labels.end(),"reference class roles must be unique and nonnegative");
        for (size_t a=ro[r];a<ro[r+1];++a) {
            require(rt[a]<T && std::binary_search(labels.begin(),labels.end(),tc[rt[a]]),"annotation target class missing from reference roles");
            require(multiplicities[a]>0 && other[a]<ci.size() && other[a]<=R,"invalid annotation multiplicity or other-class count");
            for (size_t c=0;c<rg.size();++c) if (rg[c][a]<0) modeled[a]=0;
        }
    }
    auto encode=[&](std::string_view s) {
        require(s.size()<=size_t(std::numeric_limits<int>::max()),"sequence exceeds native search length");
        std::vector<uint8_t> out; out.reserve(s.size());
        for (char a:s) { const auto c=codec().encode(a); require(c!=Codec::kInvalid && code[c]>=0,"sequence symbol absent from similarity matrix"); out.push_back(uint8_t(code[c])); }
        return out;
    };
    auto self=[&](const std::vector<uint8_t>& s) {
        double value=0; for (auto a:s) value+=similarity[a][a];
        if (!std::isfinite(value)) throw std::overflow_error("self similarity overflow"); return value;
    };
    auto key=[](const std::vector<uint8_t>& s) { return std::string(s.begin(),s.end()); };
    std::vector<std::vector<uint8_t>> refs(R),qs(Q); std::vector<double> rself(R),qself(Q);
    std::unordered_map<std::string,uint32_t> identities;
    for (size_t r=0;r<R;++r) { refs[r]=encode(ref_seq(uint32_t(r))); rself[r]=self(refs[r]); require(identities.emplace(key(refs[r]),uint32_t(r)).second,"alignment max index references must be unique sequences"); }
    std::vector<int64_t> exact_ref(Q,-1);
    for (size_t i=0;i<Q;++i) { qs[i]=encode(queries[i]); qself[i]=self(qs[i]); if (auto found=identities.find(key(qs[i]));found!=identities.end()) exact_ref[i]=found->second; }
    AlignmentMaxResult out; const double nan=std::numeric_limits<double>::quiet_NaN();
    out.scores.assign(M,nan); out.raw_similarity.assign(M,nan); out.normalized_similarity.assign(M,nan);
    out.pfit.assign(M,nan); out.informativeness.assign(M,nan); out.gap_counts.assign(M,0);
    out.best_ref_ids.assign(M,-1); out.best_lane_ids.assign(M,-1);
    out.counts.assign(M,0); out.modeled_counts.assign(M,0); out.aa_counts.assign(M,0); out.exact_counts.assign(M,0);
    out.candidate_counts.assign(Q,0); out.nw_pairs.assign(Q,0); out.query_supported.assign(Q,1);
    struct Worker { Searcher search; std::vector<Hit> hits; Rows rows; std::vector<int64_t> slots; std::vector<uint64_t> aa_seen;
        Worker(const Index& index,size_t targets):search(index),slots(targets,-1),aa_seen(targets,0) {} };
    parallel_for(Q,threads,[&]{return Worker(*this,T);},[&](size_t i,Worker& w) {
        for (size_t c=0;c<qg.size();++c) if (qg[c][i]<0) out.query_supported[i]=0;
        for (size_t k=qo[i];k<qo[i+1];++k) w.slots[qt[k]]=int64_t(k);
        if (count_unit::bounded_eligible(p,max_indels) && !queries[i].empty() && queries[i].size()<=63) {
            count_unit::Pattern pattern(queries[i],codec()); w.hits.clear();
            for (uint32_t r=0;r<R;++r) {
                const int d=count_unit::bounded_distance(pattern,ref_seq(r),codec(),p,max_indels);
                // Only ref_id is consumed below; no operation trace is reconstructed.
                if (d<=p.max_total_edits) w.hits.push_back(Hit{r,d,0,0,0});
            }
        } else w.search.search_into(queries[i],p,w.hits,max_indels);
        out.candidate_counts[i]=w.hits.size();
        const int64_t er=exact_ref[i]; const size_t qb=er<0?0:co[size_t(er)],qe=er<0?0:co[size_t(er)+1];
        const uint64_t N=ci.size()-(exclude_exact?qe-qb:0);
        for (const auto& hit:w.hits) {
            const uint32_t r=hit.ref_id; const bool exact=int64_t(r)==er; bool align=false;
            for (size_t a=ro[r];a<ro[r+1];++a) {
                const uint32_t t=rt[a]; const int64_t k=w.slots[t]; if (k<0) continue;
                if (exact) add(out.exact_counts[size_t(k)],multiplicities[a]);
                if (exclude_exact && exact) continue;
                add(out.counts[size_t(k)],multiplicities[a]);
                if (w.aa_seen[t]!=uint64_t(r)+1) { ++out.aa_counts[size_t(k)]; w.aa_seen[t]=uint64_t(r)+1; }
                if (modeled[a] && out.query_supported[i]) { add(out.modeled_counts[size_t(k)],multiplicities[a]); align=true; }
            }
            if (!align) continue; // Filter requested targets and model availability BEFORE NW.
            const Cell features=nw(qs[i],refs[r],similarity,gap,w.rows); ++out.nw_pairs[i];
            const double normalized=features.raw-std::max(qself[i],rself[r]);
            if (!std::isfinite(features.score) || !std::isfinite(normalized)) throw std::overflow_error("NW similarity overflow");
            for (size_t a=ro[r];a<ro[r+1];++a) {
                const uint32_t t=rt[a]; const int64_t slot=w.slots[t]; if (slot<0 || !modeled[a]) continue;
                const size_t k=size_t(slot); bool query_other=false;
                for (size_t c=qb;c<qe;++c) if (ci[c]!=tc[t]) query_other=true;
                const uint64_t removed=exclude_exact && query_other?1:0;
                require(other[a]>=removed,"punctured other-class counter underflow");
                const uint64_t numerator=other[a]-removed+1;
                require(N>=numerator,"punctured reference-role denominator is too small");
                const double information=-std::log10(double(numerator)/double(N));
                double eta=coefficients[0];
                for (size_t c=0;c<offsets.size();++c) eta+=offsets[c][size_t(qg[c][i])][size_t(rg[c][a])];
                eta+=coefficients[1]*normalized+coefficients[2]*features.gaps;
                if (!std::isfinite(eta)) throw std::overflow_error("linear predictor overflow");
                const double probability=-std::expm1(-std::exp(eta)),score=probability*information;
                if (out.best_lane_ids[k]<0 || score>out.scores[k] || (score==out.scores[k] && int64_t(a)<out.best_lane_ids[k])) {
                    out.scores[k]=score; out.best_ref_ids[k]=r; out.best_lane_ids[k]=int64_t(a);
                    out.raw_similarity[k]=features.raw; out.normalized_similarity[k]=normalized; out.gap_counts[k]=features.gaps;
                    out.pfit[k]=probability; out.informativeness[k]=information;
                }
            }
        }
        for (size_t k=qo[i];k<qo[i+1];++k) {
            if (out.query_supported[i] && !out.counts[k]) out.scores[k]=0;
            w.slots[qt[k]]=-1; w.aa_seen[qt[k]]=0;
        }
    });
    return out;
}
} // namespace seqtree
