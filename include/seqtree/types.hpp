#pragma once
#include <cstdint>
#include <string>

namespace seqtree {

// Which symbol set the sequences use. AminoAcid order matches BLOSUM62 so a
// substitution matrix can be indexed directly by symbol code.
enum class Alphabet : uint8_t { AminoAcid, Nucleotide, NucleotideIUPAC };

// Which search driver to use. Auto picks seqtm for substitution-only / small-k
// indel work and seqtrie for matrix-weighted budgets (see searcher.cpp).
enum class Engine : uint8_t { Auto, SeqTrie, SeqTm };

// AllHits/TopHit are global (whole-query) matches; Local returns, per reference,
// the best-scoring fixed-width window of the query (trim/shift, e.g. class-II core).
enum class Mode : uint8_t { AllHits, TopHit, Local };

class SubstitutionMatrix;  // defined in seqtree.hpp
class PositionalMatrix;    // defined in seqtree.hpp; per-position penalties pen(pos,a,b)

// All limits are caps; 0 on a per-type cap means "zero of that type allowed".
// max_total_edits == 0 means "no separate total cap" (derive from per-type sum).
// max_score_penalty <= 0 means "no explicit budget" (derive from edit caps in
// unit-cost mode, unbounded in matrix mode).
struct SearchParams {
    Engine engine = Engine::Auto;
    Mode   mode   = Mode::AllHits;

    uint16_t max_substitutions = 0;
    uint16_t max_insertions    = 0;
    uint16_t max_deletions     = 0;
    uint16_t max_total_edits   = 0;

    int32_t  max_score_penalty = 0;            // explicit budget; <=0 => derived
    const SubstitutionMatrix* matrix = nullptr; // null => unit cost
    // Optional per-position penalties; used only on the substitution/Hamming path
    // (no indels), where query position is unambiguous. Masked positions (weight 0)
    // are free and do not count as substitutions. Overrides `matrix` when set.
    const PositionalMatrix* pos_matrix = nullptr;
    int32_t  gap_open   = 1;                    // linear gap cost per indel (v1)
    int32_t  gap_extend = 1;                    // reserved for affine gaps (roadmap)

    uint32_t max_hits = 0;                      // 0 => unlimited
};

// Payload-agnostic. Consumers map ref_id back to their own payload + filter.
// score is a non-negative penalty (0 == exact match). n_* are exact for seqtm,
// 0 / derived-on-demand for seqtrie.
struct Hit {
    uint32_t ref_id = 0;
    int32_t  score  = 0;
    uint16_t n_subs = 0;
    uint16_t n_ins  = 0;
    uint16_t n_dels = 0;
};

struct Alignment {
    std::string aligned_query;
    std::string aligned_ref;
    std::string ops;   // per column: 'M' match, 'S' substitution, 'I' insertion, 'D' deletion
    int32_t     score = 0;
};

}  // namespace seqtree
