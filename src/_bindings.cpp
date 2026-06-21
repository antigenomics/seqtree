#include "seqtree/seqtree.hpp"
#include "seqtree/kmer_index.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

namespace py = pybind11;
using namespace seqtree;

namespace {

// Python-facing params: strings at the edge, translated to the C++ struct per call.
struct PyParams {
    int max_subs = 0, max_ins = 0, max_dels = 0, max_total_edits = 0;
    long max_penalty = 0;
    int gap_open = 1, gap_extend = 1;
    std::string matrix;          // named builtin: "" (unit/identity), blosum62/pam250/pam100/structural
    std::optional<SubstitutionMatrix> matrix_obj;  // explicit/custom matrix (overrides name)
    std::optional<PositionalMatrix> pos_matrix_obj;  // per-position penalties (Hamming path)
    std::string engine = "auto"; // auto | seqtrie | seqtm
    std::string mode = "all";    // all | top | local
};

std::string lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c) { return std::tolower(c); });
    return s;
}

Alphabet parse_alphabet(const std::string& a) {
    std::string l = lower(a);
    if (l == "aa" || l == "amino" || l == "protein") return Alphabet::AminoAcid;
    if (l == "nt" || l == "dna" || l == "nucleotide") return Alphabet::Nucleotide;
    if (l == "nt_iupac" || l == "iupac") return Alphabet::NucleotideIUPAC;
    throw py::value_error("unknown alphabet '" + a + "' (use 'aa', 'nt', or 'iupac')");
}

Engine parse_engine(const std::string& e) {
    std::string l = lower(e);
    if (l == "auto") return Engine::Auto;
    if (l == "seqtrie") return Engine::SeqTrie;
    if (l == "seqtm") return Engine::SeqTm;
    throw py::value_error("unknown engine '" + e + "' (use 'auto', 'seqtrie', or 'seqtm')");
}

Mode parse_mode(const std::string& m) {
    std::string l = lower(m);
    if (l == "all") return Mode::AllHits;
    if (l == "top") return Mode::TopHit;
    if (l == "local") return Mode::Local;
    throw py::value_error("unknown mode '" + m + "' (use 'all', 'top', or 'local')");
}

// Symbols in codec code order for an alphabet (custom matrices must match this order).
std::string alphabet_symbols(Alphabet a) {
    Codec c(a);
    std::string s;
    for (uint8_t i = 0; i < c.size(); ++i) s += c.decode(i);
    return s;
}

// Built-in matrix names. "identity" is the unit matrix (any alphabet); the rest are
// amino-acid only. Keep this list in sync with the SubstitutionMatrix factories.
constexpr const char* kMatrixNames = "'identity', 'BLOSUM62', 'PAM250', 'PAM100', 'structural'";

bool is_matrix_name(const std::string& l) {
    return l.empty() || l == "identity" || l == "blosum62" || l == "pam250"
        || l == "pam100" || l == "structural";
}

SubstitutionMatrix named_matrix(const std::string& l, Alphabet a) {
    if (l == "identity") return SubstitutionMatrix::unit(Codec(a).size());
    if (a != Alphabet::AminoAcid)
        throw py::value_error(l + " requires the amino-acid alphabet");
    if (l == "blosum62") return SubstitutionMatrix::blosum62();
    if (l == "pam250") return SubstitutionMatrix::pam250();
    if (l == "pam100") return SubstitutionMatrix::pam100();
    return SubstitutionMatrix::structural();  // l == "structural"
}

// Returns nullopt for unit cost; throws for an unknown name or alphabet mismatch.
// An explicit matrix object (custom or built via SubstitutionMatrix factories) wins
// over the named builtin; we only check that its size matches the alphabet.
std::optional<SubstitutionMatrix> make_matrix(const PyParams& pp, Alphabet a) {
    if (pp.matrix_obj) {
        if (pp.matrix_obj->size() != Codec(a).size())
            throw py::value_error("matrix size does not match the alphabet");
        return *pp.matrix_obj;
    }
    if (pp.matrix.empty()) return std::nullopt;
    std::string l = lower(pp.matrix);
    if (!is_matrix_name(l))
        throw py::value_error("unknown matrix '" + pp.matrix + "' (use '', " +
                              kMatrixNames + ", or a SubstitutionMatrix)");
    return named_matrix(l, a);
}

// Accept either a builtin name (see kMatrixNames) or a SubstitutionMatrix.
void set_matrix(PyParams& p, const py::object& m) {
    p.matrix.clear();
    p.matrix_obj.reset();
    if (m.is_none()) return;
    if (py::isinstance<py::str>(m)) {
        std::string name = m.cast<std::string>();
        if (!is_matrix_name(lower(name)))
            throw py::value_error("unknown matrix '" + name + "' (use '', " +
                                  kMatrixNames + ", or a SubstitutionMatrix)");
        p.matrix = std::move(name);
    } else if (py::isinstance<SubstitutionMatrix>(m)) {
        p.matrix_obj = m.cast<SubstitutionMatrix>();
    } else {
        throw py::type_error("matrix must be a name string or a SubstitutionMatrix");
    }
}

SearchParams to_cpp(const PyParams& pp, const SubstitutionMatrix* mat) {
    SearchParams p;
    p.engine = parse_engine(pp.engine);
    p.mode = parse_mode(pp.mode);
    p.max_substitutions = uint16_t(pp.max_subs);
    p.max_insertions = uint16_t(pp.max_ins);
    p.max_deletions = uint16_t(pp.max_dels);
    p.max_total_edits = uint16_t(pp.max_total_edits);
    p.max_score_penalty = int32_t(pp.max_penalty);
    p.gap_open = pp.gap_open;
    p.gap_extend = pp.gap_extend;
    p.matrix = mat;
    p.pos_matrix = pp.pos_matrix_obj ? &*pp.pos_matrix_obj : nullptr;
    return p;
}

py::list hits_to_list(const std::vector<Hit>& hits) {
    py::list out(hits.size());
    for (size_t i = 0; i < hits.size(); ++i) out[i] = py::cast(hits[i]);
    return out;
}

py::list py_search(Index& idx, const std::string& q, const PyParams& pp) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    Searcher s(idx);
    return hits_to_list(s.search(q, cp));
}

py::list py_search_top(Index& idx, const std::string& q, const PyParams& pp, int k) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    cp.mode = Mode::TopHit;
    cp.max_hits = uint32_t(k < 1 ? 1 : k);
    Searcher s(idx);
    return hits_to_list(s.search(q, cp));
}

py::list py_search_batch(const Index& idx, const std::vector<std::string>& queries,
                         const PyParams& pp, int threads) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    std::vector<std::vector<Hit>> results;
    {
        py::gil_scoped_release release;  // pure C++ region, no Python objects touched
        results = idx.search_batch(queries, cp, threads);
    }
    py::list out(results.size());
    for (size_t i = 0; i < results.size(); ++i) out[i] = hits_to_list(results[i]);
    return out;
}

std::vector<uint64_t> py_collisions_batch(const Index& idx, const std::vector<std::string>& queries,
                                          const PyParams& pp, int threads) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    py::gil_scoped_release release;
    return idx.collisions_batch(queries, cp, threads);
}

py::list py_pairwise_batch(const std::vector<std::string>& a, const std::vector<std::string>& b,
                           const PyParams& pp, const std::string& alphabet, int threads) {
    Alphabet alph = parse_alphabet(alphabet);
    auto mat = make_matrix(pp, alph);
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    std::vector<std::vector<Hit>> results;
    {
        py::gil_scoped_release release;
        results = pairwise_batch(a, b, alph, cp, threads);
    }
    py::list out(results.size());
    for (size_t i = 0; i < results.size(); ++i) out[i] = hits_to_list(results[i]);
    return out;
}

Alignment py_align(const Index& idx, uint32_t ref_id, const std::string& q, const PyParams& pp) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    return idx.align(q, ref_id, cp);
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "seqtree: fuzzy biological-sequence search (C++ core)";

    py::class_<SubstitutionMatrix>(m, "SubstitutionMatrix",
                                   "Non-negative substitution penalties (penalty(a,a)==0). Build a "
                                   "named builtin (``blosum62``/``pam250``/``pam100``/``structural``, "
                                   "or ``unit`` for identity) or a custom one from a similarity grid "
                                   "whose row/column order matches ``amino_acids()`` (or "
                                   "``alphabet_symbols(alphabet)``).")
        .def_static("blosum62", &SubstitutionMatrix::blosum62)
        .def_static("pam250", &SubstitutionMatrix::pam250)
        .def_static("pam100", &SubstitutionMatrix::pam100)
        .def_static("structural", &SubstitutionMatrix::structural)
        .def_static("unit", &SubstitutionMatrix::unit, py::arg("size"))
        .def_static(
            "from_similarity",
            [](const std::vector<std::vector<int32_t>>& grid) {
                size_t n = grid.size();
                if (n == 0 || n > 32) throw py::value_error("matrix size must be 1..32");
                std::vector<int32_t> flat;
                flat.reserve(n * n);
                for (const auto& row : grid) {
                    if (row.size() != n) throw py::value_error("similarity matrix must be square");
                    flat.insert(flat.end(), row.begin(), row.end());
                }
                return SubstitutionMatrix::from_similarity(uint8_t(n), flat.data());
            },
            py::arg("grid"),
            "Build from a square similarity grid (higher == more similar), converted to "
            "penalties via max(sim[a,a], sim[b,b]) - sim[a,b]. Row/column order must match the "
            "target alphabet's symbol order (see ``amino_acids()``).")
        .def("size", &SubstitutionMatrix::size)
        .def("__repr__", [](const SubstitutionMatrix& s) {
            return "SubstitutionMatrix(size=" + std::to_string(s.size()) + ")";
        });

    py::class_<PositionalMatrix>(m, "PositionalMatrix",
                                 "Per-position penalties pen(pos, a, b) over a fixed frame width. "
                                 "Build from a base SubstitutionMatrix and per-position integer "
                                 "weights: weight 0 masks the position (free, not counted as a "
                                 "substitution -- e.g. an anchor); >1 up-weights it (e.g. a TCR "
                                 "hotspot). Used on the seqtm Hamming path when width == query "
                                 "length.")
        .def_static("from_weights", &PositionalMatrix::from_weights,
                    py::arg("base"), py::arg("weights"),
                    "pen[pos][a][b] = weights[pos] * base.penalty(a, b); weight 0 masks the "
                    "position. len(weights) is the frame width.")
        .def("size", &PositionalMatrix::size)
        .def("width", &PositionalMatrix::width)
        .def("masked", &PositionalMatrix::masked, py::arg("pos"))
        .def("penalty", &PositionalMatrix::penalty, py::arg("pos"), py::arg("a"), py::arg("b"))
        .def("__repr__", [](const PositionalMatrix& p) {
            return "PositionalMatrix(size=" + std::to_string(p.size()) +
                   ", width=" + std::to_string(p.width()) + ")";
        });

    m.def("alphabet_symbols", [](const std::string& a) { return alphabet_symbols(parse_alphabet(a)); },
          py::arg("alphabet") = "aa",
          "Symbols in code order for an alphabet; custom matrices must follow this order.");
    m.def("amino_acids", [] { return alphabet_symbols(Alphabet::AminoAcid); },
          "The amino-acid symbol order used by the built-in matrices and custom AA matrices.");

    py::class_<PyParams>(m, "SearchParams",
                         "Search scope and budget. Scope: max_subs/max_ins/max_dels (exact, "
                         "seqtm) and max_total_edits. Budget: max_penalty with an optional "
                         "matrix (identity/BLOSUM62/PAM250/PAM100/structural) and gap costs. engine is 'auto'|'seqtrie'|'seqtm', "
                         "mode is 'all'|'top'.")
        .def(py::init([](int max_subs, int max_ins, int max_dels, int max_total_edits,
                         long max_penalty, py::object matrix, int gap_open, int gap_extend,
                         std::string engine, std::string mode) {
                 PyParams p;
                 p.max_subs = max_subs; p.max_ins = max_ins; p.max_dels = max_dels;
                 p.max_total_edits = max_total_edits; p.max_penalty = max_penalty;
                 set_matrix(p, matrix); p.gap_open = gap_open; p.gap_extend = gap_extend;
                 parse_engine(engine); parse_mode(mode);  // validate eagerly
                 p.engine = std::move(engine); p.mode = std::move(mode);
                 return p;
             }),
             py::arg("max_subs") = 0, py::arg("max_ins") = 0, py::arg("max_dels") = 0,
             py::arg("max_total_edits") = 0, py::arg("max_penalty") = 0,
             py::arg("matrix") = "", py::arg("gap_open") = 1, py::arg("gap_extend") = 1,
             py::arg("engine") = "auto", py::arg("mode") = "all")
        .def_readwrite("max_subs", &PyParams::max_subs)
        .def_readwrite("max_ins", &PyParams::max_ins)
        .def_readwrite("max_dels", &PyParams::max_dels)
        .def_readwrite("max_total_edits", &PyParams::max_total_edits)
        .def_readwrite("max_penalty", &PyParams::max_penalty)
        .def_property(
            "matrix",
            [](const PyParams& p) -> py::object {
                if (p.matrix_obj) return py::cast(*p.matrix_obj);
                return py::cast(p.matrix);
            },
            [](PyParams& p, const py::object& m) { set_matrix(p, m); })
        .def_property(
            "pos_matrix",
            [](const PyParams& p) -> py::object {
                if (p.pos_matrix_obj) return py::cast(*p.pos_matrix_obj);
                return py::none();
            },
            [](PyParams& p, const py::object& m) {
                if (m.is_none()) p.pos_matrix_obj.reset();
                else if (py::isinstance<PositionalMatrix>(m)) p.pos_matrix_obj = m.cast<PositionalMatrix>();
                else throw py::type_error("pos_matrix must be a PositionalMatrix or None");
            })
        .def_readwrite("gap_open", &PyParams::gap_open)
        .def_readwrite("gap_extend", &PyParams::gap_extend)
        .def_property("engine", [](const PyParams& p) { return p.engine; },
                      [](PyParams& p, std::string v) { parse_engine(v); p.engine = std::move(v); })
        .def_property("mode", [](const PyParams& p) { return p.mode; },
                      [](PyParams& p, std::string v) { parse_mode(v); p.mode = std::move(v); });

    py::class_<Hit>(m, "Hit",
                    "A search result. Payload-agnostic: map ``ref_id`` back to your own "
                    "payload downstream. ``score`` is a non-negative penalty (0 == exact). "
                    "``n_subs``/``n_ins``/``n_dels`` are exact for the seqtm engine and 0 for "
                    "seqtrie. Iterable as ``(ref_id, score, n_subs, n_ins, n_dels)``.")
        .def_readonly("ref_id", &Hit::ref_id)
        .def_readonly("score", &Hit::score)
        .def_readonly("n_subs", &Hit::n_subs)
        .def_readonly("n_ins", &Hit::n_ins)
        .def_readonly("n_dels", &Hit::n_dels)
        .def("__iter__", [](const Hit& h) {
            return py::iter(py::make_tuple(h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels));
        })
        .def("__repr__", [](const Hit& h) {
            return "Hit(ref_id=" + std::to_string(h.ref_id) + ", score=" + std::to_string(h.score) +
                   ", n_subs=" + std::to_string(h.n_subs) + ", n_ins=" + std::to_string(h.n_ins) +
                   ", n_dels=" + std::to_string(h.n_dels) + ")";
        });

    py::class_<Alignment>(m, "Alignment",
                          "Global alignment of a query to a reference. ``ops`` has one char per "
                          "column: 'M' match, 'S' substitution, 'I' insertion, 'D' deletion.")
        .def_readonly("aligned_query", &Alignment::aligned_query)
        .def_readonly("aligned_ref", &Alignment::aligned_ref)
        .def_readonly("ops", &Alignment::ops)
        .def_readonly("score", &Alignment::score)
        .def("__repr__", [](const Alignment& a) {
            return "Alignment(score=" + std::to_string(a.score) + ", ops='" + a.ops + "')";
        });

    py::class_<Index>(m, "Index",
                      "Immutable search index over a set of reference sequences. Build once, "
                      "then query concurrently; reference id is the position in ``refs``.")
        .def_static(
            "build",
            [](std::vector<std::string> refs, const std::string& alphabet) {
                return Index::build(std::move(refs), parse_alphabet(alphabet));
            },
            py::arg("refs"), py::arg("alphabet") = "aa",
            "Build an index. ``alphabet`` is 'aa', 'nt', or 'iupac'. Raises ValueError "
            "on a symbol outside the alphabet.")
        .def("__len__", &Index::size)
        .def("ref_seq", [](const Index& i, uint32_t id) { return std::string(i.ref_seq(id)); },
             py::arg("ref_id"), "Return the reference sequence string for a reference id.")
        .def("search", &py_search, py::arg("query"), py::arg("params"),
             "Return all hits for one query within the scope/budget in ``params``.")
        .def("search_top", &py_search_top, py::arg("query"), py::arg("params"), py::arg("k") = 1,
             "Return up to ``k`` best (lowest-score) hits for one query.")
        .def("search_batch", &py_search_batch, py::arg("queries"), py::arg("params"),
             py::arg("threads") = 0,
             "Search many queries in parallel (releases the GIL). ``threads=0`` uses all "
             "cores. Returns one hit list per query, in input order.")
        .def("align", &py_align, py::arg("ref_id"), py::arg("query"), py::arg("params"),
             "Compute a global alignment between ``query`` and a reference, on demand.")
        .def("collisions_batch", &py_collisions_batch, py::arg("queries"), py::arg("params"),
             py::arg("threads") = 0,
             "Per-query count of seqtm collisions: how often a reference was re-reached via a "
             "different edit path during branch-and-bound (0 for seqtrie / substitution-only).")
        .def("save", &Index::save, py::arg("path"),
             "Serialize the index to a flat binary file for fast reload.")
        .def_static("load", &Index::load, py::arg("path"),
                    "Load an index previously written with save(); raises on a corrupt/old file.");

    m.def("pairwise_batch", &py_pairwise_batch, py::arg("a"), py::arg("b"), py::arg("params"),
          py::arg("alphabet") = "aa", py::arg("threads") = 0,
          "Batch-vs-batch search. Indexes the larger set internally and streams the smaller; "
          "results are a-major (one hit list per a[i]) with Hit.ref_id pointing into b.");

    py::class_<Candidate>(m, "Candidate",
                          "A seed-and-gather hit: peptide_id, shared_kmers (distinct query k-mers "
                          "that hit it), best_score. Iterable as (peptide_id, shared_kmers, best_score).")
        .def_readonly("peptide_id", &Candidate::peptide_id)
        .def_readonly("shared_kmers", &Candidate::shared_kmers)
        .def_readonly("best_score", &Candidate::best_score)
        .def("__iter__", [](const Candidate& c) {
            return py::iter(py::make_tuple(c.peptide_id, c.shared_kmers, c.best_score));
        })
        .def("__repr__", [](const Candidate& c) {
            return "Candidate(peptide_id=" + std::to_string(c.peptide_id) +
                   ", shared_kmers=" + std::to_string(c.shared_kmers) +
                   ", best_score=" + std::to_string(c.best_score) + ")";
        });

    py::class_<KmerIndex>(m, "KmerIndex",
                          "Seed-and-extend k-mer index for homology. Build from per-peptide k-mer "
                          "lists (anchor-masked upstream) + optional allele tags; seed_and_gather "
                          "fuzzy-matches query k-mers and merges posting lists into ranked "
                          "candidates entirely in C++ (GIL released).")
        .def_static(
            "build",
            [](const std::vector<std::vector<std::string>>& kmers, const std::string& alphabet,
               const std::vector<uint32_t>& allele_ids) {
                return KmerIndex::build(kmers, parse_alphabet(alphabet), allele_ids);
            },
            py::arg("kmers_per_peptide"), py::arg("alphabet") = "aa",
            py::arg("allele_ids") = std::vector<uint32_t>{})
        .def("num_peptides", &KmerIndex::num_peptides)
        .def("num_kmers", &KmerIndex::num_kmers)
        .def("__len__", &KmerIndex::num_peptides)
        .def(
            "seed_and_gather",
            [](const KmerIndex& ki, const std::vector<std::vector<std::string>>& qk,
               const PyParams& pp, uint32_t min_shared, int64_t allele_filter, int threads) {
                auto mat = make_matrix(pp, ki.alphabet());
                SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
                std::vector<std::vector<Candidate>> res;
                {
                    py::gil_scoped_release release;
                    res = ki.seed_and_gather(qk, cp, min_shared, allele_filter, threads);
                }
                py::list out(res.size());
                for (size_t i = 0; i < res.size(); ++i) {
                    py::list inner(res[i].size());
                    for (size_t j = 0; j < res[i].size(); ++j) inner[j] = py::cast(res[i][j]);
                    out[i] = inner;
                }
                return out;
            },
            py::arg("query_kmers"), py::arg("params"), py::arg("min_shared") = 1,
            py::arg("allele_filter") = -1, py::arg("threads") = 0,
            "For each query (its k-mer list) return ranked Candidates with >= min_shared shared "
            "k-mers; allele_filter >= 0 restricts to that allele tag.")
        .def("save", &KmerIndex::save, py::arg("path"))
        .def_static("load", &KmerIndex::load, py::arg("path"));
}
