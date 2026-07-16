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
    std::string mode = "all";    // all | top
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
    throw py::value_error("unknown mode '" + m + "' (use 'all' or 'top')");
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
constexpr const char* kMatrixNames =
    "'identity', 'BLOSUM62', 'BLOSUM45', 'BLOSUM80', 'PAM250', 'PAM100', 'structural'";

bool is_matrix_name(const std::string& l) {
    return l.empty() || l == "identity" || l == "blosum62" || l == "blosum45"
        || l == "blosum80" || l == "pam250" || l == "pam100" || l == "structural";
}

SubstitutionMatrix named_matrix(const std::string& l, Alphabet a) {
    if (l == "identity") return SubstitutionMatrix::unit(Codec(a).size());
    if (a != Alphabet::AminoAcid)
        throw py::value_error(l + " requires the amino-acid alphabet");
    if (l == "blosum62") return SubstitutionMatrix::blosum62();
    if (l == "blosum45") return SubstitutionMatrix::blosum45();
    if (l == "blosum80") return SubstitutionMatrix::blosum80();
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

// Owns an N*K int32 block and lends it out through the CPython buffer protocol. seqtree has no
// runtime dependencies, so we cannot hand back a numpy array; a buffer lets numpy (or plain
// memoryview) wrap the same memory with no copy.
struct ScoreMatrix {
    std::vector<int32_t> data;
    size_t rows = 0, cols = 0;
};

AlignMode parse_align_mode(const std::string& m) {
    std::string l = lower(m);
    if (l == "global" || l == "nw" || l == "needleman-wunsch") return AlignMode::Global;
    if (l == "local" || l == "sw" || l == "smith-waterman") return AlignMode::Local;
    throw py::value_error("unknown mode '" + m + "' (use 'global' or 'local')");
}

ScoreMatrix py_gapblock_matrix(const std::vector<std::string>& queries,
                               const std::vector<std::string>& refs, const std::string& alphabet,
                               const std::optional<SubstitutionMatrix>& matrix, int32_t gap_open,
                               int32_t gap_extend, const std::vector<int32_t>& prior,
                               uint32_t prior_width, int threads) {
    Alphabet alph = parse_alphabet(alphabet);
    ScoreMatrix out;
    out.rows = queries.size();
    out.cols = refs.size();
    {
        py::gil_scoped_release release;
        out.data = gapblock_matrix(queries, refs, alph, matrix ? &*matrix : nullptr, gap_open,
                                   gap_extend, prior, prior_width, threads);
    }
    return out;
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
        .def_static("blosum45", &SubstitutionMatrix::blosum45)
        .def_static("blosum80", &SubstitutionMatrix::blosum80)
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
            "non-negative penalties via the Gram / squared-distance transform "
            "s[a,a] + s[b,b] - 2*s[a,b] (clamped at 0). Row/column order must match the "
            "target alphabet's symbol order (see ``amino_acids()``).")
        .def("size", &SubstitutionMatrix::size)
        .def(
            "similarity",
            [](const SubstitutionMatrix& self, const std::string& a, const std::string& b) {
                if (a.size() != 1 || b.size() != 1)
                    throw py::value_error("similarity() takes two single amino-acid characters");
                static const std::string aa = alphabet_symbols(Alphabet::AminoAcid);
                auto ia = aa.find(a[0]), ib = aa.find(b[0]);
                if (ia == std::string::npos || ib == std::string::npos)
                    throw py::value_error("unknown amino acid; expected one of " + aa);
                // Without this the 24-symbol AA index runs off the end of a smaller matrix
                // (e.g. unit(4) for nucleotides) and returns heap garbage. penalty() has always
                // checked; similarity() must too.
                if (ia >= self.size() || ib >= self.size())
                    throw py::value_error("residue out of range for this matrix's alphabet");
                return self.similarity(uint8_t(ia), uint8_t(ib));
            },
            py::arg("a"), py::arg("b"),
            "Raw log-odds similarity (signed). penalty() is the non-negative Gram "
            "transform of this; the transform is lossy, so both are kept.")
        .def("scale", &SubstitutionMatrix::scale,
             "Median penalty over all mismatched symbol pairs -- this matrix's natural unit. "
             "Gap costs must be on this scale: BLOSUM62 has scale() == 14, so the default "
             "gap_open of 1 makes gaps ~14x cheaper than substitutions and the aligner gaps "
             "rather than substitutes. Use ``gap_open = 1-2 * m.scale()``.")
        .def(
            "penalty",
            [](const SubstitutionMatrix& self, const std::string& a, const std::string& b) {
                if (a.size() != 1 || b.size() != 1)
                    throw py::value_error("penalty() takes two single amino-acid characters");
                static const std::string aa = alphabet_symbols(Alphabet::AminoAcid);
                auto ia = aa.find(a[0]), ib = aa.find(b[0]);
                if (ia == std::string::npos || ib == std::string::npos)
                    throw py::value_error("unknown amino acid; expected one of " + aa);
                if (ia >= self.size() || ib >= self.size())
                    throw py::value_error("residue out of range for this matrix's alphabet");
                return self.penalty(uint8_t(ia), uint8_t(ib));
            },
            py::arg("a"), py::arg("b"),
            "Gram-distance substitution penalty between two amino acids: 0 when identical, "
            "larger when more dissimilar (s(a,a)+s(b,b)-2 s(a,b)). Characters use the "
            "``amino_acids()`` order.")
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
                    "position. len(weights) is the frame width. NOTE: penalty(a, a) == 0 for "
                    "every base matrix, so a weight scales MISMATCH cost only -- it is a "
                    "mismatch-tolerance profile, not an information/match weighting.")
        .def_static("from_tables", &PositionalMatrix::from_tables,
                    py::arg("size"), py::arg("width"), py::arg("data"),
                    py::arg("masked") = std::vector<uint8_t>{},
                    "Full per-position PSSM. ``data`` is row-major [width][size][size]; "
                    "``masked`` is an optional length-``width`` flag array (non-zero == free "
                    "position). Use this to give different regions different matrices, e.g. a "
                    "germline-flank matrix and an N-region core matrix in one frame.")
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

    py::class_<ScoreMatrix>(m, "ScoreMatrix", py::buffer_protocol(),
                            "A read-only (n_queries, n_refs) int32 penalty matrix, row-major. "
                            "Exposes the buffer protocol, so ``numpy.asarray(sm)`` and "
                            "``memoryview(sm)`` both wrap it without copying. Index it with "
                            "``sm[i, k]`` or pull one row with ``sm.row(i)``.")
        .def_buffer([](ScoreMatrix& s) {
            return py::buffer_info(s.data.data(), sizeof(int32_t),
                                   py::format_descriptor<int32_t>::format(), 2, {s.rows, s.cols},
                                   {sizeof(int32_t) * s.cols, sizeof(int32_t)});
        })
        .def_property_readonly("shape",
                               [](const ScoreMatrix& s) { return py::make_tuple(s.rows, s.cols); })
        .def("__len__", [](const ScoreMatrix& s) { return s.rows; })
        .def(
            "row",
            [](const ScoreMatrix& s, size_t i) {
                if (i >= s.rows) throw py::index_error("row out of range");
                return std::vector<int32_t>(s.data.begin() + i * s.cols,
                                            s.data.begin() + (i + 1) * s.cols);
            },
            py::arg("i"), "Row i as a list of penalties, one per reference.")
        .def("__getitem__",
             [](const ScoreMatrix& s, std::pair<size_t, size_t> ik) {
                 if (ik.first >= s.rows || ik.second >= s.cols)
                     throw py::index_error("index out of range");
                 return s.data[ik.first * s.cols + ik.second];
             })
        .def("__repr__", [](const ScoreMatrix& s) {
            return "ScoreMatrix(" + std::to_string(s.rows) + ", " + std::to_string(s.cols) + ")";
        });

    m.def(
        "align_score",
        [](const std::string& q, const std::string& r, const SubstitutionMatrix& mat,
           const std::string& mode, int32_t gap_open, int32_t gap_extend,
           const std::string& alphabet) {
            return align_score(q, r, mat, parse_alphabet(alphabet), parse_align_mode(mode),
                               gap_open, gap_extend);
        },
        py::arg("query"), py::arg("ref"), py::arg("matrix"), py::arg("mode") = "global",
        py::arg("gap_open") = 11, py::arg("gap_extend") = 1, py::arg("alphabet") = "aa",
        "Optimal similarity score. 'global' is Needleman-Wunsch, 'local' Smith-Waterman; "
        "gap_open == gap_extend gives linear gaps. Gap costs are positive magnitudes.");

    m.def(
        "align_pair",
        [](const std::string& q, const std::string& r, const SubstitutionMatrix& mat,
           const std::string& mode, int32_t gap_open, int32_t gap_extend,
           const std::string& alphabet) {
            return align_pair(q, r, mat, parse_alphabet(alphabet), parse_align_mode(mode),
                              gap_open, gap_extend);
        },
        py::arg("query"), py::arg("ref"), py::arg("matrix"), py::arg("mode") = "global",
        py::arg("gap_open") = 11, py::arg("gap_extend") = 1, py::arg("alphabet") = "aa",
        "As align_score, but also returns the aligned strings and ops (Alignment.score is the "
        "similarity, not a penalty).");

    m.def(
        "align_score_matrix",
        [](const std::vector<std::string>& q, const std::vector<std::string>& r,
           const SubstitutionMatrix& mat, const std::string& mode, int32_t gap_open,
           int32_t gap_extend, const std::string& alphabet, int threads) {
            Alphabet a = parse_alphabet(alphabet);
            AlignMode md = parse_align_mode(mode);
            ScoreMatrix out;
            out.rows = q.size();
            out.cols = r.size();
            {
                py::gil_scoped_release release;
                out.data = align_score_matrix(q, r, mat, a, md, gap_open, gap_extend, threads);
            }
            return out;
        },
        py::arg("queries"), py::arg("refs"), py::arg("matrix"), py::arg("mode") = "global",
        py::arg("gap_open") = 11, py::arg("gap_extend") = 1, py::arg("alphabet") = "aa",
        py::arg("threads") = 0,
        "Dense (n_queries, n_refs) similarity matrix, GIL released.");

    m.def(
        "align_dist_matrix",
        [](const std::vector<std::string>& q, const std::vector<std::string>& r,
           const SubstitutionMatrix& mat, const std::string& mode, int32_t gap_open,
           int32_t gap_extend, const std::string& alphabet, int threads) {
            Alphabet a = parse_alphabet(alphabet);
            AlignMode md = parse_align_mode(mode);
            ScoreMatrix out;
            out.rows = q.size();
            out.cols = r.size();
            {
                py::gil_scoped_release release;
                out.data = align_dist_matrix(q, r, mat, a, md, gap_open, gap_extend, threads);
            }
            return out;
        },
        py::arg("queries"), py::arg("refs"), py::arg("matrix"), py::arg("mode") = "global",
        py::arg("gap_open") = 11, py::arg("gap_extend") = 1, py::arg("alphabet") = "aa",
        py::arg("threads") = 0,
        "Dense (n_queries, n_refs) distance matrix d = s(a,a) + s(b,b) - 2*s(a,b), the "
        "sequence-level Gram transform of the alignment scores. Non-negative, zero on identity.");

    m.def("gapblock_matrix", &py_gapblock_matrix, py::arg("queries"), py::arg("refs"),
          py::arg("alphabet") = "aa", py::arg("matrix") = std::nullopt, py::arg("gap_open") = 1,
          py::arg("gap_extend") = 1, py::arg("prior") = std::vector<int32_t>{},
          py::arg("prior_width") = 0, py::arg("threads") = 0,
          "Exhaustive single-gap-block penalties for every (query, ref) pair, GIL released. "
          "`prior` is the gap prior flattened to [m][d][i]; see seqtree.gapblock.score_matrix, "
          "which builds it for you.");

    m.def("hamming", &hamming, py::arg("a"), py::arg("b"),
          "Hamming distance: the number of positions at which two EQUAL-length sequences differ. "
          "Raises ValueError on a length mismatch. Case-sensitive, byte-for-byte.");

    m.def("levenshtein", &levenshtein, py::arg("a"), py::arg("b"),
          "Levenshtein (edit) distance: the minimum number of single-character insertions, "
          "deletions, and substitutions to turn `a` into `b`, each cost 1. Case-sensitive.");

    m.def(
        "hamming_matrix",
        [](const std::vector<std::string>& a, const std::vector<std::string>& b, int threads) {
            ScoreMatrix out;
            out.rows = a.size();
            out.cols = b.size();
            {
                py::gil_scoped_release release;
                out.data = hamming_matrix(a, b, threads);
            }
            return out;
        },
        py::arg("a"), py::arg("b"), py::arg("threads") = 0,
        "Dense (len(a), len(b)) int32 Hamming-distance matrix, GIL released. Raises ValueError "
        "if any pair has mismatched lengths.");

    m.def(
        "levenshtein_matrix",
        [](const std::vector<std::string>& a, const std::vector<std::string>& b, int threads) {
            ScoreMatrix out;
            out.rows = a.size();
            out.cols = b.size();
            {
                py::gil_scoped_release release;
                out.data = levenshtein_matrix(a, b, threads);
            }
            return out;
        },
        py::arg("a"), py::arg("b"), py::arg("threads") = 0,
        "Dense (len(a), len(b)) int32 Levenshtein-distance matrix, GIL released.");

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
