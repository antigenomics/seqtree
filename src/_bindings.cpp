#include "seqtree/seqtree.hpp"

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
    std::string matrix;          // "" => unit cost
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

// Returns nullopt for unit cost; throws for an unknown name or alphabet mismatch.
std::optional<SubstitutionMatrix> make_matrix(const std::string& name, Alphabet a) {
    if (name.empty()) return std::nullopt;
    std::string l = lower(name);
    if (l == "blosum62") {
        if (a != Alphabet::AminoAcid)
            throw py::value_error("BLOSUM62 requires the amino-acid alphabet");
        return SubstitutionMatrix::blosum62();
    }
    throw py::value_error("unknown matrix '" + name + "' (use '' or 'BLOSUM62')");
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
    return p;
}

py::list hits_to_list(const std::vector<Hit>& hits) {
    py::list out(hits.size());
    for (size_t i = 0; i < hits.size(); ++i) out[i] = py::cast(hits[i]);
    return out;
}

py::list py_search(Index& idx, const std::string& q, const PyParams& pp) {
    auto mat = make_matrix(pp.matrix, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    Searcher s(idx);
    return hits_to_list(s.search(q, cp));
}

py::list py_search_top(Index& idx, const std::string& q, const PyParams& pp, int k) {
    auto mat = make_matrix(pp.matrix, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    cp.mode = Mode::TopHit;
    cp.max_hits = uint32_t(k < 1 ? 1 : k);
    Searcher s(idx);
    return hits_to_list(s.search(q, cp));
}

py::list py_search_batch(const Index& idx, const std::vector<std::string>& queries,
                         const PyParams& pp, int threads) {
    auto mat = make_matrix(pp.matrix, idx.alphabet());
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

py::list py_pairwise_batch(const std::vector<std::string>& a, const std::vector<std::string>& b,
                           const PyParams& pp, const std::string& alphabet, int threads) {
    Alphabet alph = parse_alphabet(alphabet);
    auto mat = make_matrix(pp.matrix, alph);
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
    auto mat = make_matrix(pp.matrix, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    return idx.align(q, ref_id, cp);
}

}  // namespace

PYBIND11_MODULE(_core, m) {
    m.doc() = "seqtree: fuzzy biological-sequence search (C++ core)";

    py::class_<PyParams>(m, "SearchParams",
                         "Search scope and budget. Scope: max_subs/max_ins/max_dels (exact, "
                         "seqtm) and max_total_edits. Budget: max_penalty with an optional "
                         "matrix ('BLOSUM62') and gap costs. engine is 'auto'|'seqtrie'|'seqtm', "
                         "mode is 'all'|'top'.")
        .def(py::init([](int max_subs, int max_ins, int max_dels, int max_total_edits,
                         long max_penalty, std::string matrix, int gap_open, int gap_extend,
                         std::string engine, std::string mode) {
                 PyParams p;
                 p.max_subs = max_subs; p.max_ins = max_ins; p.max_dels = max_dels;
                 p.max_total_edits = max_total_edits; p.max_penalty = max_penalty;
                 p.matrix = std::move(matrix); p.gap_open = gap_open; p.gap_extend = gap_extend;
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
        .def_readwrite("matrix", &PyParams::matrix)
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
             "Compute a global alignment between ``query`` and a reference, on demand.");

    m.def("pairwise_batch", &py_pairwise_batch, py::arg("a"), py::arg("b"), py::arg("params"),
          py::arg("alphabet") = "aa", py::arg("threads") = 0,
          "Batch-vs-batch search. Indexes the larger set internally and streams the smaller; "
          "results are a-major (one hit list per a[i]) with Hit.ref_id pointing into b.");
}
