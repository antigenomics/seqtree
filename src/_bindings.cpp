#include "seqtree/seqtree.hpp"
#include "seqtree/kmer_index.hpp"
#include "seqtree/text_index.hpp"

#include <nanobind/nanobind.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/pair.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/string_view.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/unique_ptr.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

namespace nb = nanobind;
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
    throw nb::value_error(("unknown alphabet '" + a + "' (use 'aa', 'nt', or 'iupac')").c_str());
}

const char* alphabet_name(Alphabet a) {
    switch (a) {
        case Alphabet::Nucleotide:      return "nt";
        case Alphabet::NucleotideIUPAC: return "iupac";
        default:                        return "aa";
    }
}

Engine parse_engine(const std::string& e) {
    std::string l = lower(e);
    if (l == "auto") return Engine::Auto;
    if (l == "seqtrie") return Engine::SeqTrie;
    if (l == "seqtm") return Engine::SeqTm;
    throw nb::value_error(("unknown engine '" + e + "' (use 'auto', 'seqtrie', or 'seqtm')").c_str());
}

Mode parse_mode(const std::string& m) {
    std::string l = lower(m);
    if (l == "all") return Mode::AllHits;
    if (l == "top") return Mode::TopHit;
    throw nb::value_error(("unknown mode '" + m + "' (use 'all' or 'top')").c_str());
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
        throw nb::value_error((l + " requires the amino-acid alphabet").c_str());
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
            throw nb::value_error("matrix size does not match the alphabet");
        return *pp.matrix_obj;
    }
    if (pp.matrix.empty()) return std::nullopt;
    std::string l = lower(pp.matrix);
    if (!is_matrix_name(l))
        throw nb::value_error(("unknown matrix '" + pp.matrix + "' (use '', " +
                               kMatrixNames + ", or a SubstitutionMatrix)").c_str());
    return named_matrix(l, a);
}

// Accept either a builtin name (see kMatrixNames) or a SubstitutionMatrix.
void set_matrix(PyParams& p, const nb::object& m) {
    p.matrix.clear();
    p.matrix_obj.reset();
    if (m.is_none()) return;
    if (nb::isinstance<nb::str>(m)) {
        std::string name = nb::cast<std::string>(m);
        if (!is_matrix_name(lower(name)))
            throw nb::value_error(("unknown matrix '" + name + "' (use '', " +
                                   kMatrixNames + ", or a SubstitutionMatrix)").c_str());
        p.matrix = std::move(name);
    } else if (nb::isinstance<SubstitutionMatrix>(m)) {
        p.matrix_obj = nb::cast<SubstitutionMatrix>(m);
    } else {
        throw nb::type_error("matrix must be a name string or a SubstitutionMatrix");
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

nb::list hits_to_list(const std::vector<Hit>& hits) {
    nb::list out;
    for (const Hit& h : hits) out.append(nb::cast(h));
    return out;
}

nb::list py_search(Index& idx, const std::string& q, const PyParams& pp) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    Searcher s(idx);
    return hits_to_list(s.search(q, cp));
}

nb::list py_search_top(Index& idx, const std::string& q, const PyParams& pp, int k) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    cp.mode = Mode::TopHit;
    cp.max_hits = uint32_t(k < 1 ? 1 : k);
    Searcher s(idx);
    return hits_to_list(s.search(q, cp));
}

nb::list py_search_batch(const Index& idx, const std::vector<std::string>& queries,
                         const PyParams& pp, int threads) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    std::vector<std::vector<Hit>> results;
    {
        nb::gil_scoped_release release;  // pure C++ region, no Python objects touched
        results = idx.search_batch(queries, cp, threads);
    }
    nb::list out;
    for (const auto& hits : results) out.append(hits_to_list(hits));
    return out;
}

std::vector<uint64_t> py_collisions_batch(const Index& idx, const std::vector<std::string>& queries,
                                          const PyParams& pp, int threads) {
    auto mat = make_matrix(pp, idx.alphabet());
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    nb::gil_scoped_release release;
    return idx.collisions_batch(queries, cp, threads);
}

nb::list py_pairwise_batch(const std::vector<std::string>& a, const std::vector<std::string>& b,
                           const PyParams& pp, const std::string& alphabet, int threads) {
    Alphabet alph = parse_alphabet(alphabet);
    auto mat = make_matrix(pp, alph);
    SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
    std::vector<std::vector<Hit>> results;
    {
        nb::gil_scoped_release release;
        results = pairwise_batch(a, b, alph, cp, threads);
    }
    nb::list out;
    for (const auto& hits : results) out.append(hits_to_list(hits));
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

// nanobind dropped pybind11's def_buffer, so the buffer slots are wired by hand and handed to
// nb::class_ as type_slots. Keeping the protocol (rather than moving to nb::ndarray) is
// deliberate: seqtree declares no runtime dependencies, and both `numpy.asarray(x)` and
// `memoryview(x)` are documented, tested guarantees.
//
// `shape` and `strides` must outlive the getbuffer call, so both live in one heap array that
// travels in `internal` and is freed by the release slot: [0, ndim) is shape, [ndim, 2*ndim)
// is strides. `fmt` must outlive it too, so it points at storage the exporting object owns.
void fill_view(PyObject* obj, Py_buffer* view, int flags, void* buf, int ndim,
               const Py_ssize_t* shape, const Py_ssize_t* strides, Py_ssize_t itemsize,
               const char* fmt) {
    auto* dims = new Py_ssize_t[2 * size_t(ndim)];
    std::copy(shape, shape + ndim, dims);
    std::copy(strides, strides + ndim, dims + ndim);
    Py_ssize_t len = itemsize;
    for (int d = 0; d < ndim; ++d) len *= shape[d];
    view->buf = buf;
    view->len = len;
    view->readonly = 1;
    view->itemsize = itemsize;
    view->format = (flags & PyBUF_FORMAT) == PyBUF_FORMAT ? const_cast<char*>(fmt) : nullptr;
    view->ndim = ndim;
    view->shape = dims;
    view->strides = dims + ndim;
    view->suboffsets = nullptr;
    view->internal = dims;
    view->obj = Py_NewRef(obj);  // keeps the exporter (and its data) alive for the view
}

void release_view(PyObject*, Py_buffer* view) {
    delete[] static_cast<Py_ssize_t*>(view->internal);
    view->internal = nullptr;
}

int score_matrix_getbuffer(PyObject* obj, Py_buffer* view, int flags) {
    const ScoreMatrix* s = nb::inst_ptr<ScoreMatrix>(obj);
    const Py_ssize_t shape[2] = {Py_ssize_t(s->rows), Py_ssize_t(s->cols)};
    const Py_ssize_t strides[2] = {Py_ssize_t(sizeof(int32_t) * s->cols),
                                   Py_ssize_t(sizeof(int32_t))};
    fill_view(obj, view, flags, const_cast<int32_t*>(s->data.data()), 2, shape, strides,
              Py_ssize_t(sizeof(int32_t)), "i");
    return 0;
}

PyType_Slot kScoreMatrixSlots[] = {
    { Py_bf_getbuffer, (void*)score_matrix_getbuffer },
    { Py_bf_releasebuffer, (void*)release_view },
    { 0, nullptr }
};

// A read-only 1-D typed window onto memory a TextResult owns. Same reasoning as ScoreMatrix:
// seqtree has no runtime dependencies, so the way out is the buffer protocol, not an ndarray.
// `owner` is a strong reference to the TextResult, so a view outlives the dict it came from.
struct ArrayView {
    nb::object  owner;
    const void* data = nullptr;
    size_t      n = 0;
    size_t      itemsize = 0;
    // NUL-terminated struct format code. Per instance, not a shared scratch buffer: two live
    // views of different dtypes must each keep reporting their own format.
    char        fmt[2] = {'i', 0};
};

int array_view_getbuffer(PyObject* obj, Py_buffer* view, int flags) {
    const ArrayView* a = nb::inst_ptr<ArrayView>(obj);
    const Py_ssize_t shape[1] = {Py_ssize_t(a->n)};
    const Py_ssize_t strides[1] = {Py_ssize_t(a->itemsize)};
    fill_view(obj, view, flags, const_cast<void*>(a->data), 1, shape, strides,
              Py_ssize_t(a->itemsize), a->fmt);
    return 0;
}

PyType_Slot kArrayViewSlots[] = {
    { Py_bf_getbuffer, (void*)array_view_getbuffer },
    { Py_bf_releasebuffer, (void*)release_view },
    { 0, nullptr }
};

template <class T>
ArrayView view_of(const nb::object& owner, const std::vector<T>& v, char fmt) {
    return ArrayView{owner, v.data(), v.size(), sizeof(T), {fmt, 0}};
}

// One hit as a Python object. Built only for the query actually asked for, never for the whole
// batch -- the flat arrays are what a 445k-query run should read.
struct PyTextHit {
    uint32_t ref_id = 0, offset = 0;
    uint16_t n_subs = 0, n_ins = 0, n_dels = 0, length = 0;
    int32_t  score = 0;
    std::vector<std::tuple<uint16_t, std::string, std::string>> mismatches;
};

nb::dict text_result_arrays(const nb::object& self) {
    const TextResult& r = nb::cast<const TextResult&>(self);
    nb::dict d;
    d["query_begin"] = view_of(self, r.query_begin, 'I');
    d["ref_id"] = view_of(self, r.ref_id, 'I');
    d["offset"] = view_of(self, r.offset, 'I');
    d["n_subs"] = view_of(self, r.n_subs, 'H');
    d["n_ins"] = view_of(self, r.n_ins, 'H');
    d["n_dels"] = view_of(self, r.n_dels, 'H');
    d["length"] = view_of(self, r.length, 'H');
    d["score"] = view_of(self, r.score, 'i');
    d["mm_begin"] = view_of(self, r.mm_begin, 'I');
    d["mm_pos"] = view_of(self, r.mm_pos, 'H');
    d["mm_query_aa"] = view_of(self, r.mm_query_aa, 'c');
    d["mm_text_aa"] = view_of(self, r.mm_text_aa, 'c');
    d["truncated"] = view_of(self, r.truncated, 'B');
    d["group_begin"] = view_of(self, r.group_begin, 'I');
    d["group_id"] = view_of(self, r.group_id, 'I');
    d["group_min_subs"] = view_of(self, r.group_min_subs, 'H');
    d["group_n_hits"] = view_of(self, r.group_n_hits, 'I');
    return d;
}

AlignMode parse_align_mode(const std::string& m) {
    std::string l = lower(m);
    if (l == "global" || l == "nw" || l == "needleman-wunsch") return AlignMode::Global;
    if (l == "local" || l == "sw" || l == "smith-waterman") return AlignMode::Local;
    throw nb::value_error(("unknown mode '" + m + "' (use 'global' or 'local')").c_str());
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
        nb::gil_scoped_release release;
        out.data = gapblock_matrix(queries, refs, alph, matrix ? &*matrix : nullptr, gap_open,
                                   gap_extend, prior, prior_width, threads);
    }
    return out;
}

}  // namespace

NB_MODULE(_core, m) {
    m.doc() = "seqtree: fuzzy biological-sequence search (C++ core)";

    nb::class_<SubstitutionMatrix>(m, "SubstitutionMatrix",
                                   "Non-negative substitution penalties (penalty(a,a)==0). Build a "
                                   "named builtin (``blosum62``/``pam250``/``pam100``/``structural``, "
                                   "or ``unit`` for identity) or a custom one from a similarity grid "
                                   "whose row/column order matches ``amino_acids()`` (or "
                                   "``alphabet_symbols(alphabet)``).")
        .def_static("blosum62", &SubstitutionMatrix::blosum62,
                    "BLOSUM62, the general-purpose default. Amino acids only; "
                    "``scale()`` is 14.")
        .def_static("blosum45", &SubstitutionMatrix::blosum45,
                    "BLOSUM45. More permissive than BLOSUM62 -- built from more divergent "
                    "alignments, so distant substitutions cost less. Amino acids only.")
        .def_static("blosum80", &SubstitutionMatrix::blosum80,
                    "BLOSUM80. Stricter than BLOSUM62 -- built from closer alignments, so it "
                    "separates near-identical sequences more sharply. Amino acids only.")
        .def_static("pam250", &SubstitutionMatrix::pam250,
                    "PAM250, an evolutionary model tuned for distant homology. "
                    "Amino acids only.")
        .def_static("pam100", &SubstitutionMatrix::pam100,
                    "PAM100, the same model at a shorter evolutionary distance -- stricter "
                    "than PAM250. Amino acids only.")
        .def_static("structural", &SubstitutionMatrix::structural,
                    "A structure-derived matrix, grouping residues by side-chain shape and "
                    "charge rather than by observed substitution frequency. Amino acids only.")
        .def_static("unit", &SubstitutionMatrix::unit, nb::arg("size"),
                    "Identity: every mismatch costs 1, every match 0, so a score is a plain "
                    "edit count. ``size`` is the alphabet size (24 for ``aa``, 4 for ``nt``); "
                    "use ``len(alphabet_symbols(alphabet))``. The only matrix valid for "
                    "nucleotides.")
        .def_static(
            "from_similarity",
            [](const std::vector<std::vector<int32_t>>& grid) {
                size_t n = grid.size();
                if (n == 0 || n > 32) throw nb::value_error("matrix size must be 1..32");
                std::vector<int32_t> flat;
                flat.reserve(n * n);
                for (const auto& row : grid) {
                    if (row.size() != n) throw nb::value_error("similarity matrix must be square");
                    flat.insert(flat.end(), row.begin(), row.end());
                }
                return SubstitutionMatrix::from_similarity(uint8_t(n), flat.data());
            },
            nb::arg("grid"),
            "Build from a square similarity grid (higher == more similar), converted to "
            "non-negative penalties via the Gram / squared-distance transform "
            "s[a,a] + s[b,b] - 2*s[a,b] (clamped at 0). Row/column order must match the "
            "target alphabet's symbol order (see ``amino_acids()``).")
        .def("size", &SubstitutionMatrix::size,
             "Alphabet size this matrix scores over -- 24 for amino acids, 4 for nucleotides. A matrix must match the alphabet of the index it is used with.")
        .def(
            "similarity",
            [](const SubstitutionMatrix& self, const std::string& a, const std::string& b) {
                if (a.size() != 1 || b.size() != 1)
                    throw nb::value_error("similarity() takes two single amino-acid characters");
                static const std::string aa = alphabet_symbols(Alphabet::AminoAcid);
                auto ia = aa.find(a[0]), ib = aa.find(b[0]);
                if (ia == std::string::npos || ib == std::string::npos)
                    throw nb::value_error(("unknown amino acid; expected one of " + aa).c_str());
                // Without this the 24-symbol AA index runs off the end of a smaller matrix
                // (e.g. unit(4) for nucleotides) and returns heap garbage. penalty() has always
                // checked; similarity() must too.
                if (ia >= self.size() || ib >= self.size())
                    throw nb::value_error("residue out of range for this matrix's alphabet");
                return self.similarity(uint8_t(ia), uint8_t(ib));
            },
            nb::arg("a"), nb::arg("b"),
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
                    throw nb::value_error("penalty() takes two single amino-acid characters");
                static const std::string aa = alphabet_symbols(Alphabet::AminoAcid);
                auto ia = aa.find(a[0]), ib = aa.find(b[0]);
                if (ia == std::string::npos || ib == std::string::npos)
                    throw nb::value_error(("unknown amino acid; expected one of " + aa).c_str());
                if (ia >= self.size() || ib >= self.size())
                    throw nb::value_error("residue out of range for this matrix's alphabet");
                return self.penalty(uint8_t(ia), uint8_t(ib));
            },
            nb::arg("a"), nb::arg("b"),
            "Gram-distance substitution penalty between two amino acids: 0 when identical, "
            "larger when more dissimilar (s(a,a)+s(b,b)-2 s(a,b)). Characters use the "
            "``amino_acids()`` order.")
        .def("__repr__", [](const SubstitutionMatrix& s) {
            return "SubstitutionMatrix(size=" + std::to_string(s.size()) + ")";
        });

    nb::class_<PositionalMatrix>(m, "PositionalMatrix",
                                 "Per-position penalties pen(pos, a, b) over a fixed frame width. "
                                 "Build from a base SubstitutionMatrix and per-position integer "
                                 "weights: weight 0 masks the position (free, not counted as a "
                                 "substitution -- e.g. an anchor); >1 up-weights it (e.g. a TCR "
                                 "hotspot). Used on the seqtm Hamming path when width == query "
                                 "length.")
        .def_static("from_weights", &PositionalMatrix::from_weights,
                    nb::arg("base"), nb::arg("weights"),
                    "pen[pos][a][b] = weights[pos] * base.penalty(a, b); weight 0 masks the "
                    "position. len(weights) is the frame width. NOTE: penalty(a, a) == 0 for "
                    "every base matrix, so a weight scales MISMATCH cost only -- it is a "
                    "mismatch-tolerance profile, not an information/match weighting.")
        .def_static("from_tables", &PositionalMatrix::from_tables,
                    nb::arg("size"), nb::arg("width"), nb::arg("data"),
                    nb::arg("masked") = std::vector<uint8_t>{},
                    "Full per-position PSSM. ``data`` is row-major [width][size][size]; "
                    "``masked`` is an optional length-``width`` flag array (non-zero == free "
                    "position). Use this to give different regions different matrices, e.g. a "
                    "germline-flank matrix and an N-region core matrix in one frame.")
        .def("size", &PositionalMatrix::size, "Alphabet size this matrix scores over.")
        .def("width", &PositionalMatrix::width,
             "Frame width in positions. A query must be exactly this long for the positional "
             "path to apply.")
        .def("masked", &PositionalMatrix::masked, nb::arg("pos"),
             "True if ``pos`` is free: mismatches there cost nothing and are not counted as "
             "substitutions. This is how an anchor position is excluded from the ball.")
        .def("penalty", &PositionalMatrix::penalty, nb::arg("pos"), nb::arg("a"), nb::arg("b"),
             "Penalty for substituting symbol ``a`` with ``b`` at position ``pos``. Symbols are "
             "single-character strings. Always >= 0, and 0 when ``a == b`` or ``pos`` is "
             "masked.")
        .def("__repr__", [](const PositionalMatrix& p) {
            return "PositionalMatrix(size=" + std::to_string(p.size()) +
                   ", width=" + std::to_string(p.width()) + ")";
        });

    m.def("alphabet_symbols", [](const std::string& a) { return alphabet_symbols(parse_alphabet(a)); },
          nb::arg("alphabet") = "aa",
          "Symbols in code order for an alphabet; custom matrices must follow this order.");
    m.def("amino_acids", [] { return alphabet_symbols(Alphabet::AminoAcid); },
          "The amino-acid symbol order used by the built-in matrices and custom AA matrices.");

    nb::class_<PyParams>(m, "SearchParams",
                         "Search scope and budget. Scope: max_subs/max_ins/max_dels (exact, "
                         "seqtm) and max_total_edits. Budget: max_penalty with an optional "
                         "matrix (identity/BLOSUM62/PAM250/PAM100/structural) and gap costs. engine is 'auto'|'seqtrie'|'seqtm', "
                         "mode is 'all'|'top'.")
        .def("__init__",
             [](PyParams* self, int max_subs, int max_ins, int max_dels, int max_total_edits,
                long max_penalty, nb::object matrix, int gap_open, int gap_extend,
                std::string engine, std::string mode) {
                 PyParams p;
                 p.max_subs = max_subs; p.max_ins = max_ins; p.max_dels = max_dels;
                 p.max_total_edits = max_total_edits; p.max_penalty = max_penalty;
                 set_matrix(p, matrix); p.gap_open = gap_open; p.gap_extend = gap_extend;
                 parse_engine(engine); parse_mode(mode);  // validate eagerly
                 p.engine = std::move(engine); p.mode = std::move(mode);
                 new (self) PyParams(std::move(p));
             },
             nb::arg("max_subs") = 0, nb::arg("max_ins") = 0, nb::arg("max_dels") = 0,
             nb::arg("max_total_edits") = 0, nb::arg("max_penalty") = 0,
             nb::arg("matrix") = "", nb::arg("gap_open") = 1, nb::arg("gap_extend") = 1,
             nb::arg("engine") = "auto", nb::arg("mode") = "all")
        .def_rw("max_subs", &PyParams::max_subs,
                "Maximum substitutions. Enforced exactly by ``seqtm``; ``seqtrie`` ignores it.")
        .def_rw("max_ins", &PyParams::max_ins,
                "Maximum insertions (residues in the query absent from the reference). "
                "Enforced exactly by ``seqtm``; ``seqtrie`` ignores it. Indels widen the "
                "search frontier far more than substitutions do.")
        .def_rw("max_dels", &PyParams::max_dels,
                "Maximum deletions (reference residues absent from the query). Enforced "
                "exactly by ``seqtm``; ``seqtrie`` ignores it.")
        .def_rw("max_total_edits", &PyParams::max_total_edits,
                "Cap on the total edit count, independent of the per-type caps rather than "
                "clamped by them. 0 means no separate total cap -- the sum of the three "
                "per-type caps applies. This is the only edit limit ``seqtrie`` honours.")
        .def_rw("max_penalty", &PyParams::max_penalty,
                "Cap on the accumulated substitution/gap penalty. 0 means unset. Required "
                "when ``engine='seqtrie'`` is paired with a matrix: seqtrie prunes on penalty "
                "alone, so without a finite budget it walks the entire index.")
        .def_prop_rw(
            "matrix",
            [](const PyParams& p) -> nb::object {
                if (p.matrix_obj) return nb::cast(*p.matrix_obj);
                return nb::cast(p.matrix);
            },
            [](PyParams& p, const nb::object& m) { set_matrix(p, m); },
            "Substitution matrix: a :class:`SubstitutionMatrix`, or the name of a builtin "
            "(``'blosum62'``, ``'blosum45'``, ``'blosum80'``, ``'pam250'``, ``'pam100'``, "
            "``'structural'``, ``'identity'``). ``''`` (the default) means unit cost, where a "
            "score is a plain edit count. Setting a matrix makes ``gap_open`` matter: pass "
            "``2 * matrix.scale()``, or gaps come out ~14x cheaper than substitutions.")
        .def_prop_rw(
            "pos_matrix",
            [](const PyParams& p) -> nb::object {
                if (p.pos_matrix_obj) return nb::cast(*p.pos_matrix_obj);
                return nb::none();
            },
            [](PyParams& p, const nb::object& m) {
                if (m.is_none()) p.pos_matrix_obj.reset();
                else if (nb::isinstance<PositionalMatrix>(m)) p.pos_matrix_obj = nb::cast<PositionalMatrix>(m);
                else throw nb::type_error("pos_matrix must be a PositionalMatrix or None");
            },
            ":class:`PositionalMatrix` giving per-position penalties, or None. Applies only on "
            "the ``seqtm`` Hamming path and only when its width equals the query length. "
            "Setting it forces ``seqtm`` regardless of ``engine``.")
        .def_rw("gap_open", &PyParams::gap_open,
                "Cost of opening a gap. A run of length L costs "
                "``gap_open + (L - 1) * gap_extend``. The default of 1 is only right for unit "
                "cost -- with any real matrix use ``2 * matrix.scale()``.")
        .def_rw("gap_extend", &PyParams::gap_extend,
                "Cost of each additional residue in a gap run. Equal to ``gap_open`` means "
                "linear gaps; there is no separate mode flag.")
        .def_prop_rw("engine", [](const PyParams& p) { return p.engine; },
                      [](PyParams& p, std::string v) { parse_engine(v); p.engine = std::move(v); },
                      "``'seqtm'``, ``'seqtrie'``, or ``'auto'`` (the default). ``'auto'`` "
                      "always resolves to ``seqtm``: it is the only engine that enforces the "
                      "per-type caps and reports an edit breakdown. Ask for ``'seqtrie'`` by "
                      "name when a score budget is the entire specification.")
        .def_prop_rw("mode", [](const PyParams& p) { return p.mode; },
                      [](PyParams& p, std::string v) { parse_mode(v); p.mode = std::move(v); },
                      "``'all'`` (the default) returns every hit inside the scope; ``'top'`` "
                      "sorts by ``(score, ref_id)`` and keeps the best ``max_hits``.");

    nb::class_<Hit>(m, "Hit",
                    "A search result. Payload-agnostic: map ``ref_id`` back to your own "
                    "payload downstream. ``score`` is a non-negative penalty (0 == exact). "
                    "``n_subs``/``n_ins``/``n_dels`` are exact for the seqtm engine and 0 for "
                    "seqtrie. Iterable as ``(ref_id, score, n_subs, n_ins, n_dels)``.")
        .def_ro("ref_id", &Hit::ref_id,
                "Index of the matched reference in the list passed to :meth:`Index.build`.")
        .def_ro("score", &Hit::score,
                "Accumulated penalty, always >= 0, with 0 an exact match. Under unit cost this "
                "is the edit distance; with a matrix it is the summed substitution and gap "
                "penalty. Lower is better.")
        .def_ro("n_subs", &Hit::n_subs,
                "Substitutions in the best path to this reference. Exact under ``seqtm``, "
                "always 0 under ``seqtrie``, which cannot see edit types.")
        .def_ro("n_ins", &Hit::n_ins,
                "Insertions in the best path. Exact under ``seqtm``, always 0 under "
                "``seqtrie``.")
        .def_ro("n_dels", &Hit::n_dels,
                "Deletions in the best path. Exact under ``seqtm``, always 0 under "
                "``seqtrie``.")
        .def("__iter__", [](const Hit& h) {
            return nb::iter(nb::make_tuple(h.ref_id, h.score, h.n_subs, h.n_ins, h.n_dels));
        })
        .def("__repr__", [](const Hit& h) {
            return "Hit(ref_id=" + std::to_string(h.ref_id) + ", score=" + std::to_string(h.score) +
                   ", n_subs=" + std::to_string(h.n_subs) + ", n_ins=" + std::to_string(h.n_ins) +
                   ", n_dels=" + std::to_string(h.n_dels) + ")";
        });

    nb::class_<Alignment>(m, "Alignment",
                          "Global alignment of a query to a reference. ``ops`` has one char per "
                          "column: 'M' match, 'S' substitution, 'I' insertion, 'D' deletion.")
        .def_ro("aligned_query", &Alignment::aligned_query,
                "The query with ``-`` inserted at deletion columns. Same length as "
                "``aligned_ref`` and ``ops``.")
        .def_ro("aligned_ref", &Alignment::aligned_ref,
                "The reference with ``-`` inserted at insertion columns.")
        .def_ro("ops", &Alignment::ops,
                "One character per alignment column: ``M`` match, ``S`` substitution, ``I`` "
                "insertion, ``D`` deletion.")
        .def_ro("score", &Alignment::score,
                "Alignment penalty, >= 0 and 0 for an exact match -- the same convention as "
                ":attr:`Hit.score`, and the opposite of :func:`seqtree.pairwise.score`, which "
                "returns a signed similarity.")
        .def("__repr__", [](const Alignment& a) {
            return "Alignment(score=" + std::to_string(a.score) + ", ops='" + a.ops + "')";
        });

    nb::class_<Index>(m, "Index",
                      "Immutable search index over a set of reference sequences. Build once, "
                      "then query concurrently; reference id is the position in ``refs``.")
        .def_static(
            "build",
            [](std::vector<std::string> refs, const std::string& alphabet) {
                return Index::build(std::move(refs), parse_alphabet(alphabet));
            },
            nb::arg("refs"), nb::arg("alphabet") = "aa",
            "Build an index. ``alphabet`` is 'aa', 'nt', or 'iupac'. Raises ValueError "
            "on a symbol outside the alphabet.")
        .def("__len__", &Index::size)
        .def("ref_seq", [](const Index& i, uint32_t id) { return std::string(i.ref_seq(id)); },
             nb::arg("ref_id"), "Return the reference sequence string for a reference id.")
        .def("search", &py_search, nb::arg("query"), nb::arg("params"),
             "Return all hits for one query within the scope/budget in ``params``.")
        .def("search_top", &py_search_top, nb::arg("query"), nb::arg("params"), nb::arg("k") = 1,
             "Return up to ``k`` best (lowest-score) hits for one query.")
        .def("search_batch", &py_search_batch, nb::arg("queries"), nb::arg("params"),
             nb::arg("threads") = 0,
             "Search many queries in parallel (releases the GIL). ``threads=0`` uses all "
             "cores. Returns one hit list per query, in input order.")
        .def("align", &py_align, nb::arg("ref_id"), nb::arg("query"), nb::arg("params"),
             "Compute a global alignment between ``query`` and a reference, on demand.")
        .def("collisions_batch", &py_collisions_batch, nb::arg("queries"), nb::arg("params"),
             nb::arg("threads") = 0,
             "Per-query count of seqtm collisions: how often a reference was re-reached via a "
             "different edit path during branch-and-bound (0 for seqtrie / substitution-only).")
        .def("save", &Index::save, nb::arg("path"),
             "Serialize the index to a flat binary file for fast reload.")
        .def_static("load", &Index::load, nb::arg("path"),
                    "Load an index previously written with save(); raises on a corrupt/old file.");

    m.def("pairwise_batch", &py_pairwise_batch, nb::arg("a"), nb::arg("b"), nb::arg("params"),
          nb::arg("alphabet") = "aa", nb::arg("threads") = 0,
          "Batch-vs-batch search. Indexes the larger set internally and streams the smaller; "
          "results are a-major (one hit list per a[i]) with Hit.ref_id pointing into b.");

    nb::class_<ScoreMatrix>(m, "ScoreMatrix", nb::type_slots(kScoreMatrixSlots),
                            "A read-only (n_queries, n_refs) int32 penalty matrix, row-major. "
                            "Exposes the buffer protocol, so ``numpy.asarray(sm)`` and "
                            "``memoryview(sm)`` both wrap it without copying. Index it with "
                            "``sm[i, k]`` or pull one row with ``sm.row(i)``.")
        .def_prop_ro("shape",
                     [](const ScoreMatrix& s) { return nb::make_tuple(s.rows, s.cols); },
                     "``(n_rows, n_cols)`` -- queries by references, the same shape "
                     "``numpy.asarray(sm).shape`` reports.")
        .def("__len__", [](const ScoreMatrix& s) { return s.rows; },
             "Number of rows, i.e. queries.")
        .def(
            "row",
            [](const ScoreMatrix& s, size_t i) {
                if (i >= s.rows) throw nb::index_error("row out of range");
                return std::vector<int32_t>(s.data.begin() + i * s.cols,
                                            s.data.begin() + (i + 1) * s.cols);
            },
            nb::arg("i"), "Row i as a list of penalties, one per reference.")
        .def("__getitem__",
             [](const ScoreMatrix& s, std::pair<size_t, size_t> ik) {
                 if (ik.first >= s.rows || ik.second >= s.cols)
                     throw nb::index_error("index out of range");
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
        nb::arg("query"), nb::arg("ref"), nb::arg("matrix"), nb::arg("mode") = "global",
        nb::arg("gap_open") = 11, nb::arg("gap_extend") = 1, nb::arg("alphabet") = "aa",
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
        nb::arg("query"), nb::arg("ref"), nb::arg("matrix"), nb::arg("mode") = "global",
        nb::arg("gap_open") = 11, nb::arg("gap_extend") = 1, nb::arg("alphabet") = "aa",
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
                nb::gil_scoped_release release;
                out.data = align_score_matrix(q, r, mat, a, md, gap_open, gap_extend, threads);
            }
            return out;
        },
        nb::arg("queries"), nb::arg("refs"), nb::arg("matrix"), nb::arg("mode") = "global",
        nb::arg("gap_open") = 11, nb::arg("gap_extend") = 1, nb::arg("alphabet") = "aa",
        nb::arg("threads") = 0,
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
                nb::gil_scoped_release release;
                out.data = align_dist_matrix(q, r, mat, a, md, gap_open, gap_extend, threads);
            }
            return out;
        },
        nb::arg("queries"), nb::arg("refs"), nb::arg("matrix"), nb::arg("mode") = "global",
        nb::arg("gap_open") = 11, nb::arg("gap_extend") = 1, nb::arg("alphabet") = "aa",
        nb::arg("threads") = 0,
        "Dense (n_queries, n_refs) distance matrix d = s(a,a) + s(b,b) - 2*s(a,b), the "
        "sequence-level Gram transform of the alignment scores. Non-negative, zero on identity.");

    m.def("gapblock_matrix", &py_gapblock_matrix, nb::arg("queries"), nb::arg("refs"),
          nb::arg("alphabet") = "aa", nb::arg("matrix") = std::nullopt, nb::arg("gap_open") = 1,
          nb::arg("gap_extend") = 1, nb::arg("prior") = std::vector<int32_t>{},
          nb::arg("prior_width") = 0, nb::arg("threads") = 0,
          "Exhaustive single-gap-block penalties for every (query, ref) pair, GIL released. "
          "`prior` is the gap prior flattened to [m][d][i]; see seqtree.gapblock.score_matrix, "
          "which builds it for you.");

    m.def("hamming", &hamming, nb::arg("a"), nb::arg("b"),
          "Hamming distance: the number of positions at which two EQUAL-length sequences differ. "
          "Raises ValueError on a length mismatch. Case-sensitive, byte-for-byte.");

    m.def("levenshtein", &levenshtein, nb::arg("a"), nb::arg("b"),
          "Levenshtein (edit) distance: the minimum number of single-character insertions, "
          "deletions, and substitutions to turn `a` into `b`, each cost 1. Case-sensitive.");

    m.def(
        "hamming_matrix",
        [](const std::vector<std::string>& a, const std::vector<std::string>& b, int threads) {
            ScoreMatrix out;
            out.rows = a.size();
            out.cols = b.size();
            {
                nb::gil_scoped_release release;
                out.data = hamming_matrix(a, b, threads);
            }
            return out;
        },
        nb::arg("a"), nb::arg("b"), nb::arg("threads") = 0,
        "Dense (len(a), len(b)) int32 Hamming-distance matrix, GIL released. Raises ValueError "
        "if any pair has mismatched lengths.");

    m.def(
        "levenshtein_matrix",
        [](const std::vector<std::string>& a, const std::vector<std::string>& b, int threads) {
            ScoreMatrix out;
            out.rows = a.size();
            out.cols = b.size();
            {
                nb::gil_scoped_release release;
                out.data = levenshtein_matrix(a, b, threads);
            }
            return out;
        },
        nb::arg("a"), nb::arg("b"), nb::arg("threads") = 0,
        "Dense (len(a), len(b)) int32 Levenshtein-distance matrix, GIL released.");

    nb::class_<Candidate>(m, "Candidate",
                          "A seed-and-gather hit: peptide_id, shared_kmers (distinct query k-mers "
                          "that hit it), best_score. Iterable as (peptide_id, shared_kmers, best_score).")
        .def_ro("peptide_id", &Candidate::peptide_id,
                "Index of the candidate peptide in the list the index was built from.")
        .def_ro("shared_kmers", &Candidate::shared_kmers,
                "How many distinct query k-mers reached this peptide. The seed evidence: "
                "candidates are ranked on it first.")
        .def_ro("best_score", &Candidate::best_score,
                "Lowest penalty over the k-mer matches that reached this peptide. A tiebreak "
                "within a ``shared_kmers`` group, not a full-length alignment score.")
        .def("__iter__", [](const Candidate& c) {
            return nb::iter(nb::make_tuple(c.peptide_id, c.shared_kmers, c.best_score));
        })
        .def("__repr__", [](const Candidate& c) {
            return "Candidate(peptide_id=" + std::to_string(c.peptide_id) +
                   ", shared_kmers=" + std::to_string(c.shared_kmers) +
                   ", best_score=" + std::to_string(c.best_score) + ")";
        });

    nb::class_<KmerIndex>(m, "KmerIndex",
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
            nb::arg("kmers_per_peptide"), nb::arg("alphabet") = "aa",
            nb::arg("allele_ids") = std::vector<uint32_t>{},
            "Build the index. ``kmers_per_peptide[i]`` is the k-mer list for peptide ``i`` -- "
            "produce it upstream (``seqtree.layout.kmers`` masks anchors first), since which "
            "k-mers represent a peptide is a domain decision. ``allele_ids`` is an optional "
            "parallel tag per peptide, used by ``seed_and_gather``'s ``allele_filter``; pass "
            "``[]`` for none. The k-mers need not all be the same width.")
        .def("num_peptides", &KmerIndex::num_peptides,
             "Number of peptides indexed -- the length of ``kmers_per_peptide``. Also "
             "``len(index)``.")
        .def("num_kmers", &KmerIndex::num_kmers,
             "Number of distinct k-mers across all peptides. The postings lists total more "
             "than this whenever a k-mer occurs in several peptides.")
        .def("__len__", &KmerIndex::num_peptides)
        .def(
            "seed_and_gather",
            [](const KmerIndex& ki, const std::vector<std::vector<std::string>>& qk,
               const PyParams& pp, uint32_t min_shared, int64_t allele_filter, int threads) {
                auto mat = make_matrix(pp, ki.alphabet());
                SearchParams cp = to_cpp(pp, mat ? &*mat : nullptr);
                std::vector<std::vector<Candidate>> res;
                {
                    nb::gil_scoped_release release;
                    res = ki.seed_and_gather(qk, cp, min_shared, allele_filter, threads);
                }
                nb::list out;
                for (const auto& cands : res) {
                    nb::list inner;
                    for (const Candidate& c : cands) inner.append(nb::cast(c));
                    out.append(inner);
                }
                return out;
            },
            nb::arg("query_kmers"), nb::arg("params"), nb::arg("min_shared") = 1,
            nb::arg("allele_filter") = -1, nb::arg("threads") = 0,
            "For each query (its k-mer list) return ranked Candidates with >= min_shared shared "
            "k-mers; allele_filter >= 0 restricts to that allele tag.")
        .def("save", &KmerIndex::save, nb::arg("path"),
             "Write the index to ``path``. Writes a temporary and renames it into place, so a "
             "concurrent reader never sees a partial file.")
        .def_static("load", &KmerIndex::load, nb::arg("path"),
                    "Read an index written by :meth:`save`. Raises on a truncated file or a "
                    "format version this build does not know.");

    nb::class_<ArrayView>(m, "ArrayView", nb::type_slots(kArrayViewSlots),
                          "A read-only 1-D view over one of a TextResult's arrays. Exposes the "
                          "buffer protocol, so ``numpy.asarray(v)`` and ``memoryview(v)`` wrap "
                          "it without copying, and it keeps its TextResult alive.")
        .def("__len__", [](const ArrayView& a) { return a.n; })
        .def_prop_ro("format", [](const ArrayView& a) { return std::string(a.fmt); },
                     "The :mod:`struct` format code of the elements -- ``'i'`` for int32, "
                     "``'I'`` for uint32, ``'H'`` for uint16, ``'b'`` for int8. The same code "
                     "``memoryview(v).format`` reports.")
        .def("__repr__", [](const ArrayView& a) {
            return "ArrayView(len=" + std::to_string(a.n) + ", format='" + std::string(a.fmt) +
                   "')";
        });

    nb::class_<PyTextHit>(m, "TextHit",
                          "One match of a query against the text. ``offset`` is the start within "
                          "record ``ref_id``; ``n_subs`` is the Hamming distance; ``score`` is "
                          "the substitution score (0 unless a matrix was passed). "
                          "``mismatches`` lists ``(pos, query_aa, text_aa)`` -- the PAIR, so a "
                          "caller ranking by chemistry can tell L->I from L->D without "
                          "re-fetching the window. Iterable as "
                          "``(ref_id, offset, n_subs, score)``.")
        .def_ro("ref_id", &PyTextHit::ref_id,
                "Index of the record this occurrence falls in, in the list passed to "
                ":meth:`TextIndex.build`.")
        .def_ro("offset", &PyTextHit::offset,
                "Start of the occurrence within that record, 0-based. The matched substring is "
                "``ix.ref_seq(ref_id)[offset:offset + len(query)]``.")
        .def_ro("n_subs", &PyTextHit::n_subs,
                "Substituted positions, always <= the ``max_subs`` that was searched. 0 with "
                "``n_ins`` and ``n_dels`` also 0 is an exact occurrence.")
        .def_ro("n_ins", &PyTextHit::n_ins,
                "Insertions: query residues with nothing opposite them in the text. Always 0 "
                "unless the search passed ``max_indels > 0``.")
        .def_ro("n_dels", &PyTextHit::n_dels,
                "Deletions: text residues with nothing opposite them in the query. Always 0 "
                "unless the search passed ``max_indels > 0``.")
        .def_ro("length", &PyTextHit::length,
                "Residues the match spans in the text, i.e. "
                "``ref_seq(ref_id)[offset : offset + length]``. Equal to the query length "
                "whenever there are no indels; ``len(query) + n_dels - n_ins`` in general.")
        .def_ro("score", &PyTextHit::score,
                "Summed similarity over the mismatched positions under the ``matrix=`` passed "
                "to the search, or 0 when none was. Scoring never changes which hits return.")
        .def_ro("mismatches", &PyTextHit::mismatches,
                "``[(pos, query_symbol, text_symbol), ...]``, one per mismatch, ``pos`` "
                "0-based within the query. Empty for an exact hit.")
        .def("__iter__", [](const PyTextHit& h) {
            return nb::iter(nb::make_tuple(h.ref_id, h.offset, h.n_subs, h.score));
        })
        .def("__repr__", [](const PyTextHit& h) {
            return "TextHit(ref_id=" + std::to_string(h.ref_id) +
                   ", offset=" + std::to_string(h.offset) +
                   ", n_subs=" + std::to_string(h.n_subs) +
                   ", score=" + std::to_string(h.score) + ")";
        });

    nb::class_<TextResult>(m, "TextResult",
                           "Results of a TextIndex batch, held as flat parallel arrays rather "
                           "than one object per hit -- a 445k-query run comes back as a handful "
                           "of arrays. ``len(res)`` is the query count and ``res[i]`` builds the "
                           "TextHit list for query i on demand, so iterating pairs up with the "
                           "query list. For the whole batch use ``arrays()`` (zero-copy views) "
                           "or ``to_numpy()``. Hits are ordered (n_subs, ref_id, offset), stable "
                           "across runs and thread counts.")
        .def("__len__", [](const TextResult& r) { return r.num_queries(); })
        .def_prop_ro("num_hits", [](const TextResult& r) { return r.num_hits(); },
                     "Total hits across every query.")
        .def(
            "__getitem__",
            [](const TextResult& r, Py_ssize_t i) {
                const Py_ssize_t n = Py_ssize_t(r.num_queries());
                if (i < 0) i += n;
                if (i < 0 || i >= n) throw nb::index_error("query index out of range");
                std::vector<PyTextHit> out;
                for (uint32_t h = r.query_begin[size_t(i)]; h < r.query_begin[size_t(i) + 1]; ++h) {
                    PyTextHit hit{r.ref_id[h], r.offset[h],  r.n_subs[h], r.n_ins[h],
                                  r.n_dels[h], r.length[h], r.score[h],  {}};
                    for (uint32_t j = r.mm_begin[h]; j < r.mm_begin[h + 1]; ++j)
                        hit.mismatches.emplace_back(r.mm_pos[j], std::string(1, r.mm_query_aa[j]),
                                                    std::string(1, r.mm_text_aa[j]));
                    out.push_back(std::move(hit));
                }
                return out;
            },
            nb::arg("query"), "The hits for one query, as TextHit objects.")
        .def(
            "groups",
            [](const TextResult& r, Py_ssize_t i) {
                const Py_ssize_t n = Py_ssize_t(r.num_queries());
                if (i < 0) i += n;
                if (i < 0 || i >= n) throw nb::index_error("query index out of range");
                if (r.group_begin.empty())
                    throw nb::value_error("this result was not computed with group_by=True");
                std::vector<std::tuple<uint32_t, uint16_t, uint32_t>> out;
                for (uint32_t g = r.group_begin[size_t(i)]; g < r.group_begin[size_t(i) + 1]; ++g)
                    out.emplace_back(r.group_id[g], r.group_min_subs[g], r.group_n_hits[g]);
                return out;
            },
            nb::arg("query"),
            "``(group_id, min_subs, n_hits)`` per group reached by query i, sorted by group id. "
            "More than one row means the nearest parents disagree -- the tie is a first-class "
            "output rather than something each caller re-derives.")
        .def_prop_ro(
            "truncated",
            [](const TextResult& r) { return std::vector<uint8_t>(r.truncated); },
            "One flag per query: 1 if max_hits capped it. A cap that is invisible is a recall "
            "bug wearing a performance costume, so it is always reported.")
        .def("arrays", &text_result_arrays,
             "Every underlying array as a zero-copy ArrayView, keyed by name: query_begin, "
             "ref_id, offset, n_subs, score, mm_begin, mm_pos, mm_query_aa, mm_text_aa, "
             "truncated, group_begin, group_id, group_min_subs, group_n_hits.")
        .def(
            "to_numpy",
            [](const nb::object& self) {
                // Imported here, not at module scope: seqtree declares no runtime dependencies
                // and numpy stays an optional accessor.
                nb::object asarray = nb::module_::import_("numpy").attr("asarray");
                nb::dict out;
                for (auto item : text_result_arrays(self)) out[item.first] = asarray(item.second);
                return out;
            },
            "The same arrays as numpy views, sharing memory with this result. Requires numpy; "
            "it is imported on the call, never at import time.")
        .def("__repr__", [](const TextResult& r) {
            return "TextResult(queries=" + std::to_string(r.num_queries()) +
                   ", hits=" + std::to_string(r.num_hits()) + ")";
        });

    nb::class_<TextIndex>(m, "TextIndex",
                          "Exact k-mismatch (Hamming) search over a CONCATENATED reference text "
                          "-- a proteome, a genome, a transcript set. Unlike Index, which builds "
                          "a trie over reference *strings* and so needs one index per query "
                          "length, ``k`` here belongs to the index: ONE build answers every "
                          "length and every ``max_subs``. Full length, no gaps, no score in the "
                          "predicate; the answer is exact, not a heuristic.")
        .def_static(
            "build",
            [](const std::vector<std::string>& refs, const std::string& alphabet, uint8_t k,
               const std::vector<uint32_t>& group_ids) {
                return TextIndex::build(refs, parse_alphabet(alphabet), k, group_ids);
            },
            nb::arg("refs"), nb::arg("alphabet") = "aa", nb::arg("k") = 4,
            nb::arg("group_ids") = std::vector<uint32_t>{},
            "Build from whole records (NOT windows). ``k`` is the seed width; the table is "
            "direct-addressed, so alphabet_size**k buckets are allocated (24**4 = 331,776 for "
            "amino acids) and larger k is refused. ``group_ids`` optionally labels each record "
            "-- gene ids, species, clusters -- so hits can be folded onto them; seqtree does "
            "not know what a group means. Raises ValueError on a symbol outside the alphabet, "
            "naming the record.")
        .def("__len__", &TextIndex::num_refs)
        .def_prop_ro("num_refs", &TextIndex::num_refs, "Number of records in the text.")
        .def_prop_ro("num_residues", &TextIndex::num_residues,
                     "Total residues, excluding the inter-record separators.")
        .def_prop_ro("num_unknown", &TextIndex::num_unknown,
                     "Text residues outside the alphabet (36 U in the human proteome, 33 in "
                     "mouse). They are kept as holes -- no hit may cross one -- and reported "
                     "here rather than dropped silently. A *query* containing one is refused.")
        .def_prop_ro("k", &TextIndex::k, "Seed width this index was built with.")
        .def_prop_ro("alphabet",
                     [](const TextIndex& t) { return alphabet_name(t.alphabet()); },
                     "'aa', 'nt' or 'iupac'.")
        .def("ref_seq", &TextIndex::ref_seq, nb::arg("ref_id"),
             "Return record ``ref_id`` as a string.")
        .def(
            "search_batch",
            [](const TextIndex& ix, const std::vector<std::string>& queries, uint16_t max_subs,
               uint16_t max_indels, bool exclude_exact, bool best_only, bool group_by,
               uint32_t max_hits, const std::optional<SubstitutionMatrix>& matrix, int threads) {
                TextQueryOpts o;
                o.max_subs = max_subs;
                o.max_indels = max_indels;
                o.exclude_exact = exclude_exact;
                o.best_only = best_only;
                o.group_by = group_by;
                o.max_hits = max_hits;
                o.matrix = matrix ? &*matrix : nullptr;
                TextResult res;
                {
                    nb::gil_scoped_release release;  // pure C++, no Python objects touched
                    res = ix.search_batch(queries, o, threads);
                }
                return res;
            },
            nb::arg("queries"), nb::arg("max_subs") = 0, nb::arg("max_indels") = 0,
            nb::arg("exclude_exact") = false, nb::arg("best_only") = false,
            nb::arg("group_by") = false, nb::arg("max_hits") = 0,
            nb::arg("matrix") = std::nullopt, nb::arg("threads") = 0,
            "Find every position matching each query within ``max_subs`` substitutions and "
            "``max_indels`` insertions/deletions (releases the GIL; ``threads=0`` uses all "
            "cores). The two caps are independent: ``max_subs=2, max_indels=1`` accepts two "
            "substitutions AND one gap, not three edits of any kind.\n\n"
            "``max_indels=0`` is the pure Hamming predicate: every hit is exactly "
            "``len(query)`` residues wide and carries per-column ``mismatches``. Above 0 the "
            "match length varies, so read :attr:`TextHit.length` rather than assuming the query "
            "length, and ``mismatches`` is empty -- the counts are reported, the individual "
            "substituted columns are not. Indel search also needs every seed block to be exact, "
            "so a query must be at least ``(max_subs + max_indels + 1) * k`` long; a shorter one "
            "raises rather than being answered incompletely.\n\n"
            "``best_only`` walks the distance upward and stops at the first shell with any hit, "
            "returning ALL of it. ``exclude_exact`` drops 0-mismatch hits. ``max_hits`` caps a "
            "query after sorting, so the best hits survive, and sets ``truncated``. ``matrix`` "
            "only SCORES hits the predicate already accepted -- it never changes which are "
            "returned, and it is ignored when ``max_indels > 0``.")
        .def("save", &TextIndex::save, nb::arg("path"),
             "Write a flat, mmap-able index file.")
        .def_static("load", &TextIndex::load, nb::arg("path"), nb::arg("mmap") = true,
                    "Load an index written by save(). With ``mmap`` the file is mapped rather "
                    "than read, so several processes share one copy of the pages.")
        .def("__repr__", [](const TextIndex& t) {
            return "TextIndex(refs=" + std::to_string(t.num_refs()) +
                   ", residues=" + std::to_string(t.num_residues()) +
                   ", k=" + std::to_string(t.k()) + ")";
        });
}
