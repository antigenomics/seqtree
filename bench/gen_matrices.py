#!/usr/bin/env python3
"""Emit src/{pam250,pam100,structural}.inc — verified, validated substitution matrices.

* PAM250 / PAM100: EMBOSS ``EPAM250`` / ``EPAM100`` (standard NCBI log-odds), in the
  same 24-symbol order as BLOSUM62 (A R N D C Q E G H I L K M F P S T W Y V B Z X *).
* structural: a Miyazawa–Jernigan **interaction-strength** similarity. Each residue's
  interaction strength ``q(a) = mean_b e(a,b)`` is derived from the MJ residue–residue
  contact potential (Miyazawa & Jernigan, J Mol Biol 1996; the "MJ" potential tabulated
  by Keskin et al. 1998), then ``sim(a,b) = 10·(1 − |q̂(a) − q̂(b)|)`` on min-max-
  normalised strengths. This orders the 20 residues from strong interactors (hydrophobic
  F W C L Y M I V) to weak ones (polar/charged S Q D E K) — the strong/weak-interactor
  axis used in TCR-recognition models (Košmrlj et al., PNAS 2008,
  doi:10.1073/pnas.0808081105). Ambiguity codes B/Z/X are the standard combos
  (B = mean(N,D), Z = mean(Q,E), X = mean of all 20 in embedding space) and ``*`` is a
  maximally dissimilar self-matching stop.

All grids are validated square + symmetric + against documented anchors before writing.
Run once; the .inc files are committed.

  python bench/gen_matrices.py
"""
from pathlib import Path

ORDER = "A R N D C Q E G H I L K M F P S T W Y V B Z X *".split()
N = len(ORDER)
IDX = {a: i for i, a in enumerate(ORDER)}

# --- EMBOSS PAM grids (row order == ORDER) -----------------------------------
PAM250 = {
    "A": [2, -2, 0, 0, -2, 0, 0, 1, -1, -1, -2, -1, -1, -3, 1, 1, 1, -6, -3, 0, 0, 0, 0, -8],
    "R": [-2, 6, 0, -1, -4, 1, -1, -3, 2, -2, -3, 3, 0, -4, 0, 0, -1, 2, -4, -2, -1, 0, -1, -8],
    "N": [0, 0, 2, 2, -4, 1, 1, 0, 2, -2, -3, 1, -2, -3, 0, 1, 0, -4, -2, -2, 2, 1, 0, -8],
    "D": [0, -1, 2, 4, -5, 2, 3, 1, 1, -2, -4, 0, -3, -6, -1, 0, 0, -7, -4, -2, 3, 3, -1, -8],
    "C": [-2, -4, -4, -5, 12, -5, -5, -3, -3, -2, -6, -5, -5, -4, -3, 0, -2, -8, 0, -2, -4, -5, -3, -8],
    "Q": [0, 1, 1, 2, -5, 4, 2, -1, 3, -2, -2, 1, -1, -5, 0, -1, -1, -5, -4, -2, 1, 3, -1, -8],
    "E": [0, -1, 1, 3, -5, 2, 4, 0, 1, -2, -3, 0, -2, -5, -1, 0, 0, -7, -4, -2, 3, 3, -1, -8],
    "G": [1, -3, 0, 1, -3, -1, 0, 5, -2, -3, -4, -2, -3, -5, 0, 1, 0, -7, -5, -1, 0, 0, -1, -8],
    "H": [-1, 2, 2, 1, -3, 3, 1, -2, 6, -2, -2, 0, -2, -2, 0, -1, -1, -3, 0, -2, 1, 2, -1, -8],
    "I": [-1, -2, -2, -2, -2, -2, -2, -3, -2, 5, 2, -2, 2, 1, -2, -1, 0, -5, -1, 4, -2, -2, -1, -8],
    "L": [-2, -3, -3, -4, -6, -2, -3, -4, -2, 2, 6, -3, 4, 2, -3, -3, -2, -2, -1, 2, -3, -3, -1, -8],
    "K": [-1, 3, 1, 0, -5, 1, 0, -2, 0, -2, -3, 5, 0, -5, -1, 0, 0, -3, -4, -2, 1, 0, -1, -8],
    "M": [-1, 0, -2, -3, -5, -1, -2, -3, -2, 2, 4, 0, 6, 0, -2, -2, -1, -4, -2, 2, -2, -2, -1, -8],
    "F": [-3, -4, -3, -6, -4, -5, -5, -5, -2, 1, 2, -5, 0, 9, -5, -3, -3, 0, 7, -1, -4, -5, -2, -8],
    "P": [1, 0, 0, -1, -3, 0, -1, 0, 0, -2, -3, -1, -2, -5, 6, 1, 0, -6, -5, -1, -1, 0, -1, -8],
    "S": [1, 0, 1, 0, 0, -1, 0, 1, -1, -1, -3, 0, -2, -3, 1, 2, 1, -2, -3, -1, 0, 0, 0, -8],
    "T": [1, -1, 0, 0, -2, -1, 0, 0, -1, 0, -2, 0, -1, -3, 0, 1, 3, -5, -3, 0, 0, -1, 0, -8],
    "W": [-6, 2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4, 0, -6, -2, -5, 17, 0, -6, -5, -6, -4, -8],
    "Y": [-3, -4, -2, -4, 0, -4, -4, -5, 0, -1, -1, -4, -2, 7, -5, -3, -3, 0, 10, -2, -3, -4, -2, -8],
    "V": [0, -2, -2, -2, -2, -2, -2, -1, -2, 4, 2, -2, 2, -1, -1, -1, 0, -6, -2, 4, -2, -2, -1, -8],
    "B": [0, -1, 2, 3, -4, 1, 3, 0, 1, -2, -3, 1, -2, -4, -1, 0, 0, -5, -3, -2, 3, 2, -1, -8],
    "Z": [0, 0, 1, 3, -5, 3, 3, 0, 2, -2, -3, 0, -2, -5, 0, 0, -1, -6, -4, -2, 2, 3, -1, -8],
    "X": [0, -1, 0, -1, -3, -1, -1, -1, -1, -1, -1, -1, -1, -2, -1, 0, 0, -4, -2, -1, -1, -1, -1, -8],
    "*": [-8] * 23 + [1],
}
PAM100 = {
    "A": [4, -3, -1, -1, -3, -2, 0, 1, -3, -2, -3, -3, -2, -5, 1, 1, 1, -7, -4, 0, -1, -1, -1, -9],
    "R": [-3, 7, -2, -4, -5, 1, -3, -5, 1, -3, -5, 2, -1, -6, -1, -1, -3, 1, -6, -4, -3, -1, -2, -9],
    "N": [-1, -2, 5, 3, -5, -1, 1, -1, 2, -3, -4, 1, -4, -5, -2, 1, 0, -5, -2, -3, 4, 0, -1, -9],
    "D": [-1, -4, 3, 5, -7, 0, 4, -1, -1, -4, -6, -1, -5, -8, -3, -1, -2, -9, -6, -4, 4, 3, -2, -9],
    "C": [-3, -5, -5, -7, 9, -8, -8, -5, -4, -3, -8, -8, -7, -7, -4, -1, -4, -9, -1, -3, -6, -8, -5, -9],
    "Q": [-2, 1, -1, 0, -8, 6, 2, -3, 3, -4, -2, 0, -2, -7, -1, -2, -2, -7, -6, -3, 0, 5, -2, -9],
    "E": [0, -3, 1, 4, -8, 2, 5, -1, -1, -3, -5, -1, -4, -8, -2, -1, -2, -9, -5, -3, 3, 4, -2, -9],
    "G": [1, -5, -1, -1, -5, -3, -1, 5, -4, -5, -6, -3, -4, -6, -2, 0, -2, -9, -7, -3, -1, -2, -2, -9],
    "H": [-3, 1, 2, -1, -4, 3, -1, -4, 7, -4, -3, -2, -4, -3, -1, -2, -3, -4, -1, -3, 1, 1, -2, -9],
    "I": [-2, -3, -3, -4, -3, -4, -3, -5, -4, 6, 1, -3, 1, 0, -4, -3, 0, -7, -3, 3, -3, -3, -2, -9],
    "L": [-3, -5, -4, -6, -8, -2, -5, -6, -3, 1, 6, -4, 3, 0, -4, -4, -3, -3, -3, 0, -5, -4, -3, -9],
    "K": [-3, 2, 1, -1, -8, 0, -1, -3, -2, -3, -4, 5, 0, -7, -3, -1, -1, -6, -6, -4, 0, -1, -2, -9],
    "M": [-2, -1, -4, -5, -7, -2, -4, -4, -4, 1, 3, 0, 9, -1, -4, -3, -1, -6, -5, 1, -4, -2, -2, -9],
    "F": [-5, -6, -5, -8, -7, -7, -8, -6, -3, 0, 0, -7, -1, 8, -6, -4, -5, -1, 4, -3, -6, -7, -4, -9],
    "P": [1, -1, -2, -3, -4, -1, -2, -2, -1, -4, -4, -3, -4, -6, 7, 0, -1, -7, -7, -3, -3, -1, -2, -9],
    "S": [1, -1, 1, -1, -1, -2, -1, 0, -2, -3, -4, -1, -3, -4, 0, 4, 2, -3, -4, -2, 0, -2, -1, -9],
    "T": [1, -3, 0, -2, -4, -2, -2, -2, -3, 0, -3, -1, -1, -5, -1, 2, 5, -7, -4, 0, -1, -2, -1, -9],
    "W": [-7, 1, -5, -9, -9, -7, -9, -9, -4, -7, -3, -6, -6, -1, -7, -3, -7, 12, -2, -9, -6, -8, -6, -9],
    "Y": [-4, -6, -2, -6, -1, -6, -5, -7, -1, -3, -3, -6, -5, 4, -7, -4, -4, -2, 9, -4, -4, -6, -4, -9],
    "V": [0, -4, -3, -4, -3, -3, -3, -3, -3, 3, 0, -4, 1, -3, -3, -2, 0, -9, -4, 5, -4, -3, -2, -9],
    "B": [-1, -3, 4, 4, -6, 0, 3, -1, 1, -3, -5, 0, -4, -6, -3, 0, -1, -6, -4, -4, 4, 2, -2, -9],
    "Z": [-1, -1, 0, 3, -8, 5, 4, -2, 1, -3, -4, -1, -2, -7, -1, -2, -2, -8, -6, -3, 2, 5, -2, -9],
    "X": [-1, -2, -1, -2, -5, -2, -2, -2, -2, -2, -3, -2, -2, -4, -2, -1, -1, -6, -4, -2, -2, -2, -2, -9],
    "*": [-9] * 23 + [1],
}


# --- Miyazawa–Jernigan contact energies e(a,b) (row order == MJ_AA) ----------
# The "MJ" residue–residue contact potential (Miyazawa & Jernigan, J Mol Biol 1996),
# as tabulated by Keskin et al. (Protein Sci 1998) and redistributed in the tcren
# package (data/MJ_Keskin_potentials.csv). More negative = a stronger contact. The
# source is perfectly symmetric, so the e(a,b) = (e[a][b] + e[b][a]) / 2 in
# structural_grid() is a no-op guard, not a repair.
#
# A-N / N-A is transcribed as 0.15 but you will not find an "A,N" row in the CSV: that
# pair is corrupted upstream. The CSV emits the lower triangle in the order A R N D C Q
# E G H I L K M F P S T W Y V, and the 4th slot -- which must be N-A -- reads "V,1"
# (mirrored at "1,V"), where "1" is a mangled residue symbol. 0.15 is that row's value.
# It is not a duplicate V-N: the V row separately and correctly lists N = 0.12.
MJ_AA = list("ACDEFGHIKLMNPQRSTVWY")
MJ_CONTACT = {
    "A": [-0.12, -0.33, 0.27, 0.38, -0.36, -0.08, 0.07, -0.37, 0.41, -0.38, -0.27, 0.15, 0.15, 0.22, 0.24, 0.10, 0.04, -0.32, -0.27, -0.20],
    "C": [-0.33, -1.19, 0.12, 0.20, -0.67, -0.31, -0.36, -0.64, 0.33, -0.65, -0.61, -0.01, -0.18, -0.07, 0.08, -0.13, -0.15, -0.59, -0.66, -0.39],
    "D": [0.27, 0.12, 0.29, 0.44, 0.18, 0.11, -0.10, 0.22, -0.01, 0.27, 0.30, 0.02, 0.33, 0.24, -0.24, 0.10, 0.11, 0.36, 0.07, -0.07],
    "E": [0.38, 0.20, 0.44, 0.46, 0.14, 0.32, 0.00, 0.17, -0.06, 0.17, 0.12, 0.12, 0.37, 0.27, -0.22, 0.18, 0.16, 0.26, 0.00, -0.08],
    "F": [-0.36, -0.67, 0.18, 0.14, -0.88, -0.19, -0.34, -0.73, 0.19, -0.80, -0.83, -0.01, -0.19, -0.11, -0.05, -0.12, -0.15, -0.67, -0.68, -0.58],
    "G": [-0.08, -0.31, 0.11, 0.32, -0.19, -0.29, 0.00, -0.13, 0.29, -0.16, -0.17, -0.01, 0.02, 0.13, 0.09, -0.01, -0.04, -0.15, -0.25, -0.22],
    "H": [0.07, -0.36, -0.10, 0.00, -0.34, 0.00, -0.40, -0.13, 0.38, -0.18, -0.29, 0.00, 0.01, 0.15, 0.05, 0.04, -0.03, -0.06, -0.37, -0.30],
    "I": [-0.37, -0.64, 0.22, 0.17, -0.73, -0.13, -0.13, -0.74, 0.24, -0.81, -0.66, 0.14, -0.05, -0.01, 0.00, 0.03, -0.15, -0.67, -0.60, -0.49],
    "K": [0.41, 0.33, -0.01, -0.06, 0.19, 0.29, 0.38, 0.24, 0.76, 0.22, 0.29, 0.22, 0.47, 0.28, 0.66, 0.36, 0.33, 0.29, 0.09, -0.05],
    "L": [-0.38, -0.65, 0.27, 0.17, -0.80, -0.16, -0.18, -0.81, 0.22, -0.84, -0.70, 0.04, -0.12, -0.04, -0.04, -0.02, -0.15, -0.74, -0.62, -0.55],
    "M": [-0.27, -0.61, 0.30, 0.12, -0.83, -0.17, -0.29, -0.66, 0.29, -0.70, -0.70, 0.04, -0.13, -0.06, 0.03, 0.05, -0.11, -0.51, -0.73, -0.56],
    "N": [0.15, -0.01, 0.02, 0.12, -0.01, -0.01, 0.00, 0.14, 0.22, 0.04, 0.04, -0.06, 0.18, 0.06, 0.10, 0.09, 0.04, 0.12, -0.10, -0.11],
    "P": [0.15, -0.18, 0.33, 0.37, -0.19, 0.02, 0.01, -0.05, 0.47, -0.12, -0.13, 0.18, 0.11, 0.17, 0.17, 0.20, 0.13, -0.05, -0.37, -0.25],
    "Q": [0.22, -0.07, 0.24, 0.27, -0.11, 0.13, 0.15, -0.01, 0.28, -0.04, -0.06, 0.06, 0.17, 0.20, 0.09, 0.22, 0.12, 0.08, -0.02, -0.14],
    "R": [0.24, 0.08, -0.24, -0.22, -0.05, 0.09, 0.05, 0.00, 0.66, -0.04, 0.03, 0.10, 0.17, 0.09, 0.19, 0.16, 0.11, 0.08, -0.21, -0.25],
    "S": [0.10, -0.13, 0.10, 0.18, -0.12, -0.01, 0.04, 0.03, 0.36, -0.02, 0.05, 0.09, 0.20, 0.22, 0.16, 0.05, 0.04, 0.04, -0.01, -0.08],
    "T": [0.04, -0.15, 0.11, 0.16, -0.15, -0.04, -0.03, -0.15, 0.33, -0.15, -0.11, 0.04, 0.13, 0.12, 0.11, 0.04, 0.03, -0.07, -0.02, -0.09],
    "V": [-0.32, -0.59, 0.36, 0.26, -0.67, -0.15, -0.06, -0.67, 0.29, -0.74, -0.51, 0.12, -0.05, 0.08, 0.08, 0.04, -0.07, -0.65, -0.51, -0.38],
    "W": [-0.27, -0.66, 0.07, 0.00, -0.68, -0.25, -0.37, -0.60, 0.09, -0.62, -0.73, -0.10, -0.37, -0.02, -0.21, -0.01, -0.02, -0.51, -0.64, -0.49],
    "Y": [-0.20, -0.39, -0.07, -0.08, -0.58, -0.22, -0.30, -0.49, -0.05, -0.55, -0.56, -0.11, -0.25, -0.14, -0.25, -0.08, -0.09, -0.38, -0.49, -0.45],
}


def grid_from_rows(rows):
    return [[rows[a][j] for j in range(N)] for a in ORDER]


def validate(M, hi, lo, w_diag):
    assert all(len(r) == N for r in M), "non-square"
    for i in range(N):
        for j in range(N):
            assert M[i][j] == M[j][i], f"asymmetric at {ORDER[i]},{ORDER[j]}"
    flat = [v for r in M for v in r]
    assert max(flat) == hi and min(flat) == lo, f"anchors: max={max(flat)} min={min(flat)}"
    assert M[IDX["W"]][IDX["W"]] == w_diag, "W-W diagonal mismatch"


def structural_grid():
    """MJ interaction-strength similarity on 0..10: q(a)=mean_b e(a,b) on the symmetrised
    MJ contact potential, min-max normalised; sim(a,b)=10·(1−|q̂(a)−q̂(b)|). Equal-strength
    residues score 10; strong (hydrophobic) and weak (polar/charged) interactors separate."""
    aa20 = MJ_AA
    j = {a: i for i, a in enumerate(aa20)}
    e = lambda a, b: (MJ_CONTACT[a][j[b]] + MJ_CONTACT[b][j[a]]) / 2     # symmetrise the source
    q = {a: sum(e(a, b) for b in aa20) / len(aa20) for a in aa20}        # per-residue interaction strength
    lo, hi = min(q.values()), max(q.values())
    qn = {a: (q[a] - lo) / (hi - lo) for a in aa20}
    sim = {(a, b): 10 * (1 - abs(qn[a] - qn[b])) for a in aa20 for b in aa20}

    # Virtual residues as weighted combos in the (Gram) embedding: bilinear extension.
    combos = {"B": {"N": 0.5, "D": 0.5}, "Z": {"Q": 0.5, "E": 0.5},
              "X": {a: 1 / 20 for a in aa20}}

    def s(a, b):  # similarity for any pair of base/virtual symbols (not '*')
        wa = combos.get(a, {a: 1.0})
        wb = combos.get(b, {b: 1.0})
        return sum(wi * wj * sim[(ai, bj)] for ai, wi in wa.items() for bj, wj in wb.items())

    M = [[0] * N for _ in range(N)]
    syms = [a for a in ORDER if a != "*"]
    for a in syms:
        for b in syms:
            M[IDX[a]][IDX[b]] = round(s(a, b))
    # stop '*': maximally dissimilar (0 on the 0..10 scale) to all, identical to itself.
    star = IDX["*"]
    for i in range(N):
        M[star][i] = M[i][star] = 0
    M[star][star] = 10
    return M


def emit(name, var, M, header):
    lines = [f"// {header}",
             f"// Row/column order: {' '.join(ORDER)}",
             "// Generated by bench/gen_matrices.py -- do not edit by hand.",
             f"static const int32_t {var}Size = {N};",
             f"static const int32_t {var}[{N} * {N}] = {{",
             "//   " + "   ".join(f"{a:>2}" for a in ORDER)]
    for i, a in enumerate(ORDER):
        lines.append("    " + ", ".join(f"{v:>3}" for v in M[i]) + ",")
    lines.append("};")
    Path(f"src/{name}.inc").write_text("\n".join(lines) + "\n")
    print(f"wrote src/{name}.inc")


def main():
    pam250 = grid_from_rows(PAM250)
    validate(pam250, hi=17, lo=-8, w_diag=17)
    emit("pam250", "kPam250", pam250, "PAM250 (EMBOSS EPAM250, NCBI log-odds; similarity).")

    pam100 = grid_from_rows(PAM100)
    validate(pam100, hi=12, lo=-9, w_diag=12)
    emit("pam100", "kPam100", pam100, "PAM100 (EMBOSS EPAM100, NCBI log-odds; similarity).")

    structural = structural_grid()
    # Anchors: the 20 standard AAs (and '*') self-match at 10; B/Z/X are virtual combos
    # so their diagonal is < 10 by construction. Whole grid stays symmetric and in 0..10.
    standard = [a for a in ORDER if a not in ("B", "Z", "X")]
    assert all(structural[IDX[a]][IDX[a]] == 10 for a in standard), "structural diagonal != 10"
    assert 0 <= min(v for r in structural for v in r) and max(v for r in structural for v in r) == 10
    for i in range(N):
        for j in range(N):
            assert structural[i][j] == structural[j][i], "structural not symmetric"
    emit("structural", "kStructural", structural,
         "Miyazawa-Jernigan interaction-strength similarity (strong/weak interactors; 0..10).")


if __name__ == "__main__":
    main()
