#!/usr/bin/env python3
"""Emit src/{pam250,pam100,structural}.inc — verified, validated substitution matrices.

* PAM250 / PAM100: EMBOSS ``EPAM250`` / ``EPAM100`` (standard NCBI log-odds), in the
  same 24-symbol order as BLOSUM62 (A R N D C Q E G H I L K M F P S T W Y V B Z X *).
* structural: the TeXshade physicochemical-similarity matrix (sidechain volume +
  hydropathy; CTAN texshade.sty), parsed from the installed package. TeXshade defines
  the 20 standard amino acids on a 0..10 similarity scale; the ambiguity codes B/Z/X
  and stop ``*`` are added by the standard convention (B = mean(N,D), Z = mean(Q,E),
  X = mean of all 20 in embedding space, ``*`` = a maximally dissimilar self-matching
  stop). One asymmetric entry in the TeXshade source (S/M: 3 vs 5) is symmetrised to 4.

All grids are validated square + symmetric + against documented anchors before writing.
Run once; the .inc files are committed.

  python bench/gen_matrices.py
"""
import re
import subprocess
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
    sty = subprocess.run(["kpsewhich", "texshade.sty"], capture_output=True, text=True).stdout.strip()
    text = Path(sty).read_text()
    block = text[text.index("%%%% structural"):text.index("%%%% PAM250")]
    sim = {(a, b): int(v) for a, b, v in re.findall(r"\\cons([A-Z])([A-Z])\{(-?\d+)\}", block)}
    aa20 = sorted({a for a, _ in sim})
    assert len(aa20) == 20, f"expected 20 residues, got {len(aa20)}"
    # Symmetrise the lone asymmetric source entry (S/M) by averaging.
    for a in aa20:
        for b in aa20:
            if sim[(a, b)] != sim[(b, a)]:
                m = round((sim[(a, b)] + sim[(b, a)]) / 2)
                sim[(a, b)] = sim[(b, a)] = m

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
         "TeXshade structural similarity (sidechain volume + hydropathy; 0..10).")


if __name__ == "__main__":
    main()
