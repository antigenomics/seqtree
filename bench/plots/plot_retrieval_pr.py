#!/usr/bin/env python3
"""Plot the retrieval precision/recall table — consumes a TSV, produces an SVG.

This is a *plot step*: it reads the table written by
``bench/tables/gen_retrieval_table.py`` and renders an SVG with gnuplot. It does
not run seqtree, so plotting never re-measures anything (the numbers come only
from the committed table).

  python bench/plots/plot_retrieval_pr.py            # tables/retrieval_pr.tsv -> figures/retrieval_pr.svg
"""
import argparse
import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def read_table(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or line.startswith("penalty"):
            continue
        k, tp, fp, fn, prec, rec, f1 = line.split("\t")
        rows.append((int(k), float(prec), float(rec), float(f1)))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default=str(REPO / "bench/tables/retrieval_pr.tsv"))
    ap.add_argument("--out", default=str(REPO / "bench/figures/retrieval_pr.svg"))
    args = ap.parse_args()

    rows = read_table(Path(args.table))
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    work = out.parent

    # Panel 1 data: precision vs recall (PR curve). Panel 2: P/R/F1 vs penalty.
    (work / "_pr_curve.tsv").write_text(
        "recall\tprecision\n" + "\n".join(f"{r:.6f}\t{p:.6f}" for _, p, r, _ in rows) + "\n")
    (work / "_pr_thresh.tsv").write_text(
        "penalty\tprecision\trecall\tf1\n"
        + "\n".join(f"{k}\t{p:.6f}\t{r:.6f}\t{f:.6f}" for k, p, r, f in rows) + "\n")

    f1max = max(rows, key=lambda x: x[3])
    gp = f"""set terminal svg size 720,760 font 'Helvetica,12' background rgb 'white'
set output '{out.name}'
set datafile separator "\\t"
set grid
set multiplot layout 2,1
set title 'BLOSUM62-weighted retrieval of Hamming-near neighbours'
set xlabel 'recall'
set ylabel 'precision'
set xrange [0:1.02]
set yrange [0.8:1.02]
set key bottom left
plot '_pr_curve.tsv' using 1:2 with linespoints lw 2 pt 7 ps 0.5 lc rgb '#1f77b4' title 'precision-recall'
set title 'precision / recall / F1 vs penalty budget'
set xlabel 'max BLOSUM62 penalty'
set ylabel 'score'
set xrange [*:*]
set yrange [0:1.05]
set key outside right top
set label 1 'best F1 = {f1max[3]:.3f} @ penalty {f1max[0]}' at graph 0.05,0.15 front
plot '_pr_thresh.tsv' using 1:2 with lines lw 2 dt 2 lc rgb '#d62728' title 'precision', \\
     '' using 1:3 with lines lw 2 dt 3 lc rgb '#1f77b4' title 'recall', \\
     '' using 1:4 with lines lw 2 lc rgb '#2ca02c' title 'F1'
unset multiplot
"""
    (work / "retrieval_pr.gp").write_text(gp)
    subprocess.run(["gnuplot", "retrieval_pr.gp"], cwd=work, check=True)
    (work / "_pr_curve.tsv").unlink()
    (work / "_pr_thresh.tsv").unlink()
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
