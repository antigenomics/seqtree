#!/usr/bin/env python3
"""Generate N OLGA human TRB CDR3 amino-acid sequences to a gzipped cache file.

Used to build OLGA-generation null backgrounds for the E-value benchmark.

    python bench/gen_olga.py --n 1000000 --out bench/cache/olga_1M.txt.gz
"""
import argparse
import gzip
import os
import subprocess
import sys
import tempfile

AA = set("ACDEFGHIKLMNPQRSTVWY")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    cli = os.path.join(os.path.dirname(sys.executable), "olga-generate_sequences")
    if not os.path.exists(cli):
        cli = "olga-generate_sequences"
    tmp = os.path.join(tempfile.mkdtemp(), "olga.tsv")
    # over-generate slightly so the unique count meets the target after dedup
    subprocess.run([cli, "--humanTRB", "-n", str(int(args.n * 1.05)), "-o", tmp], check=True)

    seen = set()
    out = []
    with open(tmp) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                s = f[1].strip().upper()
                if s and s not in seen and all(c in AA for c in s):
                    seen.add(s)
                    out.append(s)
                    if len(out) >= args.n:
                        break
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with gzip.open(args.out, "wt", encoding="utf-8") as o:
        o.write("\n".join(out) + "\n")
    os.remove(tmp)
    print(f"wrote {len(out)} unique CDR3aa to {args.out}")


if __name__ == "__main__":
    main()
