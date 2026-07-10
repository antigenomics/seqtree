#!/usr/bin/env python3
"""Generate the bundled control subset shipped with the wheel.

Streams the human TRB amino-acid control from isalgo/airr_control and writes unique productive
CDR3aa, one per line, gzip-compressed.

The upstream table is sorted by clonotype abundance. Taking the first N rows -- which this script
used to do -- returns the N most expanded public clones, which is 25.8x more self-similar than a
uniform sample of the same size and inflates mean control ball mass 3.1x at a BLOSUM62 budget of
40. That violates ass:indep in appendix/evalue.tex ("the unique clonotypes of C are i.i.d. from
P0") and quietly deflates every E-value. Reservoir-sample instead: one full pass, every unique
clonotype equally likely.

    python bench/gen_control.py --n 250000 --seed 0
"""
import argparse
import gzip
import io
import random
import urllib.request

URL = "https://huggingface.co/datasets/isalgo/airr_control/resolve/main/human.trb.aa.vdjtools.tsv.gz"
OUT = "python/seqtree/data/control_human_trb_aa.txt.gz"
AA = set("ACDEFGHIKLMNPQRSTVWY")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=250_000)
    ap.add_argument("--url", default=URL)
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    seen = set()
    order = []
    n_seen = 0
    req = urllib.request.Request(args.url, headers={"User-Agent": "seqtree-gen-control"})
    with urllib.request.urlopen(req) as resp:
        with gzip.GzipFile(fileobj=resp) as gz:
            text = io.TextIOWrapper(gz, encoding="utf-8")
            header = text.readline().rstrip("\n").split("\t")
            cols = {c.lower(): i for i, c in enumerate(header)}
            ci = cols.get("cdr3aa", cols.get("cdr3.aa"))
            if ci is None:
                raise SystemExit(f"no cdr3aa column in header: {header}")
            for line in text:
                f = line.rstrip("\n").split("\t")
                if len(f) <= ci:
                    continue
                s = f[ci].strip().upper()
                if not s or s in seen or not all(c in AA for c in s):
                    continue          # drop '_' (out of frame), '*' (stop), ambiguous residues
                seen.add(s)
                n_seen += 1           # reservoir over UNIQUE clonotypes, not over rows
                if len(order) < args.n:
                    order.append(s)
                else:
                    j = rng.randrange(n_seen)
                    if j < args.n:
                        order[j] = s

    order.sort()                      # order carries no information; make the asset deterministic
    with gzip.open(args.out, "wt", encoding="utf-8") as out:
        out.write("\n".join(order) + "\n")
    print(f"wrote {len(order):,} unique productive CDR3aa (uniform over {n_seen:,}) to {args.out}")


if __name__ == "__main__":
    main()
