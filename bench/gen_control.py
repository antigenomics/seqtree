#!/usr/bin/env python3
"""Generate the bundled control subset shipped with the wheel.

Streams the human TRB amino-acid control from isalgo/airr_control (reading only
enough of the 424 MB gzip to collect the subset) and writes unique CDR3aa, one
per line, gzip-compressed. Re-run to refresh the asset.

    python bench/gen_control.py --n 250000
"""
import argparse
import gzip
import io
import urllib.request

URL = "https://huggingface.co/datasets/isalgo/airr_control/resolve/main/human.trb.aa.vdjtools.tsv.gz"
OUT = "python/seqtree/data/control_human_trb_aa.txt.gz"
AA = set("ACDEFGHIKLMNPQRSTVWY")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=250_000)
    ap.add_argument("--url", default=URL)
    ap.add_argument("--out", default=OUT)
    args = ap.parse_args()

    seen = set()
    order = []
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
                if s and s not in seen and all(c in AA for c in s):
                    seen.add(s)
                    order.append(s)
                    if len(order) >= args.n:
                        break

    with gzip.open(args.out, "wt", encoding="utf-8") as out:
        out.write("\n".join(order) + "\n")
    print(f"wrote {len(order)} unique CDR3aa to {args.out}")


if __name__ == "__main__":
    main()
