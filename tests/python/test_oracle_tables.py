"""Oracle regression: regenerating the benchmark table from the compiled seqtree in
this repo must reproduce the committed table exactly. A diff means either a real
behaviour change (update the oracle: ``bench/tables/make_tables.sh``) or a bug."""
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PRODUCER = REPO / "bench" / "tables" / "gen_retrieval_table.py"
ORACLE = REPO / "bench" / "tables" / "retrieval_pr.tsv"


def _parse(text):
    header, rows = None, []
    for line in text.splitlines():
        if line.startswith("#"):
            header = line
        elif line.startswith("penalty"):
            continue
        elif line.strip():
            f = line.split("\t")
            rows.append((int(f[0]), int(f[1]), int(f[2]), int(f[3]),
                         float(f[4]), float(f[5]), float(f[6])))
    return header, rows


@pytest.mark.skipif(not PRODUCER.exists() or not ORACLE.exists(), reason="bench tables absent")
def test_retrieval_table_matches_oracle():
    fresh = subprocess.run([sys.executable, str(PRODUCER)], cwd=REPO,
                           capture_output=True, text=True, check=True).stdout
    fh, fr = _parse(fresh)
    oh, orr = _parse(ORACLE.read_text())

    assert fh == oh, f"table parameters changed:\n  fresh:  {fh}\n  oracle: {oh}"
    assert len(fr) == len(orr), "row count differs from oracle"
    for a, b in zip(fr, orr):
        assert a[:4] == b[:4], f"integer counts differ: {a[:4]} vs oracle {b[:4]}"
        for x, y in zip(a[4:], b[4:]):  # precision / recall / f1
            assert abs(x - y) < 1e-9, f"metric differs: {a} vs oracle {b}"
