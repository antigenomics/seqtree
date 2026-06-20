"""Performance regression: build/search time and peak RSS on a fixed seeded workload
must stay near the committed baseline. Time thresholds are generous (CI runners vary
several-fold from the dev machine); the RSS threshold is tight because memory is
machine-independent — it tracks the index data structures, so it catches real
memory regressions. Refresh the baseline with ``bench/tables/make_tables.sh --perf``."""
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PRODUCER = REPO / "bench" / "tables" / "gen_perf_table.py"
BASELINE = REPO / "bench" / "tables" / "perf_baseline.tsv"

# metric -> max allowed (measured / baseline). Time is loose, memory tight.
TOLERANCE = {"build_ms": 10.0, "search_ms": 8.0, "peak_rss_mb": 1.6}

# Timing is noisy on shared runners, so this runs only in the dedicated benchmarks
# job (RUN_PERF=1), not across the whole test matrix.
pytestmark = pytest.mark.skipif(not os.getenv("RUN_PERF"), reason="set RUN_PERF=1")


def _metrics(text):
    out = {}
    for line in text.splitlines():
        if line.startswith("#") or line.startswith("metric") or not line.strip():
            continue
        k, v = line.split("\t")
        out[k] = float(v)
    return out


@pytest.mark.skipif(not PRODUCER.exists() or not BASELINE.exists(), reason="perf baseline absent")
def test_perf_within_threshold():
    fresh = _metrics(subprocess.run([sys.executable, str(PRODUCER)], cwd=REPO,
                                    capture_output=True, text=True, check=True).stdout)
    base = _metrics(BASELINE.read_text())
    regressions = []
    for k, factor in TOLERANCE.items():
        if base.get(k, 0) <= 0:
            continue
        ratio = fresh[k] / base[k]
        print(f"{k}: {fresh[k]:.1f} vs baseline {base[k]:.1f}  (x{ratio:.2f}, limit x{factor})")
        if ratio > factor:
            regressions.append(f"{k} regressed x{ratio:.2f} > x{factor}")
    assert not regressions, "; ".join(regressions)
