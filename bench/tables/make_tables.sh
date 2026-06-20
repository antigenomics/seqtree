#!/usr/bin/env bash
# Regenerate the benchmark tables from the *compiled seqtree in the current repo*.
# The retrieval table is a deterministic oracle (checked byte-equivalent in CI);
# the perf table is a machine-dependent baseline (checked within a tolerance).
#
#   bench/tables/make_tables.sh            # regenerate all tables in place
#   bench/tables/make_tables.sh --perf     # also refresh the perf baseline
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"

python "$here/gen_retrieval_table.py" --out "$here/retrieval_pr.tsv"
echo "wrote $here/retrieval_pr.tsv"

if [ "${1:-}" = "--perf" ]; then
  python "$here/gen_perf_table.py" --out "$here/perf_baseline.tsv"
  echo "wrote $here/perf_baseline.tsv"
fi
