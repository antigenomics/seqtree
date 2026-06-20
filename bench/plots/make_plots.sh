#!/usr/bin/env bash
# Render figures from the committed benchmark tables (no seqtree run, no measuring).
#   bench/plots/make_plots.sh
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
python "$here/plot_retrieval_pr.py"
