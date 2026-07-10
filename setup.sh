#!/usr/bin/env fish
# Repo-local venv + editable install. Flags:
#   --tests  also install pytest extra
#   --bench  also install benchmark extras (huggingface_hub, pandas, psutil)
set repo (dirname (status --current-filename))
cd $repo

if not test -d .venv
    python3 -m venv .venv
end
source .venv/bin/activate.fish
pip install -U pip

set extras ""
if contains -- --tests $argv
    set extras "[test]"
else if contains -- --bench $argv
    set extras "[bench]"
end

pip install -e ".$extras"

echo ""
echo "Done. C++ tests:   cmake -S . -B build -G Ninja -DSEQTREE_TESTS=ON; cmake --build build; ctest --test-dir build"
echo "      Python tests: pytest tests/python"
echo "      Benchmarks:   python bench/bench.py   (1M tier: env RUN_BENCHMARK=1 python bench/bench.py --sizes 1000000 --queries 1000000)"
