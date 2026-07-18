#!/usr/bin/env fish
# Repo-local venv + editable install, via uv. Flags:
#   --tests  also install pytest extra
#   --bench  also install benchmark extras (huggingface_hub, psutil)
#
# Needs uv (https://docs.astral.sh/uv/); `brew install uv` or the standalone installer.
set repo (dirname (status --current-filename))
cd $repo

if not command -q uv
    echo "setup.sh needs uv -- install it with 'brew install uv' or see https://docs.astral.sh/uv/"
    exit 1
end

# uv creates and manages .venv. `uv pip` finds .venv in the cwd without an active
# venv, but we activate so the pytest/cmake commands echoed below use it too.
if not test -d .venv
    uv venv
end
source .venv/bin/activate.fish

set extras ""
if contains -- --tests $argv
    set extras "[test]"
else if contains -- --bench $argv
    set extras "[bench]"
end

uv pip install -e ".$extras"

echo ""
echo "Done. C++ tests:   cmake -S . -B build -G Ninja -DSEQTREE_TESTS=ON; cmake --build build; ctest --test-dir build"
echo "      Python tests: pytest tests/python"
echo "      Benchmarks:   python bench/bench.py   (1M tier: env RUN_BENCHMARK=1 python bench/bench.py --sizes 1000000 --queries 1000000)"
