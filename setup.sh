#!/bin/sh
# Repo-local venv + editable install, via uv. POSIX sh -- runs under bash and zsh. Flags:
#   --tests  also install pytest extra
#   --bench  also install benchmark extras (huggingface_hub, psutil)
#
# Needs uv (https://docs.astral.sh/uv/); `brew install uv` or the standalone installer.
cd "$(dirname "$0")" || exit 1

if ! command -v uv >/dev/null 2>&1; then
    echo "setup.sh needs uv -- install it with 'brew install uv' or see https://docs.astral.sh/uv/"
    exit 1
fi

# uv creates and manages .venv. `uv pip` finds .venv in the cwd without an active
# venv, but we activate so the pytest/cmake commands echoed below use it too.
if [ ! -d .venv ]; then
    uv venv
fi
. .venv/bin/activate

extras=""
for arg in "$@"; do
    case "$arg" in
        --tests) extras="[test]" ;;
        --bench) extras="[bench]" ;;
    esac
done

uv pip install -e ".$extras"

echo ""
echo "Done. C++ tests:   cmake -S . -B build -G Ninja -DSEQTREE_TESTS=ON; cmake --build build; ctest --test-dir build"
echo "      Python tests: pytest tests/python"
echo "      Benchmarks:   python bench/bench.py   (1M tier: env RUN_BENCHMARK=1 python bench/bench.py --sizes 1000000 --queries 1000000)"
