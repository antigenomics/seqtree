# CLAUDE.md — working in this repo

`seqtree` is a payload-agnostic fuzzy sequence-search core: C++20 + nanobind, **zero runtime
dependencies**, no domain data. What the library *exposes* is in `skills/seqtree/SKILL.md` — read
that first, it is the API reference and the gotcha list. This file is only how to work here.

## Where things are

| | |
|---|---|
| `src/`, `include/seqtree/` | the C++ core; `src/_bindings.cpp` is the whole nanobind surface |
| `python/seqtree/` | pure-Python layers (`evalue`, `gapblock`, `distance`, `pmhc`, `layout`, `seeds`) |
| `tests/cpp/`, `tests/python/` | doctest and pytest |
| `bench/` | benchmarks; `bench/tables/` holds the committed gate baselines |
| `docs/`, `docs/design/` | Sphinx site; `docs/design/*.md` are the design docs with their correction logs |
| `appendix/evalue.tex` | the E-value theory — treat it as the spec, not the code |
| `skills/seqtree/SKILL.md` | API reference for agents |
| `ROADMAP.md` | the contract with downstream `vdjmatch` / `mhcmatch` |
| `docs/roadmap.rst` | seqtree's own internal roadmap |

## Build and test

```bash
bash setup.sh --tests --bench
cmake -S . -B build -G Ninja -DSEQTREE_TESTS=ON -DCMAKE_BUILD_TYPE=Release && cmake --build build --parallel
ctest --test-dir build --output-on-failure
pytest tests/python -q
```

- **`editable.rebuild = false`** (pyproject). Editing C++ and running pytest tests the *old*
  `_core.so` — reinstall (`uv pip install -e .`) after every C++ change or the result is a lie.
- Ninja occasionally misses a same-second edit; `touch src/foo.cpp` if it says "no work to do".
- `tests/cpp/test_codec.cpp` owns `DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN`. New test files must not
  define a main, and go in the flat `add_executable` list in `CMakeLists.txt`.
- Gates: `RUN_PERF=1 pytest tests/python/test_perf_regression.py tests/python/test_oracle_tables.py -q -s`.
  Regenerate with `bench/tables/make_tables.sh --perf` — **on an idle machine**. A baseline taken
  while another benchmark saturates the cores reads as a 1.5x regression on untouched code.
- Heavy benchmarks (proteome downloads, 1M tiers) are gated behind `RUN_BENCHMARK=1`.
- Docs are built with `-W` in CI (`sphinx-build -b html -W --keep-going docs docs/_build/html`);
  a warning is a failure. Locally: `env -C docs make html`.
- Windows CI is not optional — the `mmap` paths in `text_index.cpp` are the only platform-specific
  code in the tree.

## Conventions that aren't visible in the code

- **Version is single-sourced** from `pyproject.toml`; `__version__` reads distribution metadata.
  Never add a literal.
- **Gitflow**: feature branch → `dev` → `master`. Commit messages end with the `Co-Authored-By`
  trailer.
- **Never publish to PyPI without an explicit release** (`ROADMAP.md` §4). Pushing a GitHub Release
  fires `publish.yml` → OIDC upload. Tagging is fine; releasing is the owner's click.
- **seqtree stays generic.** No TCR germline data, no VDJdb, no clustering, no structures — those
  live downstream in `vdjmatch` / `mhcmatch` / `tcren`. Don't import domain data into `bench/`.
- **Never fabricate a citation** — verify every DOI with a tool before it goes in a doc or comment.
- `docs/design/*.md` carry a status block listing corrections and *rejected* optimisations with
  reasons. When a design decision is superseded, append a correction there rather than editing the
  body — the rejections are load-bearing (they stop the same optimisation being re-proposed).

## Open loops

- **`v1.0.0` is tagged on `master` (`fc9dde4`) but not released.** CI is green; nothing is on PyPI.
  The GitHub Release is the owner's to publish.
- **`TextIndex` at `L < 2k`** genuinely falls back to one block and a full radius-`max_subs` ball —
  documented under *Limits*, not fixed. The fix would be spaced/gapped seeds (`docs/design/text_index.md`
  §3.7), deferred by an earlier decision.
- Next design items are in `docs/roadmap.rst` (PSSM-graded d_TCR, native local alignment, Flashback
  build); downstream wrapper design is in `ROADMAP.md` §2–3 with `(TBD: owner)` markers.
