# Release audit — 17 October 2026

The requested generic batch/search APIs are implemented. No new native algorithm is needed for the audited downstream requirements.

- Native batch search and edit histograms; aggregate insertion/deletion caps.
- Float64 global-alignment target reductions, separate candidate/final scoring semantics, fixed populations under exact exclusion.
- Grouped geometry, fixed audit partitions, bounded deterministic sparse radius hits and optional-geometry compatibility at the supported nanobind floor.
- Shared CPU allocation checks; one-thread inline execution and explicit over-budget rejection.

Verification on `0852ff5`: **379 passed, 1 Linux-only platform skip** across `test_batch_threads.py`, `test_alignment_max.py`, `test_edit_histogram.py`, `test_grouped_optional_geometry.py` and `test_grouped_sparse_hits.py`. An additional 280 checks cover grouped populations, withholding, paired lanes, offsets, score bins and partitions. Native budget was two CPUs; numerical libraries used one thread. The committed `examples/01_gapped_search.py` runs and its per-query cutoff/search output was inspected. Sphinx builds with `-W`.

The local `dev` commits `0126e00`, `858056b` and `0852ff5` were three commits ahead of the remote at audit start. They are preserved and included when this audited `dev` is pushed. The changelog now covers these capabilities; public API documentation already describes the alignment reductions. [Issue 18](https://github.com/antigenomics/seqtree/issues/18) tracks integration into the release head and the remaining exact-head cross-platform/wheel gate. Publishing remains the maintainer's action; this audit creates no tag or release. Focused checks do not substitute for complete release CI.
