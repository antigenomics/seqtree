"""Background control repertoires for E-value calibration.

Small requests use the subset bundled with the package; larger ones stream from
the ``isalgo/airr_control`` HuggingFace dataset (needs ``huggingface_hub``), build
an index, and cache a serialized copy under ``~/.cache/seqtree`` for fast reuse.
Two properties are load-bearing, and both were once wrong here.

**The control must be productive.** ``appendix/evalue.tex`` (ass:match) requires the control and
the target to share the background law ``P0``. Targets are productive TCRs. VDJtools marks
out-of-frame rearrangements with ``_`` and in-frame stops with ``*``; both are non-coding. Worse,
out-of-frame junctions escape thymic selection and are the standard empirical proxy for the
*generation* law ``Pgen``, which Lemma "the two nulls differ" says is not ``P0``. 13.7% of the
mouse TRB table is out-of-frame; the human table is already productive-only.

**The control must be a uniform sample.** ``ass:indep`` assumes the unique clonotypes of the
control are i.i.d. from ``P0``. The upstream tables are sorted by clonotype abundance, so taking
the first ``size`` rows returns the most expanded clonotypes -- a public-clone head, not a sample.
``_download`` reservoir-samples uniformly over the unique clonotypes instead, and the bundled asset
is shuffled so that ``bundled[:size]`` is itself a valid sub-sample.

.. note::
   The out-of-frame clonotypes dropped here are not merely noise. They escape thymic selection, so
   they are the standard empirical proxy for the *generation* law ``Pgen`` -- exactly the analytic
   fallback Lemma "the two nulls differ" reaches for when a rare query's productive-control ball is
   empty. Recovering them needs the ``*.ntvj`` tables (the ``.aa`` table's ``_`` has already
   destroyed the residue count). Nothing calls for it yet; build it when something does.
"""
import gzip
import hashlib
import os
import random
import warnings
from importlib import resources

from ._core import Index, alphabet_symbols

#: The standard 20. A productive junction contains nothing else: not ``_`` (out of frame), not
#: ``*`` (stop), not ``X``/``B``/``Z`` (ambiguous, and ``pen(X, a)`` is barely half a mismatch).
_PRODUCTIVE_AA = frozenset("ACDEFGHIKLMNPQRSTVWY")

# name -> bundled asset filename
_BUNDLED = {"human_trb_aa": "control_human_trb_aa.txt.gz"}
# name -> (hf repo, file, cdr3-aa column)
_HF = {
    "human_trb_aa": ("isalgo/airr_control", "human.trb.aa.vdjtools.tsv.gz", "cdr3aa"),
    "human_tra_aa": ("isalgo/airr_control", "human.tra.aa.vdjtools.tsv.gz", "cdr3aa"),
    "mouse_trb_aa": ("isalgo/airr_control", "mouse.trb.aa.vdjtools.tsv.gz", "cdr3aa"),
    "mouse_tra_aa": ("isalgo/airr_control", "mouse.tra.aa.vdjtools.tsv.gz", "cdr3aa"),
}


#: Bump when the download or sampling pipeline changes what a given ``(name, size, seed)`` yields.
#: The bundled path needs no epoch -- it is fingerprinted from the asset's own bytes.
_DOWNLOAD_EPOCH = 1

#: Clonotypes in each bundled asset. Lets :func:`_cache_key` tell "served from the bundle" from
#: "served by a download" without decompressing the asset on every call; checked against the real
#: thing whenever the asset is actually read, so it cannot drift silently.
_BUNDLED_LEN = {"human_trb_aa": 250_000}


def _cache_dir(cache_dir):
    d = cache_dir or os.path.join(os.path.expanduser("~"), ".cache", "seqtree")
    os.makedirs(d, exist_ok=True)
    return d


def _asset_digest(name):
    """SHA-256 of the bundled asset's *compressed* bytes -- a few ms, no decompression."""
    h = hashlib.sha256()
    with resources.files("seqtree").joinpath("data", _BUNDLED[name]).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _cache_key(name, size, alphabet, seed):
    """A cache filename that changes whenever the index's content would.

    The old key was ``control_{name}_{size}.sqtree``, which named neither the **alphabet** nor the
    **seed** nor the **source data**. Three ways that went wrong:

    * Two calls differing only in ``seed`` -- which must draw different reservoir samples -- shared
      one cache file, so the second silently received the first's sequences.
    * The same for ``alphabet``.
    * An upgrade that changed the bundled asset kept the same filename, so a warm cache served the
      *previous release's* control. That is exactly how 0.3.0's corrected (uniform) control could
      be masked by a stale 0.2.0 (abundance-head) cache, and why its release notes had to ask
      people to delete ``~/.cache/seqtree`` by hand. Fingerprinting the asset makes a stale cache
      simply miss.

    ``seed`` enters the key only on the download path, because the bundled path takes a prefix of a
    pre-shuffled asset and ignores the seed entirely -- including it there would build byte-identical
    caches under different names.
    """
    from_bundle = name in _BUNDLED and (size is None or size <= _BUNDLED_LEN[name])
    parts = [name, str(size or "bundled"), alphabet]
    if from_bundle:
        parts.append(_asset_digest(name))
    else:
        repo, fname, _ = _HF[name]
        # There is no cheap *offline* way to learn the remote revision, and paying a network
        # round-trip on every warm-cache hit would be worse than the problem. Pin an epoch instead.
        parts += [repo, fname, f"epoch{_DOWNLOAD_EPOCH}", f"seed{seed}"]
    digest = hashlib.sha256("|".join(parts).encode()).hexdigest()[:12]
    return f"control_{name}_{size or 'bundled'}_{digest}.sqtree"


def _prune_superseded(cache_dir, name, size, keep):
    """Drop caches for the same control whose fingerprint no longer matches, and pre-fingerprint
    ones from older seqtree releases. Best-effort: a cache we cannot delete is only wasted disk."""
    prefix = f"control_{name}_{size or 'bundled'}"
    for f in os.listdir(cache_dir):
        if f == keep or not f.startswith(prefix) or not f.endswith(".sqtree"):
            continue
        rest = f[len(prefix):-len(".sqtree")]
        if rest == "" or (rest.startswith("_") and "_" not in rest[1:]):  # legacy, or a stale digest
            try:
                os.remove(os.path.join(cache_dir, f))
            except OSError:
                pass


class _NoLock:
    """Stand-in when ``filelock`` is absent. Correctness does not depend on the lock."""

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _build_lock(cache):
    """An inter-process lock around build-and-save, if ``filelock`` happens to be installed.

    This is an *optimisation*, not the fix. ``Index.save`` writes to a temporary and renames it
    into place, so a reader can never observe a half-written index whether or not this lock is
    taken. What the lock buys is work: without it, a cold-cache fan-out of N workers has every
    worker build the same 250k-clonotype index (seconds of CPU, 45 MB written) and N-1 of those
    are thrown away. With it, one builds and the rest wait and load.

    seqtree has no runtime dependencies, so ``filelock`` is used only if something else already
    put it there -- ``huggingface_hub`` pulls it in, and that is installed whenever a control is
    downloaded rather than read from the bundled asset.
    """
    try:
        from filelock import FileLock
    except ImportError:
        return _NoLock()
    return FileLock(cache + ".lock", timeout=600)


def _load_cached(cache):
    """The cached index, or ``None`` if there is nothing usable there.

    A cache is *best effort*: a file written by an older seqtree, truncated by a full disk, or
    left behind by a killed process must send us back to rebuilding it, never raise at the caller.
    """
    if not os.path.exists(cache):
        return None
    try:
        return Index.load(cache)
    except (RuntimeError, OSError):
        return None


def _read_bundled(name):
    with resources.files("seqtree").joinpath("data", _BUNDLED[name]).open("rb") as fh:
        with gzip.open(fh, "rt", encoding="utf-8") as gz:
            seqs = [line.strip() for line in gz if line.strip()]
    # _cache_key trusts this count to route bundle-vs-download without decompressing the asset on
    # every call, so it must never drift from the asset itself.
    if len(seqs) != _BUNDLED_LEN[name]:
        raise RuntimeError(f"bundled control '{name}' holds {len(seqs):,} sequences but "
                           f"_BUNDLED_LEN says {_BUNDLED_LEN[name]:,}; update it")
    return seqs


def sanitize(seqs, alphabet="aa"):
    """Keep only productive clonotypes. Drop non-coding and ambiguous ones.

    For ``"aa"`` this means the standard 20 residues and nothing else. ``_`` marks an out-of-frame
    rearrangement and cannot be repaired at the amino-acid level: VDJtools collapses a *run* of
    untranslatable positions into a single character (``replaceAll(/([atgc#~_?])+/, "_")``), so the
    residue count is already lost -- ``CAPM_QYF`` is eight characters standing for a 23-nucleotide
    junction. ``*`` marks an in-frame stop.

    .. warning::
       Do not relax this to keep ``_`` or ``*``. The control set *is* the E-value null and must
       share the target's background law. Out-of-frame junctions escape thymic selection and
       estimate ``Pgen``, not ``P0``; keeping them would swap the null for the pre-selection
       generation law and move ``M`` for mouse TRB from 694,241 to 816,382.

    Args:
        seqs: candidate sequences.
        alphabet: ``"aa"`` applies the productive 20-residue filter; other alphabets fall back to
            their symbol set.

    Returns:
        ``(kept, n_dropped)``.

    Example:
        >>> sanitize(["CASSLYEQYF", "C*A_FF", "CASSXYEQYF"])
        (['CASSLYEQYF'], 2)
    """
    ok = _PRODUCTIVE_AA if alphabet == "aa" else set(alphabet_symbols(alphabet))
    kept = [s for s in seqs if s and set(s) <= ok]
    return kept, len(seqs) - len(kept)


def _download(name, size, alphabet="aa", seed=0):
    """Stream a control table, keep productive clonotypes, and sample ``size`` of them uniformly.

    The upstream tables are sorted by clonotype abundance, so ``head -n size`` returns the most
    expanded public clones. Reservoir sampling gives every unique clonotype the same chance,
    which is what ``ass:indep`` assumes. Costs one full pass; the result is cached.
    """
    import csv

    from huggingface_hub import hf_hub_download

    repo, fname, col = _HF[name]
    path = hf_hub_download(repo_id=repo, filename=fname, repo_type="dataset")
    ok = _PRODUCTIVE_AA if alphabet == "aa" else set(alphabet_symbols(alphabet))
    rng = random.Random(seed)
    seen, out, dropped, n_seen = set(), [], 0, 0
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = [c.lower() for c in next(reader)]
        ci = header.index(col)
        for row in reader:
            if len(row) <= ci:
                continue
            s = row[ci].strip().upper()
            if not s:
                continue
            if not set(s) <= ok:
                dropped += 1          # '_' out of frame, '*' stop, ambiguous residues
                continue
            if s in seen:
                continue
            seen.add(s)
            if size is None:
                out.append(s)
                continue
            n_seen += 1               # reservoir over the UNIQUE clonotypes, not the rows
            if len(out) < size:
                out.append(s)
            else:
                j = rng.randrange(n_seen)
                if j < size:
                    out[j] = s
    if dropped:
        warnings.warn(f"control '{name}': dropped {dropped:,} non-productive or ambiguous "
                      f"clonotypes (VDJtools marks these '_' out-of-frame and '*' stop)")
    return out


def load_control(name="human_trb_aa", size=None, cache_dir=None, alphabet="aa", seed=0):
    """Build (or load from cache) an :class:`Index` over a background control set.

    The control is the E-value null. It is filtered to productive clonotypes (:func:`sanitize`)
    and, when subsampled, drawn uniformly over unique clonotypes rather than taken from the
    abundance-sorted head of the table.

    The built index is cached under ``cache_dir``. The cache is **content-addressed**: its filename
    carries a fingerprint of everything that determines the sequences -- the bundled asset's own
    bytes, or the download's source and seed -- so an upgrade that changes the control simply misses
    the old cache instead of silently serving it. Caches superseded that way are deleted.

    Safe to call from many processes at once on a cold cache; see :meth:`Index.save`.

    Args:
        name: control identifier (e.g. ``"human_trb_aa"``).
        size: number of unique clonotypes. ``None`` uses the full bundled subset;
            a value larger than the bundled subset triggers a HuggingFace download.
        cache_dir: where to store the serialized index (default ``~/.cache/seqtree``).
        alphabet: sequence alphabet for the index.
        seed: reservoir-sampling seed for the **download** path, so a given ``(name, size, seed)``
            is reproducible. The bundled path takes a prefix of a pre-shuffled asset and ignores it.

    Returns:
        An immutable ``Index`` of unique control clonotypes.
    """
    if name not in _BUNDLED and name not in _HF:
        raise ValueError(f"unknown control '{name}' (known: {sorted(set(_BUNDLED) | set(_HF))})")

    cdir = _cache_dir(cache_dir)
    cache = os.path.join(cdir, _cache_key(name, size, alphabet, seed))
    cached = _load_cached(cache)
    if cached is not None:
        return cached

    # Cold cache. Several processes may arrive here at once -- pytest-xdist, a Snakemake or
    # Nextflow fan-out, any multi-process pipeline sharing ~/.cache. Hold the lock, then look
    # again: whoever we were queued behind has very likely just built it.
    with _build_lock(cache):
        cached = _load_cached(cache)
        if cached is not None:
            return cached

        bundled = _read_bundled(name) if name in _BUNDLED else None
        if bundled is not None and (size is None or size <= len(bundled)):
            seqs = bundled if size is None else bundled[:size]
            seqs, dropped = sanitize(seqs, alphabet)
            if dropped:
                warnings.warn(f"bundled control '{name}': dropped {dropped:,} invalid sequences")
        else:
            if size is not None and size > 5_000_000:
                warnings.warn(f"downloading {size} control sequences may use several GB of memory")
            seqs = _download(name, size, alphabet, seed)

        seqs = list(dict.fromkeys(seqs))  # unique clonotypes, stable order
        idx = Index.build(seqs, alphabet=alphabet)
        try:
            idx.save(cache)  # temp file + rename: a reader never sees a partial index
            _prune_superseded(cdir, name, size, os.path.basename(cache))
        except (RuntimeError, OSError):
            pass  # cache is best-effort; a read-only or full ~/.cache must not fail the call
        return idx
