"""Background control repertoires for E-value calibration.

Small requests use the subset bundled with the package; larger ones stream from
the ``isalgo/airr_control`` HuggingFace dataset (needs ``huggingface_hub``), build
an index, and cache a serialized copy under ``~/.cache/seqtree`` for fast reuse.
Sequences are deduplicated to unique clonotypes (the E-value null assumes this).
"""
import gzip
import os
import warnings
from importlib import resources

from ._core import Index

# name -> bundled asset filename
_BUNDLED = {"human_trb_aa": "control_human_trb_aa.txt.gz"}
# name -> (hf repo, file, cdr3-aa column)
_HF = {
    "human_trb_aa": ("isalgo/airr_control", "human.trb.aa.vdjtools.tsv.gz", "cdr3aa"),
    "human_tra_aa": ("isalgo/airr_control", "human.tra.aa.vdjtools.tsv.gz", "cdr3aa"),
    "mouse_trb_aa": ("isalgo/airr_control", "mouse.trb.aa.vdjtools.tsv.gz", "cdr3aa"),
    "mouse_tra_aa": ("isalgo/airr_control", "mouse.tra.aa.vdjtools.tsv.gz", "cdr3aa"),
}


def _cache_dir(cache_dir):
    d = cache_dir or os.path.join(os.path.expanduser("~"), ".cache", "seqtree")
    os.makedirs(d, exist_ok=True)
    return d


def _read_bundled(name):
    with resources.files("seqtree").joinpath("data", _BUNDLED[name]).open("rb") as fh:
        with gzip.open(fh, "rt", encoding="utf-8") as gz:
            return [line.strip() for line in gz if line.strip()]


def _download(name, size):
    import csv

    from huggingface_hub import hf_hub_download

    repo, fname, col = _HF[name]
    path = hf_hub_download(repo_id=repo, filename=fname, repo_type="dataset")
    seen, out = set(), []
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = [c.lower() for c in next(reader)]
        ci = header.index(col)
        for row in reader:
            if len(row) <= ci:
                continue
            s = row[ci].strip().upper()
            if s and s not in seen:
                seen.add(s)
                out.append(s)
                if size and len(out) >= size:
                    break
    return out


def load_control(name="human_trb_aa", size=None, cache_dir=None, alphabet="aa"):
    """Build (or load from cache) an :class:`Index` over a background control set.

    Args:
        name: control identifier (e.g. ``"human_trb_aa"``).
        size: number of unique clonotypes. ``None`` uses the full bundled subset;
            a value larger than the bundled subset triggers a HuggingFace download.
        cache_dir: where to store the serialized index (default ``~/.cache/seqtree``).
        alphabet: sequence alphabet for the index.

    Returns:
        An immutable ``Index`` of unique control clonotypes.
    """
    if name not in _BUNDLED and name not in _HF:
        raise ValueError(f"unknown control '{name}' (known: {sorted(set(_BUNDLED) | set(_HF))})")

    cache = os.path.join(_cache_dir(cache_dir), f"control_{name}_{size or 'bundled'}.sqtree")
    if os.path.exists(cache):
        return Index.load(cache)

    bundled = _read_bundled(name) if name in _BUNDLED else None
    if bundled is not None and (size is None or size <= len(bundled)):
        seqs = bundled if size is None else bundled[:size]
    else:
        if size is not None and size > 5_000_000:
            warnings.warn(f"downloading {size} control sequences may use several GB of memory")
        seqs = _download(name, size)

    seqs = list(dict.fromkeys(seqs))  # unique clonotypes, stable order
    idx = Index.build(seqs, alphabet=alphabet)
    try:
        idx.save(cache)
    except OSError:
        pass  # cache is best-effort
    return idx
