# Shared helpers for the bench scripts. 2026-09-06
#
# These lived as six near-copies of peak-RSS and three of mutate, which drifted: two of the RSS
# copies used psutil and reported *current* RSS rather than peak, so they under-reported the
# number the others measured. stdlib `resource` gives the peak directly and is one less
# benchmark dependency.
import random
import resource
import sys

#: The 20 standard residues, in the order the bench scripts have always used.
AA = "ACDEFGHIKLMNPQRSTVWY"


def peak_rss_mb() -> float:
    """Peak resident set size of this process, in MiB.

    ``ru_maxrss`` is bytes on macOS and kibibytes on Linux -- the one portability wart, and the
    reason this is worth having in one place.
    """
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return rss / (1024 * 1024) if sys.platform == "darwin" else rss / 1024


def peak_rss_gb() -> float:
    """Peak resident set size of this process, in GiB."""
    return peak_rss_mb() / 1024


def mutate(seq: str, n_subs: int, rng: random.Random, n_indels: int = 0) -> str:
    """``seq`` with ``n_subs`` random substitutions and ``n_indels`` random indels applied.

    Substitutions are drawn uniformly from :data:`AA` and may re-draw the residue already
    there, so the realised edit distance is a lower bound on ``n_subs`` -- which is what the
    benchmarks want, since it matches how a real repertoire varies.
    """
    if not seq:
        return seq
    s = list(seq)
    for _ in range(n_subs):
        s[rng.randrange(len(s))] = rng.choice(AA)
    for _ in range(n_indels):
        j = rng.randrange(len(s))
        if rng.random() < 0.5 and len(s) > 1:
            del s[j]
        else:
            s.insert(j, rng.choice(AA))
    return "".join(s)
