"""Pluggable index layout for epitopes/CDR3: anchor specs, masking, k-mer extraction.

Homology that drives TCR cross-reactivity is a shared *central, TCR-facing* motif
(Dolton et al., Cell 2023: one HLA-A*02:01 TCR sees EAAGIGILTV / LLLGIGILVL /
NLSALGIFST via x-x-x-A/G-I/L-G-I-x-x-x), not the MHC anchors. We therefore mask
anchor positions and search anchor-masked k-mers. Anchor positions are parametrized
(presets per class, overridable); ``MASK``/``GAP`` are spare symbols in the
amino-acid alphabet so masked positions are first-class characters.
"""
from dataclasses import dataclass, field

MASK = "X"  # anchor / TCR-facing wildcard (a real symbol in the AA alphabet)
GAP = "B"   # neutral gap-pad for fixed-width class-I frames


@dataclass(frozen=True)
class AnchorSpec:
    """1-based anchor positions; negatives count from the C-terminus (-1 == last)."""
    cls: str
    anchors: tuple = ()

    def resolve(self, length: int) -> set:
        """0-based anchor indices for a peptide of the given length."""
        out = set()
        for a in self.anchors:
            idx = (a - 1) if a > 0 else (length + a)
            if 0 <= idx < length:
                out.add(idx)
        return out


# Optimal defaults per class; override per allele/class as needed.
#  mhc1: P2 + C-terminus (PΩ) are the buried anchors.
#  mhc2: register-agnostic by default (the 9-mer core register is unknown), so no
#        fixed mask -- shared central k-mers match regardless of register/flank trim.
#  cdr3: conserved C...F/W ends carry no internal anchor mask.
DEFAULTS = {
    "mhc1": AnchorSpec("mhc1", (2, -1)),
    "mhc2": AnchorSpec("mhc2", ()),
    "cdr3": AnchorSpec("cdr3", ()),
}


def spec_for(cls: str, override: AnchorSpec | None = None) -> AnchorSpec:
    if override is not None:
        return override
    if cls not in DEFAULTS:
        raise ValueError(f"unknown class '{cls}' (use 'mhc1', 'mhc2', or 'cdr3')")
    return DEFAULTS[cls]


def mask_anchors(pep: str, spec: AnchorSpec) -> str:
    """Replace anchor positions with MASK so they don't drive homology."""
    idx = spec.resolve(len(pep))
    return "".join(MASK if i in idx else c for i, c in enumerate(pep))


def kmers(pep: str, k: int, spec: AnchorSpec | None = None) -> list[str]:
    """Anchor-masked contiguous k-mers. Peptides shorter than k yield the whole
    (masked) peptide as a single token so they remain searchable."""
    s = mask_anchors(pep, spec) if spec else pep
    if len(s) < k:
        return [s] if s else []
    return [s[i:i + k] for i in range(len(s) - k + 1)]


def presentation_features(pep: str, cls: str) -> list[str]:
    """Short binding-motif signatures (anchor / pocket residues) for the reverse
    problem -- peptide -> presenting allele. These KEEP the anchors and drop the
    TCR-facing positions (the opposite of :func:`kmers`), so peptides binding the
    same allele share a signature. A peptide may yield several (class-II registers).

    class I: N-pocket P1-P3 + C-pocket P(Ω-1),PΩ -> one 5-residue signature.
    class II: core anchors P1,P4,P6,P9 over every 9-mer window -> one per register.
    """
    p = pep.strip().upper()
    L = len(p)
    if cls == "mhc1":
        if L < 4:
            return [p]
        return [p[0] + p[1] + p[2] + p[-2] + p[-1]]
    if cls == "mhc2":
        if L < 9:
            return [p]
        feats = [p[s] + p[s + 3] + p[s + 5] + p[s + 8] for s in range(L - 9 + 1)]
        return list(dict.fromkeys(feats))
    return [p]


def weight_profile(length: int, spec: AnchorSpec, mode: str = "tcr_facing",
                   hotspot: tuple = (), hotspot_weight: int = 1) -> list[int]:
    """Per-position weights for a PositionalMatrix over a length-`length` frame.

    mode='tcr_facing': anchors -> 0 (free), others -> 1 (optionally up-weight a
    central hotspot, e.g. class-I P4-P7). mode='presentation': the inverse --
    anchors -> 1, TCR-facing positions -> 0 -- for allele assignment.
    """
    anchors = spec.resolve(length)
    hot = {(h - 1) if h > 0 else (length + h) for h in hotspot}
    w = []
    for i in range(length):
        is_anchor = i in anchors
        if mode == "presentation":
            w.append(1 if is_anchor else 0)
        else:  # tcr_facing
            if is_anchor:
                w.append(0)
            elif i in hot:
                w.append(hotspot_weight)
            else:
                w.append(1)
    return w
