#!/usr/bin/env python3
"""Codon-level variant generators for CARD coding sequences.

Three kinds of edit, all defined only on a valid coding sequence (in frame,
bacterial start codon, terminal stop codon, no internal stop):

  synonymous_recode  protein preserved      -- the "positive"
  frameshift         protein destroyed      -- a "negative"
  nonsense           protein truncated      -- a "negative"

A fourth, ``random_substitute``, is the control: it changes the same number of
nucleotides as a synonymous recoding but without preserving the protein, so
that flip rates under synonymous edits can be compared against flip rates under
edits of equal size.
"""
from __future__ import annotations

import random
from collections import defaultdict

# Standard genetic code (NCBI translation table 11, bacterial).  The table is
# identical to table 1 for the 64 codons; table 11 differs only in which codons
# may initiate translation.
CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

STOPS = {"TAA", "TAG", "TGA"}
# Bacteria initiate at ATG, GTG or TTG; all three are read as methionine when
# they are the initiation codon, which is why the start codon is never recoded.
STARTS = {"ATG", "GTG", "TTG"}
BASES = ("A", "C", "G", "T")

# amino acid -> list of codons encoding it
SYNONYMS: dict[str, list[str]] = defaultdict(list)
for _codon, _aa in CODON_TABLE.items():
    SYNONYMS[_aa].append(_codon)
SYNONYMS = dict(SYNONYMS)


def codons(seq: str) -> list[str]:
    """Split a sequence into codons.  Assumes len(seq) % 3 == 0."""
    return [seq[i:i + 3] for i in range(0, len(seq), 3)]


def translate(seq: str) -> str:
    """Translate a coding sequence to its amino-acid string, stopping at the
    first stop codon.  The initiation codon is read as methionine regardless of
    whether it is ATG, GTG or TTG."""
    out = []
    for index, codon in enumerate(codons(seq)):
        if index == 0 and codon in STARTS:
            out.append("M")
            continue
        amino = CODON_TABLE.get(codon)
        if amino is None:          # codon contains a non-ACGT character
            return "".join(out) + "?"
        if amino == "*":
            break
        out.append(amino)
    return "".join(out)


def is_valid_cds(seq: str) -> bool:
    """True if the sequence is a coding sequence on which the edits below are
    well defined: in frame, bacterial start, terminal stop, no internal stop."""
    if len(seq) < 6 or len(seq) % 3:
        return False
    cods = codons(seq)
    if any(c not in CODON_TABLE for c in cods):
        return False
    if cods[0] not in STARTS:
        return False
    if cods[-1] not in STOPS:
        return False
    return not any(c in STOPS for c in cods[:-1])


def recodable_positions(seq: str) -> list[int]:
    """Interior codon positions that have at least one alternative synonymous
    codon.  The start and stop codons are excluded: the start reads as
    methionine whatever its identity, and the stop is not an amino acid."""
    cods = codons(seq)
    return [
        i for i in range(1, len(cods) - 1)
        if len(SYNONYMS[CODON_TABLE[cods[i]]]) > 1
    ]


def synonymous_recode(seq: str, p: float = 0.3, rng: random.Random | None = None) -> str:
    """Replace a fraction ``p`` of the recodable interior codons with a
    different codon for the same amino acid.  The translated protein is
    unchanged by construction."""
    rng = rng or random.Random()
    cods = codons(seq)
    positions = recodable_positions(seq)
    if not positions:
        return seq
    k = max(1, round(p * len(positions)))
    for i in rng.sample(positions, min(k, len(positions))):
        alternatives = [c for c in SYNONYMS[CODON_TABLE[cods[i]]] if c != cods[i]]
        cods[i] = rng.choice(alternatives)      # never the original codon
    return "".join(cods)


def frameshift(seq: str, rng: random.Random | None = None) -> str:
    """Delete one nucleotide, shifting the reading frame for everything
    downstream.  The deletion site is drawn from the first 80% of the sequence
    so that a meaningful stretch of protein is actually disrupted."""
    rng = rng or random.Random()
    lo, hi = 3, max(4, int(0.8 * len(seq)))     # never inside the start codon
    cut = rng.randrange(lo, hi)
    return seq[:cut] + seq[cut + 1:]


def nonsense(seq: str, rng: random.Random | None = None) -> str:
    """Replace one interior codon with a stop codon, truncating the protein.
    Restricted to the first 80% of codons: a stop near the 3' end removes only
    a few residues and is not reliably disruptive."""
    rng = rng or random.Random()
    cods = codons(seq)
    hi = max(2, int(0.8 * (len(cods) - 1)))
    i = rng.randrange(1, hi)
    cods[i] = rng.choice(sorted(STOPS))
    return "".join(cods)


def random_substitute(seq: str, n_changes: int, rng: random.Random | None = None) -> str:
    """Control edit: change ``n_changes`` nucleotides at random interior
    positions, without regard to the protein.  Used to check whether a model is
    sensitive to synonymous edits specifically or merely to edits of that size."""
    rng = rng or random.Random()
    chars = list(seq)
    lo, hi = 3, len(seq) - 3                    # leave start and stop intact
    for i in rng.sample(range(lo, hi), min(n_changes, hi - lo)):
        chars[i] = rng.choice([b for b in BASES if b != chars[i]])
    return "".join(chars)


def n_substitutions(original: str, edited: str) -> int:
    """Number of differing nucleotides between two equal-length sequences."""
    return sum(a != b for a, b in zip(original, edited))