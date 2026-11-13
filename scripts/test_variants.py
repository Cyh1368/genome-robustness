#!/usr/bin/env python3
"""Check that the variant generators preserve or destroy the protein as claimed.

Every number in this project rests on the assumption that a synonymous recoding
leaves the translated protein untouched and that a frameshift or nonsense edit
does not.  This script verifies both on real CARD sequences before any model is
run.

Usage:
    python scripts/test_variants.py data/processed/test.csv
"""
from __future__ import annotations

import csv
import random
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from variants import (
    codons,
    frameshift,
    is_valid_cds,
    nonsense,
    n_substitutions,
    random_substitute,
    recodable_positions,
    synonymous_recode,
    translate,
)


def main(path: str, n_test: int = 500, p: float = 0.3, seed: int = 0) -> None:
    rng = random.Random(seed)

    with open(path) as handle:
        rows = [(row["sequence"], int(row["label"])) for row in csv.DictReader(handle)]

    valid = [(s, y) for s, y in rows if is_valid_cds(s)]
    print(f"{path}: {len(valid)} / {len(rows)} rows are valid CDS "
          f"({len(valid) / len(rows):.1%})")
    print(f"  by label: {Counter(y for _, y in valid)}")

    sample = rng.sample(valid, min(n_test, len(valid)))
    failures = Counter()
    changed_syn, changed_nonsense_len = [], []

    for seq, _ in sample:
        protein = translate(seq)

        # 1. Synonymous recoding must leave the protein identical.
        recoded = synonymous_recode(seq, p, rng)
        if translate(recoded) != protein:
            failures["synonymous_changed_protein"] += 1
        if len(recoded) != len(seq):
            failures["synonymous_changed_length"] += 1
        if recoded == seq and recodable_positions(seq):
            failures["synonymous_no_op"] += 1
        changed_syn.append(n_substitutions(seq, recoded))

        # 2. Frameshift must change the protein and shorten the sequence by one.
        shifted = frameshift(seq, rng)
        if len(shifted) != len(seq) - 1:
            failures["frameshift_wrong_length"] += 1
        if translate(shifted) == protein:
            failures["frameshift_kept_protein"] += 1

        # 3. Nonsense must truncate the protein.
        truncated = nonsense(seq, rng)
        if len(truncated) != len(seq):
            failures["nonsense_changed_length"] += 1
        if translate(truncated) == protein:
            failures["nonsense_kept_protein"] += 1
        if len(translate(truncated)) >= len(protein):
            failures["nonsense_did_not_shorten"] += 1
        changed_nonsense_len.append(len(translate(truncated)) / max(1, len(protein)))

        # 4. Control edit must match the synonymous edit in size.
        control = random_substitute(seq, changed_syn[-1], rng)
        if n_substitutions(seq, control) != changed_syn[-1]:
            failures["control_wrong_size"] += 1

    print(f"\nTested {len(sample)} sequences (p={p}, seed={seed})")
    if failures:
        print("FAILURES:")
        for name, count in failures.most_common():
            print(f"  {name:<28} {count}")
        sys.exit(1)

    mean_syn = sum(changed_syn) / len(changed_syn)
    mean_len = sum(changed_nonsense_len) / len(changed_nonsense_len)
    median_codons = sorted(len(codons(s)) for s, _ in sample)[len(sample) // 2]
    print("All invariants hold.")
    print(f"  median sequence length        {median_codons} codons")
    print(f"  mean nucleotides changed      {mean_syn:.1f} per synonymous recoding")
    print(f"  mean protein retained         {mean_len:.1%} after a nonsense edit")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "data/processed/test.csv")