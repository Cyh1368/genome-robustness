import csv
from collections import Counter

STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"}

def status(seq):
    if len(seq) % 3: return "not_in_frame"
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons[0] not in STARTS: return "no_start"
    if codons[-1] not in STOPS: return "no_stop"
    if any(c in STOPS for c in codons[:-1]): return "internal_stop"
    return "valid"

for split in ["train", "dev", "test"]:
    counts = Counter()
    with open(f"data/processed/{split}.csv") as f:
        for row in csv.DictReader(f):
            counts[(row["label"], status(row["sequence"]))] += 1
    total = sum(counts.values())
    print(f"\n{split} (n={total})")
    for (label, st), n in sorted(counts.items()):
        print(f"  label={label:>1}  {st:<14} {n:>5}  ({n/total:.1%})")
