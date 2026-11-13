import sys
from pathlib import Path
from collections import Counter
sys.path.insert(0, "scripts")
from prepare_card_dataset import card_positives

STOPS = {"TAA", "TAG", "TGA"}
STARTS = {"ATG", "GTG", "TTG"} 

def status(seq):
    if len(seq) % 3: return "not_in_frame"
    codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    if codons[0] not in STARTS: return "no_start"
    if codons[-1] not in STOPS: return "no_stop"
    if any(c in STOPS for c in codons[:-1]): return "internal_stop"
    return "valid"

pos = card_positives(Path("data/raw"))
counts = Counter(status(s) for _, s in pos)
total = sum(counts.values())
print(f"CARD positives: {total}")
for st, n in counts.most_common():
    print(f"  {st:<14} {n:>5}  ({n/total:.1%})")
