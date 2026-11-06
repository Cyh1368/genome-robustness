#!/usr/bin/env python3
"""Build a leakage-aware CARD AMR-gene classification dataset.

Positive examples are CARD nucleotide reference sequences.  Negative examples
are coding sequences from four non-pathogenic bacterial reference genomes.  A
gene-family split keeps any individual CARD ARO term in only one partition.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
import random
import re
import shutil
import tarfile
import urllib.request
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

CARD_URL = "https://card.mcmaster.ca/latest/data"
GENOMES = {
    "ecoli_k12": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2",
    "salmonella_lt2": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/945/GCF_000006945.2_ASM694v2/GCF_000006945.2_ASM694v2",
    "pseudomonas_pao1": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/765/GCF_000006765.1_ASM676v1/GCF_000006765.1_ASM676v1",
    "bacillus_168": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/045/GCF_000009045.1_ASM904v1/GCF_000009045.1_ASM904v1",
}


def fetch(url: str, destination: Path) -> None:
    if destination.exists():
        return
    print(f"Downloading {url}")
    with urllib.request.urlopen(url) as source, destination.open("wb") as target:
        shutil.copyfileobj(source, target)


def card_positives(raw: Path) -> list[tuple[str, str]]:
    archive = raw / "card-data.tar.bz2"
    fetch(CARD_URL, archive)
    with tarfile.open(archive, "r:bz2") as bundle:
        fasta_names = [name for name in bundle.getnames() if name.startswith("./nucleotide_fasta_")]
        positives = []
        for name in fasta_names:
            with bundle.extractfile(name) as handle:
                assert handle is not None
                for record in SeqIO.parse(io.TextIOWrapper(handle), "fasta"):
                    sequence = re.sub("[^ACGT]", "", str(record.seq).upper())
                    # CARD headers encode the stable term as ``ARO:########``.
                    match = re.search(r"ARO:\d+", record.description)
                    if match is None:
                        continue
                    aro = match.group(0)
                    if 90 <= len(sequence) <= 2500:
                        positives.append((aro, sequence))
    # Keep one copy of duplicate alleles and avoid unusually common AROs dominating.
    unique = {(aro, sequence) for aro, sequence in positives}
    by_aro: dict[str, list[str]] = defaultdict(list)
    for aro, sequence in unique:
        by_aro[aro].append(sequence)
    rng = random.Random(42)
    return [(aro, sequence) for aro, seqs in by_aro.items() for sequence in rng.sample(seqs, min(20, len(seqs)))]


def genome_controls(raw: Path) -> list[str]:
    controls = []
    for name, base in GENOMES.items():
        fna = raw / f"{name}.fna.gz"
        gff = raw / f"{name}.gff.gz"
        fetch(base + "_genomic.fna.gz", fna)
        fetch(base + "_genomic.gff.gz", gff)
        records = {record.id: record.seq for record in SeqIO.parse(gzip.open(fna, "rt"), "fasta")}
        for row in csv.reader(gzip.open(gff, "rt"), delimiter="\t"):
            if len(row) != 9 or row[2] != "CDS":
                continue
            sequence = records.get(row[0])
            if sequence is None:
                continue
            start, end = int(row[3]) - 1, int(row[4])
            gene = sequence[start:end]
            if row[6] == "-":
                gene = gene.reverse_complement()
            text = re.sub("[^ACGT]", "", str(gene).upper())
            if 90 <= len(text) <= 2500:
                controls.append(text)
    return list(set(controls))


def write_split(path: Path, rows: list[tuple[str, int]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["sequence", "label"])
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    positives = card_positives(args.raw_dir)
    controls = genome_controls(args.raw_dir)
    rng = random.Random(42)
    rng.shuffle(controls)
    if len(controls) < len(positives):
        raise RuntimeError("Not enough non-AMR coding-sequence controls")
    # Assign all alleles of a CARD ARO term together to prevent allele leakage.
    groups = defaultdict(list)
    for aro, sequence in positives:
        groups[aro].append(sequence)
    aro_ids = list(groups)
    rng.shuffle(aro_ids)
    cut1, cut2 = int(.70 * len(aro_ids)), int(.85 * len(aro_ids))
    split_ids = {"train": aro_ids[:cut1], "dev": aro_ids[cut1:cut2], "test": aro_ids[cut2:]}
    offset = 0
    for split, ids in split_ids.items():
        pos = [sequence for aro in ids for sequence in groups[aro]]
        neg = controls[offset:offset + len(pos)]
        offset += len(pos)
        rows = [(sequence, 1) for sequence in pos] + [(sequence, 0) for sequence in neg]
        rng.shuffle(rows)
        write_split(args.output_dir / f"{split}.csv", rows)
        print(f"{split}: {len(pos)} CARD positives, {len(neg)} controls")


if __name__ == "__main__":
    main()
