#!/usr/bin/env python3
"""Deterministically simulate truth-labelled HiFi-like reads from phased FASTAs."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import random
from pathlib import Path


COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--hap1", required=True, type=Path)
    parser.add_argument("--hap2", required=True, type=Path)
    parser.add_argument("--coverage", type=float, required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--read-length-mean", type=float, required=True)
    parser.add_argument("--read-length-sd", type=float, required=True)
    parser.add_argument("--substitution-rate", type=float, required=True)
    parser.add_argument("--insertion-rate", type=float, required=True)
    parser.add_argument("--deletion-rate", type=float, required=True)
    parser.add_argument("--empirical-profile-json", default="")
    parser.add_argument("--fastq", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def fasta_records(path: Path):
    name = None
    chunks = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks).upper()


def mutate(sequence: str, rng: random.Random, rates: dict[str, float]) -> str:
    bases = "ACGT"
    output = []
    for base in sequence:
        if rng.random() < rates["deletion"]:
            continue
        if base in bases and rng.random() < rates["substitution"]:
            output.append(rng.choice(bases.replace(base, "")))
        else:
            output.append(base if base in bases else "N")
        if rng.random() < rates["insertion"]:
            output.append(rng.choice(bases))
    return "".join(output)


def main() -> None:
    args = parse_args()
    profile = {
        "read_length_mean": args.read_length_mean,
        "read_length_sd": args.read_length_sd,
        "substitution_rate": args.substitution_rate,
        "insertion_rate": args.insertion_rate,
        "deletion_rate": args.deletion_rate,
    }
    profile_source = "configured_rates"
    if args.empirical_profile_json:
        with Path(args.empirical_profile_json).open() as handle:
            profile.update(json.load(handle))
        profile_source = args.empirical_profile_json
    rates = {
        "substitution": float(profile["substitution_rate"]),
        "insertion": float(profile["insertion_rate"]),
        "deletion": float(profile["deletion_rate"]),
    }
    rng = random.Random(args.seed)
    args.fastq.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.fastq) + ".tmp")
    read_count = output_bases = source_bases = 0
    with gzip.open(temporary, "wt", compresslevel=6) as output:
        for hap, fasta in (("hap1", args.hap1), ("hap2", args.hap2)):
            for contig, sequence in fasta_records(fasta):
                source_bases += len(sequence)
                expected = len(sequence) * (args.coverage / 2.0) / float(profile["read_length_mean"])
                count = int(expected)
                if rng.random() < expected - count:
                    count += 1
                for _ in range(count):
                    length = max(1000, int(rng.gauss(float(profile["read_length_mean"]), float(profile["read_length_sd"]))))
                    length = min(length, len(sequence))
                    if length <= 0:
                        continue
                    start = rng.randrange(0, len(sequence) - length + 1)
                    end = start + length
                    strand = "+" if rng.random() < 0.5 else "-"
                    source = sequence[start:end]
                    if strand == "-":
                        source = source.translate(COMPLEMENT)[::-1]
                    read = mutate(source, rng, rates)
                    name = f"sim|{args.sample}|{hap}|{contig}:{start}-{end}|{strand}|{read_count:09d}"
                    output.write(f"@{name}\n{read}\n+\n{'I' * len(read)}\n")
                    read_count += 1
                    output_bases += len(read)
    temporary.replace(args.fastq)
    row = {
        "sample_id": args.sample,
        "seed": args.seed,
        "requested_diploid_coverage": args.coverage,
        "source_assembly_bases": source_bases,
        "simulated_reads": read_count,
        "simulated_bases": output_bases,
        "profile_source": profile_source,
        "read_name_schema": "sim|sample|haplotype|contig:start-end|strand|ordinal",
    }
    args.stats.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.stats) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(args.stats)


if __name__ == "__main__":
    main()
