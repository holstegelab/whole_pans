#!/usr/bin/env python3
"""Select long reads reproducibly by name without writing another FASTQ."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fastq", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--coverage", required=True)
    parser.add_argument("--genome-size", required=True, type=int)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--read-names", required=True, type=Path)
    parser.add_argument("--stats", required=True, type=Path)
    return parser.parse_args()


def open_fastq(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def records(path: Path):
    with open_fastq(path) as handle:
        while True:
            name = handle.readline()
            if not name:
                return
            sequence = handle.readline().rstrip("\n")
            plus = handle.readline()
            quality = handle.readline().rstrip("\n")
            if not quality:
                raise ValueError(f"{path}: truncated FASTQ record near {name.strip()}")
            if not name.startswith("@") or not plus.startswith("+"):
                raise ValueError(f"{path}: malformed FASTQ record near {name.strip()}")
            if len(sequence) != len(quality):
                raise ValueError(f"{path}: sequence/quality length mismatch near {name.strip()}")
            yield name[1:].split()[0], len(sequence)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    args = parse_args()
    if not args.fastq.is_file() or args.fastq.stat().st_size == 0:
        raise FileNotFoundError(f"missing or empty FASTQ: {args.fastq}")
    total_reads = 0
    total_bases = 0
    for _, length in records(args.fastq):
        total_reads += 1
        total_bases += length
    if total_bases == 0:
        raise ValueError(f"no FASTQ bases in {args.fastq}")

    if args.coverage.lower() == "full":
        fraction = 1.0
        requested_bases = total_bases
    else:
        coverage = float(args.coverage)
        if coverage <= 0:
            raise ValueError("selected coverage must be positive or 'full'")
        requested_bases = int(round(coverage * args.genome_size))
        fraction = min(1.0, requested_bases / total_bases)

    args.read_names.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.read_names) + ".tmp")
    selected_reads = 0
    selected_bases = 0
    seen = set()
    with temporary.open("w") as output:
        for name, length in records(args.fastq):
            if name in seen:
                raise ValueError(f"duplicate FASTQ read name: {name}")
            seen.add(name)
            digest = hashlib.blake2b(
                f"{args.seed}\0{name}".encode(), digest_size=8,
            ).digest()
            value = int.from_bytes(digest, "big") / 2**64
            if value < fraction:
                output.write(name + "\n")
                selected_reads += 1
                selected_bases += length
    temporary.replace(args.read_names)
    if selected_reads == 0:
        raise ValueError("deterministic selection yielded zero reads")

    row = {
        "sample_id": args.sample,
        "input_fastq": str(args.fastq),
        "input_sha256": file_sha256(args.fastq),
        "selection_seed": args.seed,
        "genome_size_bp": args.genome_size,
        "requested_coverage": args.coverage,
        "requested_bases": requested_bases,
        "sampling_fraction": f"{fraction:.12g}",
        "total_reads": total_reads,
        "total_bases": total_bases,
        "selected_reads": selected_reads,
        "selected_bases": selected_bases,
        "selected_coverage": f"{selected_bases / args.genome_size:.8f}",
        "selection_method": "blake2b_seeded_read_name_threshold",
        "read_names_sha256": file_sha256(args.read_names),
    }
    fields = list(row)
    args.stats.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.stats) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(args.stats)


if __name__ == "__main__":
    main()
