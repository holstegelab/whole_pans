#!/usr/bin/env python3
"""Apply KEEP/REVIEW/SPLIT/REMOVE decisions without modifying source FASTA."""

import argparse
import csv
import gzip
import io
from collections import defaultdict
from contextlib import contextmanager
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--decisions", required=True)
    parser.add_argument("--split-bed", required=True)
    parser.add_argument("--cleaned", required=True)
    parser.add_argument("--removed", required=True)
    parser.add_argument("--review", required=True)
    parser.add_argument("--split-map", required=True)
    return parser.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


@contextmanager
def deterministic_gzip_text(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed) as text:
                yield text


def fasta_records(path):
    header = None
    chunks = []
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].rstrip("\n")
                chunks = []
            else:
                if header is None:
                    raise ValueError(f"Sequence before first FASTA header at {path}:{line_number}")
                chunks.append(line.strip())
    if header is not None:
        yield header, "".join(chunks)


def write_record(handle, header, sequence):
    handle.write(f">{header}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    start, end = sorted(intervals)[0]
    for next_start, next_end in sorted(intervals)[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            merged.append((start, end))
            start, end = next_start, next_end
    merged.append((start, end))
    return merged


def complement_intervals(length, removed):
    kept = []
    cursor = 0
    for start, end in removed:
        if cursor < start:
            kept.append((cursor, start))
        cursor = max(cursor, end)
    if cursor < length:
        kept.append((cursor, length))
    return kept


def main():
    args = parse_args()
    with open(args.decisions, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    decisions = {row["contig"]: row["decision"] for row in rows}
    if len(decisions) != len(rows):
        raise ValueError("Duplicate contig in decision table")

    split_intervals = defaultdict(list)
    with open(args.split_bed) as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed BED line {line_number}")
            split_intervals[fields[0]].append((int(fields[1]), int(fields[2])))
    split_intervals = {key: merge_intervals(value) for key, value in split_intervals.items()}

    Path(args.split_map).parent.mkdir(parents=True, exist_ok=True)
    observed = set()
    cleaned_records = 0
    with deterministic_gzip_text(args.cleaned) as cleaned, deterministic_gzip_text(
        args.removed
    ) as removed, deterministic_gzip_text(args.review) as review, open(
        args.split_map, "w", newline=""
    ) as map_handle:
        map_fields = [
            "original_contig",
            "original_length_bp",
            "output_contig",
            "start_1based",
            "end_1based",
            "length_bp",
            "disposition",
        ]
        map_writer = csv.DictWriter(
            map_handle, fieldnames=map_fields, delimiter="\t", lineterminator="\n"
        )
        map_writer.writeheader()

        for full_header, sequence in fasta_records(args.assembly):
            contig = full_header.split()[0]
            if contig in observed:
                raise ValueError(f"Duplicate FASTA identifier: {contig}")
            observed.add(contig)
            if contig not in decisions:
                raise ValueError(f"No decision for contig: {contig}")
            decision = decisions[contig]
            length = len(sequence)

            if decision in {"KEEP", "REVIEW"}:
                write_record(cleaned, full_header, sequence)
                cleaned_records += 1
                if decision == "REVIEW":
                    write_record(review, full_header, sequence)
                map_writer.writerow(
                    {
                        "original_contig": contig,
                        "original_length_bp": length,
                        "output_contig": contig,
                        "start_1based": 1,
                        "end_1based": length,
                        "length_bp": length,
                        "disposition": decision.lower(),
                    }
                )
            elif decision == "REMOVE":
                write_record(removed, f"{contig} decision=REMOVE", sequence)
                map_writer.writerow(
                    {
                        "original_contig": contig,
                        "original_length_bp": length,
                        "output_contig": "",
                        "start_1based": 1,
                        "end_1based": length,
                        "length_bp": length,
                        "disposition": "removed",
                    }
                )
            elif decision == "SPLIT":
                intervals = split_intervals.get(contig, [])
                if not intervals:
                    raise ValueError(f"SPLIT decision without intervals: {contig}")
                if any(start < 0 or end > length or start >= end for start, end in intervals):
                    raise ValueError(f"Invalid split interval for {contig}: {intervals}")

                for number, (start, end) in enumerate(intervals, start=1):
                    output_name = f"{contig}__nonhuman_{number}__{start + 1}_{end}"
                    write_record(removed, output_name, sequence[start:end])
                    map_writer.writerow(
                        {
                            "original_contig": contig,
                            "original_length_bp": length,
                            "output_contig": output_name,
                            "start_1based": start + 1,
                            "end_1based": end,
                            "length_bp": end - start,
                            "disposition": "removed_split_interval",
                        }
                    )
                for number, (start, end) in enumerate(
                    complement_intervals(length, intervals), start=1
                ):
                    if start == end:
                        continue
                    output_name = f"{contig}__clean_part{number}__{start + 1}_{end}"
                    write_record(cleaned, output_name, sequence[start:end])
                    cleaned_records += 1
                    map_writer.writerow(
                        {
                            "original_contig": contig,
                            "original_length_bp": length,
                            "output_contig": output_name,
                            "start_1based": start + 1,
                            "end_1based": end,
                            "length_bp": end - start,
                            "disposition": "kept_split_fragment",
                        }
                    )
            else:
                raise ValueError(f"Unknown decision for {contig}: {decision}")

    missing = set(decisions) - observed
    if missing:
        raise ValueError(f"Decision contigs absent from FASTA: {sorted(missing)[:5]}")
    if cleaned_records == 0:
        raise ValueError("Filtering removed the entire assembly")
    print(f"Wrote {cleaned_records} cleaned FASTA records")


if __name__ == "__main__":
    main()
