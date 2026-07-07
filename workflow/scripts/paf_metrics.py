#!/usr/bin/env python3
"""Summarize primary/supplementary assembly-to-reference PAF alignments."""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", required=True)
    parser.add_argument("--query-stats", required=True)
    parser.add_argument("--reference-stats", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--min-mapq", type=int, default=5)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def seqkit_total(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 1 or "sum_len" not in rows[0]:
        raise ValueError(f"Expected one SeqKit data row with sum_len in {path}")
    return int(rows[0]["sum_len"].replace(",", ""))


def merged_length(intervals):
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    start, end = intervals[0]
    for next_start, next_end in intervals[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def main():
    args = parse_args()
    query_lengths = {}
    target_lengths = {}
    query_intervals = defaultdict(list)
    target_intervals = defaultdict(list)
    query_names = set()
    target_names = set()
    matches = 0
    alignment_bases = 0
    alignment_count = 0

    opener = gzip.open if args.paf.endswith(".gz") else open
    with opener(args.paf, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number}: expected at least 12 fields")

            query, query_length = fields[0], int(fields[1])
            qstart, qend = int(fields[2]), int(fields[3])
            target, target_length = fields[5], int(fields[6])
            tstart, tend = int(fields[7]), int(fields[8])
            residue_matches, block_length = int(fields[9]), int(fields[10])
            mapq = int(fields[11])

            query_lengths[query] = query_length
            target_lengths[target] = target_length
            if mapq < args.min_mapq:
                continue

            query_intervals[query].append((qstart, qend))
            target_intervals[target].append((tstart, tend))
            query_names.add(query)
            target_names.add(target)
            matches += residue_matches
            alignment_bases += block_length
            alignment_count += 1

    # PAF omits completely unaligned query and target sequences. Use full FASTA
    # totals from SeqKit so coverage cannot be inflated by those omissions.
    query_total = seqkit_total(args.query_stats)
    target_total = seqkit_total(args.reference_stats)
    query_aligned = sum(merged_length(intervals) for intervals in query_intervals.values())
    target_covered = sum(merged_length(intervals) for intervals in target_intervals.values())

    row = {
        "reference": args.reference,
        "alignment_count": alignment_count,
        "query_contigs_aligned": len(query_names),
        "reference_contigs_covered": len(target_names),
        "query_total_bp": query_total,
        "query_aligned_bp": query_aligned,
        "query_aligned_percent": 100.0 * query_aligned / query_total if query_total else 0.0,
        "reference_total_bp": target_total,
        "reference_covered_bp": target_covered,
        "reference_covered_percent": 100.0 * target_covered / target_total if target_total else 0.0,
        "alignment_identity_percent": 100.0 * matches / alignment_bases if alignment_bases else 0.0,
    }

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
