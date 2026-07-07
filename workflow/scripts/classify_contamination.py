#!/usr/bin/env python3
"""Classify contigs using non-human BLAST hits and human-reference support."""

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path

import yaml


DECISION_FIELDS = [
    "assembly_id",
    "contig",
    "contig_length_bp",
    "decision",
    "nonhuman_covered_bp",
    "nonhuman_covered_percent",
    "human_covered_bp",
    "human_covered_percent",
    "human_overlap_nonhuman_bp",
    "human_overlap_nonhuman_percent",
    "human_outside_nonhuman_bp",
    "largest_nonhuman_block_bp",
    "best_nonhuman_identity_percent",
    "nonhuman_groups",
    "nonhuman_hit_count",
    "reason",
]


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--assembly-id", required=True)
    parser.add_argument("--human-paf", action="append", required=True)
    parser.add_argument(
        "--blast",
        action="append",
        default=[],
        metavar="GROUP=PATH",
        help="Non-human group and BLAST TSV path; may be repeated",
    )
    parser.add_argument("--config", required=True)
    parser.add_argument("--decisions", required=True)
    parser.add_argument("--split-bed", required=True)
    parser.add_argument("--remove-list", required=True)
    parser.add_argument("--review-list", required=True)
    return parser.parse_args()


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fasta_lengths(path):
    lengths = {}
    current = None
    length = 0
    with open_text(path) as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                current = line[1:].strip().split()[0]
                if not current:
                    raise ValueError(f"Empty FASTA identifier at {path}:{line_number}")
                if current in lengths:
                    raise ValueError(f"Duplicate FASTA identifier: {current}")
                length = 0
            else:
                if current is None:
                    raise ValueError(f"Sequence before first FASTA header at {path}:{line_number}")
                length += len(line.strip())
    if current is not None:
        lengths[current] = length
    if not lengths:
        raise ValueError(f"No sequences found in {path}")
    return lengths


def merge_intervals(intervals, gap=0):
    if not intervals:
        return []
    merged = []
    start, end = sorted(intervals)[0]
    for next_start, next_end in sorted(intervals)[1:]:
        if next_start <= end + gap:
            end = max(end, next_end)
        else:
            merged.append((start, end))
            start, end = next_start, next_end
    merged.append((start, end))
    return merged


def interval_length(intervals):
    return sum(end - start for start, end in intervals)


def intersection_length(first, second):
    total = 0
    i = 0
    j = 0
    while i < len(first) and j < len(second):
        start = max(first[i][0], second[j][0])
        end = min(first[i][1], second[j][1])
        if start < end:
            total += end - start
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return total


def read_human_alignments(paths, min_mapq, min_identity):
    intervals = defaultdict(list)
    for path in paths:
        with open_text(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 12:
                    raise ValueError(f"Malformed PAF line at {path}:{line_number}")
                matches = int(fields[9])
                block_length = int(fields[10])
                mapq = int(fields[11])
                identity = 100.0 * matches / block_length if block_length else 0.0
                if mapq >= min_mapq and identity >= min_identity:
                    intervals[fields[0]].append((int(fields[2]), int(fields[3])))
    return intervals


def read_blast_alignments(specifications, min_identity, min_length, max_evalue):
    intervals = defaultdict(list)
    groups = defaultdict(set)
    best_identity = defaultdict(float)
    hit_count = defaultdict(int)

    for specification in specifications:
        if "=" not in specification:
            raise ValueError(f"Expected GROUP=PATH for --blast, got: {specification}")
        group, path = specification.split("=", 1)
        with open_text(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 9:
                    raise ValueError(f"Expected 9 BLAST columns at {path}:{line_number}")
                query = fields[0]
                identity = float(fields[3])
                alignment_length = int(fields[4])
                evalue = float(fields[7])
                if identity < min_identity or alignment_length < min_length or evalue > max_evalue:
                    continue
                qstart = int(fields[5])
                qend = int(fields[6])
                start = min(qstart, qend) - 1
                end = max(qstart, qend)
                intervals[query].append((start, end))
                groups[query].add(group)
                best_identity[query] = max(best_identity[query], identity)
                hit_count[query] += 1
    return intervals, groups, best_identity, hit_count


def classify(length, nonhuman, human, values):
    nonhuman_bp = interval_length(nonhuman)
    human_bp = interval_length(human)
    overlap_bp = intersection_length(nonhuman, human)
    human_outside_bp = max(0, human_bp - overlap_bp)
    nonhuman_percent = 100.0 * nonhuman_bp / length if length else 0.0
    human_percent = 100.0 * human_bp / length if length else 0.0
    overlap_percent = 100.0 * overlap_bp / nonhuman_bp if nonhuman_bp else 0.0
    largest_block = max((end - start for start, end in nonhuman), default=0)

    remove = (
        nonhuman_bp >= values["remove_min_nonhuman_bp"]
        and nonhuman_percent >= values["remove_min_nonhuman_percent"]
        and human_percent <= values["remove_max_human_percent"]
    )
    split = (
        nonhuman_bp >= values["split_min_nonhuman_bp"]
        and nonhuman_percent >= values["split_min_nonhuman_percent"]
        and nonhuman_percent < values["remove_min_nonhuman_percent"]
        and overlap_percent <= values["split_max_human_overlap_percent"]
        and human_outside_bp >= values["split_min_human_outside_bp"]
    )

    if remove:
        decision = "REMOVE"
        reason = "nearly whole contig is non-human with little human-reference support"
    elif split:
        decision = "SPLIT"
        reason = "localized non-human block is separable from human-supported sequence"
    elif nonhuman_bp >= values["review_min_nonhuman_bp"]:
        decision = "REVIEW"
        reason = "non-human evidence is partial, short, or conflicts with human support"
    else:
        decision = "KEEP"
        reason = "no non-human match above the review threshold"

    metrics = {
        "decision": decision,
        "nonhuman_covered_bp": nonhuman_bp,
        "nonhuman_covered_percent": nonhuman_percent,
        "human_covered_bp": human_bp,
        "human_covered_percent": human_percent,
        "human_overlap_nonhuman_bp": overlap_bp,
        "human_overlap_nonhuman_percent": overlap_percent,
        "human_outside_nonhuman_bp": human_outside_bp,
        "largest_nonhuman_block_bp": largest_block,
        "reason": reason,
    }
    return metrics


def main():
    args = parse_args()
    with open(args.config) as handle:
        config = yaml.safe_load(handle)
    blast_config = config["blast"]
    values = config["classification"]

    lengths = fasta_lengths(args.assembly)
    human_raw = read_human_alignments(
        args.human_paf,
        int(values["human_min_mapq"]),
        float(values["human_min_identity"]),
    )
    nonhuman_raw, groups, best_identity, hit_count = read_blast_alignments(
        args.blast,
        float(blast_config["min_identity"]),
        int(blast_config["min_hsp_length_bp"]),
        float(blast_config["evalue"]),
    )
    unknown_queries = (set(human_raw) | set(nonhuman_raw)) - set(lengths)
    if unknown_queries:
        raise ValueError(f"Alignment query IDs absent from FASTA: {sorted(unknown_queries)[:5]}")

    gap = int(values["interval_merge_gap_bp"])
    rows = []
    split_rows = []
    remove_contigs = []
    review_contigs = []

    for contig, length in lengths.items():
        human = merge_intervals(human_raw.get(contig, []))
        nonhuman = merge_intervals(nonhuman_raw.get(contig, []), gap=gap)
        metrics = classify(length, nonhuman, human, values)
        row = {
            "assembly_id": args.assembly_id,
            "contig": contig,
            "contig_length_bp": length,
            **metrics,
            "best_nonhuman_identity_percent": best_identity.get(contig, 0.0),
            "nonhuman_groups": ",".join(sorted(groups.get(contig, []))),
            "nonhuman_hit_count": hit_count.get(contig, 0),
        }
        rows.append(row)
        if metrics["decision"] == "REMOVE":
            remove_contigs.append(contig)
        elif metrics["decision"] == "REVIEW":
            review_contigs.append(contig)
        elif metrics["decision"] == "SPLIT":
            for start, end in nonhuman:
                split_rows.append((contig, start, end))

    for output_path in (args.decisions, args.split_bed, args.remove_list, args.review_list):
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)

    with open(args.decisions, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=DECISION_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with open(args.split_bed, "w") as handle:
        for contig, start, end in split_rows:
            handle.write(f"{contig}\t{start}\t{end}\n")
    with open(args.remove_list, "w") as handle:
        for contig in remove_contigs:
            handle.write(contig + "\n")
    with open(args.review_list, "w") as handle:
        for contig in review_contigs:
            handle.write(contig + "\n")

    counts = defaultdict(int)
    for row in rows:
        counts[row["decision"]] += 1
    print(
        f"{args.assembly_id}: KEEP={counts['KEEP']} REVIEW={counts['REVIEW']} "
        f"SPLIT={counts['SPLIT']} REMOVE={counts['REMOVE']}"
    )


if __name__ == "__main__":
    main()
