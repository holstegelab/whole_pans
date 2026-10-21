#!/usr/bin/env python3
"""Summarize long-read PAF/GAF mappings and calculate paired contrasts."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


METRIC_FIELDS = [
    "sample_id", "sample_role", "graph_id", "graph_membership",
    "membership_status", "mapper", "data_mode",
    "total_reads", "total_bases", "mapped_reads", "mapped_read_fraction",
    "aligned_query_bases", "aligned_base_fraction", "primary_alignments",
    "primary_identity", "clipped_bases", "clipped_base_fraction",
    "split_read_count", "split_read_fraction", "multimapped_read_count",
    "secondary_multimapping_rate", "alignment_length_p10",
    "alignment_length_median", "alignment_length_p90", "mapq_p10",
    "mapq_median", "mapq_p90", "runtime_seconds", "peak_rss_mb",
    "selection_seed", "input_sha256", "selected_bases",
]
NUMERIC_CONTRASTS = [
    "mapped_read_fraction", "aligned_base_fraction", "primary_identity",
    "clipped_base_fraction", "split_read_fraction",
    "secondary_multimapping_rate", "alignment_length_median", "mapq_median",
    "runtime_seconds", "peak_rss_mb",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    summary = sub.add_parser("summarize")
    summary.add_argument("--alignment", required=True, type=Path)
    summary.add_argument("--selection-stats", required=True, type=Path)
    summary.add_argument("--samples", required=True, type=Path)
    summary.add_argument("--membership", required=True, type=Path)
    summary.add_argument("--benchmark", required=True, type=Path)
    summary.add_argument("--sample", required=True)
    summary.add_argument("--graph", required=True)
    summary.add_argument("--mapper", required=True)
    summary.add_argument("--data-mode", default="real")
    summary.add_argument("--output", required=True, type=Path)
    merge = sub.add_parser("merge")
    merge.add_argument("--input", action="append", default=[])
    merge.add_argument("--output", required=True, type=Path)
    contrasts = sub.add_parser("contrasts")
    contrasts.add_argument("--metrics", required=True, type=Path)
    contrasts.add_argument("--output", required=True, type=Path)
    accuracy = sub.add_parser("sim-accuracy")
    accuracy.add_argument("--alignment", required=True, type=Path)
    accuracy.add_argument("--sample", required=True)
    accuracy.add_argument("--graph", required=True)
    accuracy.add_argument("--tolerance", type=int, default=10000)
    accuracy.add_argument("--output", required=True, type=Path)
    merge_accuracy = sub.add_parser("merge-accuracy")
    merge_accuracy.add_argument("--input", action="append", default=[])
    merge_accuracy.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return math.nan
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low = int(math.floor(position))
    high = int(math.ceil(position))
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def merged_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    total = 0
    start, end = sorted(intervals)[0]
    for next_start, next_end in sorted(intervals)[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def tags(fields: list[str]) -> dict[str, str]:
    result = {}
    for tag in fields[12:]:
        parts = tag.split(":", 2)
        if len(parts) == 3:
            result[parts[0]] = parts[2]
    return result


def benchmark_values(path: Path) -> tuple[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return "", ""
    rows = read_table(path)
    if not rows:
        return "", ""
    row = rows[-1]
    seconds = row.get("s", "")
    rss = row.get("max_rss", row.get("max_rss_mb", ""))
    return seconds, rss


def summarize(args: argparse.Namespace) -> dict[str, object]:
    selection_rows = read_table(args.selection_stats)
    if len(selection_rows) != 1:
        raise ValueError(f"expected one selection row in {args.selection_stats}")
    selection = selection_rows[0]
    sample_rows = {row["sample_id"]: row for row in read_table(args.samples)}
    if args.sample not in sample_rows:
        raise ValueError(f"sample {args.sample} absent from {args.samples}")
    membership_rows = [
        row for row in read_table(args.membership)
        if row["sample_id"] == args.sample and row["graph_id"] == args.graph
    ]
    if len(membership_rows) != 1:
        raise ValueError(
            f"expected one membership row for {args.sample}/{args.graph} in {args.membership}"
        )
    membership = membership_rows[0]

    query_intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    nonsecondary_counts: defaultdict[str, int] = defaultdict(int)
    secondary_reads = set()
    mapped_reads = set()
    lengths: list[float] = []
    mapqs: list[float] = []
    primary_matches = 0
    primary_blocks = 0
    primary_alignments = 0
    clipped_bases = 0
    with args.alignment.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"{args.alignment}:{line_number}: expected >=12 fields")
            read = fields[0]
            query_length = int(fields[1])
            query_start, query_end = int(fields[2]), int(fields[3])
            block_length = int(fields[10])
            mapq = int(fields[11])
            record_tags = tags(fields)
            secondary = record_tags.get("tp") == "S"
            mapped_reads.add(read)
            mapqs.append(mapq)
            lengths.append(block_length)
            if secondary:
                secondary_reads.add(read)
                continue
            nonsecondary_counts[read] += 1
            query_intervals[read].append((query_start, query_end))
            primary_alignments += 1
            primary_matches += int(fields[9])
            primary_blocks += block_length
            clipped_bases += query_start + query_length - query_end

    total_reads = int(selection["selected_reads"])
    total_bases = int(selection["selected_bases"])
    aligned_bases = sum(merged_length(value) for value in query_intervals.values())
    split_reads = sum(value > 1 for value in nonsecondary_counts.values())
    runtime, peak_rss = benchmark_values(args.benchmark)
    def ratio(numerator: float, denominator: float) -> float:
        return numerator / denominator if denominator else math.nan

    return {
        "sample_id": args.sample,
        "sample_role": sample_rows[args.sample]["role"],
        "graph_id": args.graph,
        "graph_membership": membership["is_member"],
        "membership_status": membership["membership_status"],
        "mapper": args.mapper,
        "data_mode": args.data_mode,
        "total_reads": total_reads,
        "total_bases": total_bases,
        "mapped_reads": len(mapped_reads),
        "mapped_read_fraction": ratio(len(mapped_reads), total_reads),
        "aligned_query_bases": aligned_bases,
        "aligned_base_fraction": ratio(aligned_bases, total_bases),
        "primary_alignments": primary_alignments,
        "primary_identity": ratio(primary_matches, primary_blocks),
        "clipped_bases": clipped_bases,
        "clipped_base_fraction": ratio(clipped_bases, total_bases),
        "split_read_count": split_reads,
        "split_read_fraction": ratio(split_reads, total_reads),
        "multimapped_read_count": len(secondary_reads),
        "secondary_multimapping_rate": ratio(len(secondary_reads), total_reads),
        "alignment_length_p10": percentile(lengths, 0.10),
        "alignment_length_median": percentile(lengths, 0.50),
        "alignment_length_p90": percentile(lengths, 0.90),
        "mapq_p10": percentile(mapqs, 0.10),
        "mapq_median": percentile(mapqs, 0.50),
        "mapq_p90": percentile(mapqs, 0.90),
        "runtime_seconds": runtime,
        "peak_rss_mb": peak_rss,
        "selection_seed": selection.get("selection_seed", ""),
        "input_sha256": selection.get("input_sha256", ""),
        "selected_bases": selection.get("selected_bases", ""),
    }


def write_rows(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def merge(args: argparse.Namespace) -> None:
    rows = []
    for path_text in args.input:
        rows.extend(read_table(Path(path_text)))
    rows.sort(key=lambda row: (row.get("data_mode", ""), row.get("sample_role", ""), row.get("sample_id", ""), row.get("mapper", ""), row.get("graph_id", "")))
    write_rows(args.output, rows, METRIC_FIELDS)


def contrasts(args: argparse.Namespace) -> None:
    rows = read_table(args.metrics)
    index = {
        (row["sample_id"], row["mapper"], row["graph_id"]): row
        for row in rows if row.get("data_mode") == "real"
    }
    pairs = [
        ("HPRC_minus_linear", "chm13_linear", "hprc_chm13_sv"),
        ("local_V1_minus_HPRC", "hprc_chm13_sv", "local_v1"),
        ("local_V2_minus_local_V1", "local_v1", "local_v2"),
        ("HPRC_plus_local_minus_HPRC", "hprc_chm13_sv", "hprc_plus_local"),
    ]
    output = []

    def membership_transition(left: dict | None, right: dict | None) -> str:
        if left is None or right is None:
            return "unknown_missing_mapping"
        left_status = left.get("graph_membership", "unknown")
        right_status = right.get("graph_membership", "unknown")
        labels = {"false": "nonmember", "true": "member"}
        if left_status not in labels or right_status not in labels:
            return "unknown_membership"
        return f"{labels[left_status]}_to_{labels[right_status]}"

    samples = sorted({row["sample_id"] for row in rows})
    for sample in samples:
        role = next(row["sample_role"] for row in rows if row["sample_id"] == sample)
        for contrast_id, left_graph, right_graph in pairs:
            left = index.get((sample, "minigraph", left_graph))
            right = index.get((sample, "minigraph", right_graph))
            for metric in NUMERIC_CONTRASTS:
                status = "PASS"
                value = ""
                if left is None or right is None:
                    status = "BLOCKED_missing_paired_mapping"
                else:
                    try:
                        value = float(right[metric]) - float(left[metric])
                    except (KeyError, TypeError, ValueError):
                        status = "BLOCKED_metric_unavailable"
                output.append({
                    "sample_id": sample,
                    "sample_role": role,
                    "headline_eligible": "true" if role == "held_out_primary" else "false",
                    "contrast_id": contrast_id,
                    "left_graph": left_graph,
                    "right_graph": right_graph,
                    "left_graph_membership": left.get("graph_membership", "") if left else "",
                    "right_graph_membership": right.get("graph_membership", "") if right else "",
                    "membership_transition": membership_transition(left, right),
                    "mapper": "minigraph",
                    "metric": metric,
                    "delta_right_minus_left": value,
                    "status": status,
                })
    write_rows(args.output, output, [
        "sample_id", "sample_role", "headline_eligible", "contrast_id",
        "left_graph", "right_graph", "left_graph_membership",
        "right_graph_membership", "membership_transition", "mapper", "metric",
        "delta_right_minus_left", "status",
    ])


SIM_NAME = re.compile(
    r"^sim\|(?P<sample>[^|]+)\|(?P<hap>[^|]+)\|(?P<contig>[^:|]+):(?P<start>\d+)-(?P<end>\d+)\|"
)


def target_coordinates(fields: list[str]) -> tuple[str, int, int]:
    path = fields[5]
    start, end = int(fields[7]), int(fields[8])
    if not path.startswith((">", "<")):
        return path, start, end
    match = re.search(r"(?:CHM13[.#])?(chr(?:[0-9]+|X|Y|M)):(\d+)-(\d+)", path)
    if match:
        return match.group(1), int(match.group(2)), int(match.group(3))
    chrom = re.search(r"chr(?:[0-9]+|X|Y|M)", path)
    return (chrom.group(0) if chrom else path, start, end)


def sim_accuracy(args: argparse.Namespace) -> None:
    best: dict[str, tuple[tuple[int, int], list[str]]] = {}
    with args.alignment.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or tags(fields).get("tp") == "S":
                continue
            score = (int(fields[11]), int(fields[10]))
            if fields[0] not in best or score > best[fields[0]][0]:
                best[fields[0]] = (score, fields)
    total = correct_chrom = correct_locus = wrong_locus = 0
    identity_sum = aligned_sum = clipped_sum = 0.0
    hap_assessable = hap_correct = 0
    for name, (_, fields) in best.items():
        match = SIM_NAME.match(name)
        if not match:
            continue
        total += 1
        source_chrom = match.group("contig")
        source_start, source_end = int(match.group("start")), int(match.group("end"))
        target_chrom, target_start, target_end = target_coordinates(fields)
        same_chrom = source_chrom == target_chrom or source_chrom.removeprefix("chr") == target_chrom.removeprefix("chr")
        if same_chrom:
            correct_chrom += 1
        overlap = target_end + args.tolerance >= source_start and target_start - args.tolerance <= source_end
        if same_chrom and overlap:
            correct_locus += 1
        else:
            wrong_locus += 1
        identity_sum += int(fields[9]) / int(fields[10]) if int(fields[10]) else 0
        aligned_sum += (int(fields[3]) - int(fields[2])) / int(fields[1])
        clipped_sum += (int(fields[2]) + int(fields[1]) - int(fields[3])) / int(fields[1])
        if match.group("sample") in fields[5] or match.group("hap") in fields[5]:
            hap_assessable += 1
            if match.group("sample") in fields[5] and match.group("hap") in fields[5]:
                hap_correct += 1
    denominator = total or 1
    row = {
        "sample_id": args.sample,
        "graph_id": args.graph,
        "simulated_reads_assessed": total,
        "correct_chromosome_fraction": correct_chrom / denominator,
        "correct_locus_fraction": correct_locus / denominator,
        "wrong_locus_fraction": wrong_locus / denominator,
        "correct_haplotype_fraction": hap_correct / hap_assessable if hap_assessable else "NA_not_represented_by_minigraph_paths",
        "haplotype_assessable_reads": hap_assessable,
        "mean_identity": identity_sum / denominator,
        "mean_aligned_fraction": aligned_sum / denominator,
        "mean_clipped_fraction": clipped_sum / denominator,
        "truth_basis": "read_name_source_haplotype_and_interval",
    }
    write_rows(args.output, [row], list(row))


def merge_accuracy(args: argparse.Namespace) -> None:
    fields = [
        "sample_id", "graph_id", "simulated_reads_assessed",
        "correct_chromosome_fraction", "correct_locus_fraction",
        "wrong_locus_fraction", "correct_haplotype_fraction",
        "haplotype_assessable_reads", "mean_identity", "mean_aligned_fraction",
        "mean_clipped_fraction", "truth_basis",
    ]
    rows = [row for path in args.input for row in read_table(Path(path))]
    rows.sort(key=lambda row: (row.get("sample_id", ""), row.get("graph_id", "")))
    write_rows(args.output, rows, fields)


def main() -> None:
    args = parse_args()
    if args.command == "summarize":
        write_rows(args.output, [summarize(args)], METRIC_FIELDS)
    elif args.command == "merge":
        merge(args)
    elif args.command == "contrasts":
        contrasts(args)
    elif args.command == "sim-accuracy":
        sim_accuracy(args)
    elif args.command == "merge-accuracy":
        merge_accuracy(args)


if __name__ == "__main__":
    main()
