#!/usr/bin/env python3
"""Stream Giraffe GAF into genome-wide QC and graph-specific TMEM evidence."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
import sys
from collections import Counter
from pathlib import Path
from typing import Iterable, TextIO


EVIDENCE_FIELDS = [
    "sample_id", "sample_role", "batch", "graph_id", "baseline_graph",
    "indexing_mode", "run_mode", "target_present", "target_callable_status",
    "evidence_status", "min_mapq", "min_junction_reads_per_side",
    "left_alt_junction_distinct_primary_reads",
    "right_alt_junction_distinct_primary_reads",
    "reference_traversal_distinct_primary_reads",
    "insertion_node_or_edge_distinct_primary_reads", "minimum_junction_support",
    "left_junction_support_per_million_mapped_reads",
    "right_junction_support_per_million_mapped_reads",
    "minimum_junction_support_per_million_mapped_reads", "local_allele_balance",
    "both_junctions_meet_threshold", "callability", "ambiguity",
    "mapped_read_denominator", "distinct_read_name_definition",
    "target_failure_reason", "limitations", "warnings",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--role", required=True)
    parser.add_argument("--batch", default="")
    parser.add_argument("--graph", required=True)
    parser.add_argument("--baseline-graph", required=True)
    parser.add_argument("--indexing-mode", required=True)
    parser.add_argument("--run-mode", choices=("smoke", "full"), required=True)
    parser.add_argument("--fastq-stats", required=True, type=Path)
    parser.add_argument("--target", required=True, type=Path)
    parser.add_argument("--min-mapq", type=int, default=5)
    parser.add_argument("--min-junction", type=int, default=3)
    parser.add_argument("--metrics-output", required=True, type=Path)
    parser.add_argument("--evidence-output", required=True, type=Path)
    parser.add_argument("--debug-gaf", type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, row: dict, fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerow(row)
    temporary.replace(path)


def parse_tags(fields: list[str]) -> dict[str, str]:
    result = {}
    for value in fields[12:]:
        parts = value.split(":", 2)
        if len(parts) == 3:
            result[parts[0]] = parts[2]
    return result


def path_tokens(value: str) -> list[str]:
    if not value.startswith((">", "<")):
        return []
    return [orientation + name for orientation, name in re.findall(r"([<>])([^<>]+)", value)]


def flip(token: str) -> str:
    return ("<" if token.startswith(">") else ">") + token[1:]


def reverse_edge(edge: tuple[str, str]) -> tuple[str, str]:
    return flip(edge[1]), flip(edge[0])


def parse_edges(value: str) -> set[tuple[str, str]]:
    result = set()
    for item in value.split(";"):
        parts = [part.strip() for part in item.split(",") if part.strip()]
        if len(parts) == 2 and all(part.startswith((">", "<")) for part in parts):
            edge = (parts[0], parts[1])
            result.add(edge)
            result.add(reverse_edge(edge))
    return result


def normalized_name(value: str) -> str:
    return re.sub(r"/(?:1|2)$", "", value.split()[0])


def histogram_median(histogram: Counter[int], count: int) -> float | str:
    if count == 0:
        return ""
    lower_rank = (count - 1) // 2
    upper_rank = count // 2
    seen = 0
    lower = upper = None
    for value in sorted(histogram):
        next_seen = seen + histogram[value]
        if lower is None and lower_rank < next_seen:
            lower = value
        if upper_rank < next_seen:
            upper = value
            break
        seen = next_seen
    return (float(lower) + float(upper)) / 2


def ratio(numerator: float, denominator: float) -> float | str:
    return numerator / denominator if denominator else ""


def target_definition(row: dict[str, str]) -> dict[str, object]:
    return {
        "present": row.get("target_present", "false"),
        "callable": row.get("callable_status", "not_callable_target_ambiguous"),
        "left": parse_edges(row.get("left_junction", "")),
        "right": parse_edges(row.get("right_junction", "")),
        "reference": parse_edges(row.get("reference_defining_edges", "")),
        "insertion_edges": parse_edges(row.get("insertion_edges", "")),
        "insertion_nodes": {
            value.lstrip("><") for value in row.get("insertion_nodes", "").split(",") if value
        },
        "failure": row.get("failure_reason", ""),
        "ambiguity": row.get("ambiguity", ""),
    }


def summarize_stream(
    lines: Iterable[str], total_reads: int, metadata: dict[str, object],
    target: dict[str, object], min_mapq: int, min_junction: int,
    debug_handle: TextIO | None = None,
) -> tuple[dict[str, object], dict[str, object]]:
    observed = malformed = gaf_unmapped = primary_mapped = secondary = 0
    mapq_hist: Counter[int] = Counter()
    mapq_sum = identity_sum = aligned_fraction_sum = clipped_sum = 0.0
    left_names: set[str] = set()
    right_names: set[str] = set()
    reference_names: set[str] = set()
    insertion_names: set[str] = set()
    warnings: list[str] = []

    for line_number, line in enumerate(lines, start=1):
        if debug_handle is not None:
            debug_handle.write(line)
        if not line.strip():
            continue
        observed += 1
        fields = line.rstrip("\r\n").split("\t")
        if len(fields) < 12:
            malformed += 1
            continue
        try:
            query_length = int(fields[1])
            query_start = int(fields[2])
            query_end = int(fields[3])
            matches = int(fields[9])
            block_length = int(fields[10])
            mapq = int(fields[11])
        except ValueError:
            malformed += 1
            continue
        if query_length < 0 or query_start < 0 or query_end < query_start or query_end > query_length:
            malformed += 1
            continue
        tags = parse_tags(fields)
        is_secondary = tags.get("tp") == "S"
        mapped = fields[5] != "*" and block_length > 0 and query_end > query_start
        if not mapped:
            gaf_unmapped += 1
            continue
        if is_secondary:
            secondary += 1
            continue

        primary_mapped += 1
        mapq_hist[mapq] += 1
        mapq_sum += mapq
        identity_sum += matches / block_length
        aligned_fraction_sum += (query_end - query_start) / query_length if query_length else 0
        clipped_sum += query_start + query_length - query_end

        if mapq < min_mapq or not str(target["callable"]).startswith("callable"):
            continue
        tokens = path_tokens(fields[5])
        if not tokens:
            continue
        edges = set(zip(tokens, tokens[1:]))
        names = {token[1:] for token in tokens}
        read_name = normalized_name(fields[0])
        if edges & target["left"]:
            left_names.add(read_name)
        if edges & target["right"]:
            right_names.add(read_name)
        if edges & target["reference"]:
            reference_names.add(read_name)
        if names & target["insertion_nodes"] or edges & target["insertion_edges"]:
            insertion_names.add(read_name)

    if malformed:
        warnings.append(f"malformed_gaf_rows={malformed}")
    if primary_mapped > total_reads:
        warnings.append("primary_mapped_records_exceed_independent_FASTQ_read_denominator")
    unmapped = max(total_reads - primary_mapped, 0)
    mapq_denominator = primary_mapped
    metrics = {
        **metadata,
        "total_input_reads": total_reads,
        "total_input_pairs": total_reads // 2,
        "gaf_records_observed": observed,
        "gaf_unmapped_records": gaf_unmapped,
        "primary_mapped_reads": primary_mapped,
        "unmapped_reads": unmapped,
        "mapped_read_fraction": ratio(primary_mapped, total_reads),
        "mapq_fraction_denominator": "primary_mapped_reads",
        "mapq_0_fraction": ratio(mapq_hist[0], mapq_denominator),
        "mapq_ge_5_fraction": ratio(sum(count for value, count in mapq_hist.items() if value >= 5), mapq_denominator),
        "mapq_ge_20_fraction": ratio(sum(count for value, count in mapq_hist.items() if value >= 20), mapq_denominator),
        "mapq_ge_40_fraction": ratio(sum(count for value, count in mapq_hist.items() if value >= 40), mapq_denominator),
        "mean_mapq": ratio(mapq_sum, mapq_denominator),
        "median_mapq": histogram_median(mapq_hist, mapq_denominator),
        "mean_alignment_identity": ratio(identity_sum, primary_mapped),
        "mean_query_aligned_fraction": ratio(aligned_fraction_sum, primary_mapped),
        "mean_clipped_bases_per_read": ratio(clipped_sum, primary_mapped),
        "secondary_multimapping_count": secondary,
        "secondary_multimapping_fraction": ratio(secondary, total_reads),
        "malformed_gaf_rows": malformed,
        "relevant_warnings": ";".join(warnings),
    }

    callable_status = str(target["callable"])
    if callable_status == "not_callable_graph_target_absent":
        evidence_status = "not_callable_graph_target_absent"
    elif not callable_status.startswith("callable"):
        evidence_status = "not_callable_target_ambiguous"
    elif metadata["run_mode"] == "smoke":
        evidence_status = "QC_only_smoke_depth_not_genotype_interpretable"
    else:
        evidence_status = "callable_evidence_reported"
    callable_target = callable_status.startswith("callable")
    left_count = len(left_names) if callable_target else ""
    right_count = len(right_names) if callable_target else ""
    ref_count = len(reference_names) if callable_target else ""
    insertion_count = len(insertion_names) if callable_target else ""
    minimum = min(left_count, right_count) if callable_target else ""
    million = primary_mapped / 1_000_000 if primary_mapped else 0
    balance_denominator = minimum + ref_count if callable_target else 0
    evidence = {
        **{key: metadata[key] for key in (
            "sample_id", "sample_role", "batch", "graph_id", "baseline_graph",
            "indexing_mode", "run_mode",
        )},
        "target_present": target["present"],
        "target_callable_status": callable_status,
        "evidence_status": evidence_status,
        "min_mapq": min_mapq,
        "min_junction_reads_per_side": min_junction,
        "left_alt_junction_distinct_primary_reads": left_count,
        "right_alt_junction_distinct_primary_reads": right_count,
        "reference_traversal_distinct_primary_reads": ref_count,
        "insertion_node_or_edge_distinct_primary_reads": insertion_count,
        "minimum_junction_support": minimum,
        "left_junction_support_per_million_mapped_reads": ratio(left_count, million) if callable_target else "",
        "right_junction_support_per_million_mapped_reads": ratio(right_count, million) if callable_target else "",
        "minimum_junction_support_per_million_mapped_reads": ratio(minimum, million) if callable_target else "",
        "local_allele_balance": ratio(minimum, balance_denominator) if callable_target else "",
        "both_junctions_meet_threshold": str(
            callable_target and left_count >= min_junction and right_count >= min_junction
        ).lower() if callable_target else "",
        "callability": callable_status,
        "ambiguity": target["ambiguity"],
        "mapped_read_denominator": primary_mapped,
        "distinct_read_name_definition": "normalized_GAF_query_name_with_/1_/2_suffix_removed",
        "target_failure_reason": target["failure"],
        "limitations": (
            "TMEM-positive_candidates_selected_from_prior_local-graph_signal_not_independent_truth;"
            "random_samples_blinded_unknown;internal_repetitive_sequence_alone_not_carrier_evidence;"
            + ("smoke_depth_absence_not_genotype_result" if metadata["run_mode"] == "smoke" else "real_reads_no_independent_truth")
        ),
        "warnings": ";".join(warnings),
    }
    return metrics, evidence


def main() -> None:
    args = parse_args()
    stats_rows = read_tsv(args.fastq_stats)
    target_rows = read_tsv(args.target)
    if len(stats_rows) != 1:
        raise ValueError(f"expected one FASTQ-stat row: {args.fastq_stats}")
    if len(target_rows) != 1 or target_rows[0].get("graph_id") != args.graph:
        raise ValueError(f"expected one target-audit row for {args.graph}")
    metadata = {
        "sample_id": args.sample,
        "sample_role": args.role,
        "batch": args.batch,
        "graph_id": args.graph,
        "baseline_graph": args.baseline_graph,
        "mapper": "vg_giraffe",
        "indexing_mode": args.indexing_mode,
        "run_mode": args.run_mode,
    }
    debug = None
    if args.debug_gaf:
        args.debug_gaf.parent.mkdir(parents=True, exist_ok=True)
        debug = gzip.open(args.debug_gaf, "wt")
    try:
        metrics, evidence = summarize_stream(
            sys.stdin,
            int(stats_rows[0]["total_input_reads"]),
            metadata,
            target_definition(target_rows[0]),
            args.min_mapq,
            args.min_junction,
            debug,
        )
    finally:
        if debug is not None:
            debug.close()
    metric_fields = list(metrics)
    write_tsv(args.metrics_output, metrics, metric_fields)
    write_tsv(args.evidence_output, evidence, EVIDENCE_FIELDS)


if __name__ == "__main__":
    main()
