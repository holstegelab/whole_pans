#!/usr/bin/env python3
"""Assess targeted V1/V2 support for configured PASS rescue events."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path

import yaml


PATH_NODE = re.compile(r"[<>]([^<>]+)")
CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--v1-graph", required=True, type=Path)
    parser.add_argument("--v2-graph", required=True, type=Path)
    parser.add_argument("--alignment", action="append", default=[], help="graph_id=sample_id=GAF")
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def read_table(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def graph_segments(path: Path) -> tuple[dict[str, int], dict[str, str]]:
    lengths = {}
    aliases = {}
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith("S\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            name, sequence = fields[1], fields[2]
            length = len(sequence) if sequence != "*" else 0
            for tag in fields[3:]:
                if tag.startswith("LN:i:"):
                    length = int(tag[5:])
                elif tag.startswith("SN:Z:"):
                    aliases[name] = tag[5:]
            lengths[name] = length
    return lengths, aliases


def tag_values(fields: list[str]) -> dict[str, str]:
    result = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            result[parts[0]] = parts[2]
    return result


def path_event_coordinate(
    path: str,
    source_segment: str,
    offset: int,
    lengths: dict[str, int],
    aliases: dict[str, str],
) -> int | None:
    if not path.startswith((">", "<")):
        return offset if path == source_segment else None
    cumulative = 0
    for node in PATH_NODE.findall(path):
        alias = aliases.get(node, node)
        length = lengths.get(node, 0)
        if node == source_segment or alias == source_segment:
            return cumulative + min(offset, length)
        cumulative += length
    return None


def alignment_evidence(
    path: Path,
    events: dict[str, dict],
    lengths: dict[str, int],
    aliases: dict[str, str],
    anchor_bp: int,
    novel_nodes: set[str] | None = None,
) -> dict[str, dict[str, object]]:
    evidence = defaultdict(lambda: {
        "mapped_records": 0,
        "both_anchor_records": 0,
        "alternate_traversal_records": 0,
        "secondary_records": 0,
        "max_indel_bp": [],
        "clipped_bp": [],
    })
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            record_tags = tag_values(fields)
            target_start, target_end = int(fields[7]), int(fields[8])
            query_length = int(fields[1])
            clipped = int(fields[2]) + query_length - int(fields[3])
            cigar_indels = [
                int(length) for length, operation in CIGAR.findall(record_tags.get("cg", ""))
                if operation in {"I", "D"}
            ]
            for event_id, event in events.items():
                coordinate = path_event_coordinate(
                    fields[5], event.get("graph_segment", ""),
                    int(event.get("segment_offset_0") or 0), lengths, aliases,
                )
                path_node_names = PATH_NODE.findall(fields[5])
                alt = event_id in fields[5] or any(
                    event_id in aliases.get(node, "") for node in path_node_names
                )
                if coordinate is not None and novel_nodes:
                    alt = alt or any(node in novel_nodes for node in path_node_names)
                if coordinate is None and not alt:
                    continue
                current = evidence[event_id]
                current["mapped_records"] += 1
                if record_tags.get("tp") == "S":
                    current["secondary_records"] += 1
                    continue
                if coordinate is not None and target_start <= coordinate - anchor_bp and target_end >= coordinate + anchor_bp:
                    current["both_anchor_records"] += 1
                if alt:
                    current["alternate_traversal_records"] += 1
                current["max_indel_bp"].append(max(cigar_indels, default=0))
                current["clipped_bp"].append(clipped)
    return evidence


def median(values: list[int]) -> float | str:
    if not values:
        return ""
    values = sorted(values)
    middle = len(values) // 2
    return values[middle] if len(values) % 2 else (values[middle - 1] + values[middle]) / 2


def main() -> None:
    args = parse_args()
    with args.config.open() as handle:
        config = yaml.safe_load(handle)
    rescue_config = config["v2_rescue_validation"]
    summary = {row["rescue_id"]: row for row in read_table(Path(rescue_config["rescue_summary"]))}
    catalog = {row["event_id"]: row for row in read_table(Path(rescue_config["normalized_catalog"]))}
    events = {}
    for configured in rescue_config["events"]:
        event_id = str(configured["rescue_id"])
        if event_id not in summary:
            raise ValueError(f"configured rescue absent from rescue_summary.tsv: {event_id}")
        if summary[event_id].get("status") != "PASS":
            raise ValueError(f"configured rescue is not PASS: {event_id} ({summary[event_id].get('status')})")
        if event_id not in catalog:
            raise ValueError(f"configured rescue absent from normalized source catalog: {event_id}")
        if str(summary[event_id].get("svtype")) != str(configured.get("expected_svtype")):
            raise ValueError(f"{event_id}: configured SV type does not match rescue summary")
        if abs(abs(int(summary[event_id]["svlen"])) - abs(int(configured["expected_svlen"]))) > 10:
            raise ValueError(f"{event_id}: configured SV length does not match rescue summary")
        events[event_id] = {**configured, **catalog[event_id], **summary[event_id]}

    v1_id = rescue_config.get("v1_graph_id", "local_v1")
    v2_id = rescue_config.get("v2_graph_id", "local_v2")
    graph_data = {v1_id: graph_segments(args.v1_graph), v2_id: graph_segments(args.v2_graph)}
    v2_novel_nodes = set(graph_data[v2_id][0]) - set(graph_data[v1_id][0])
    observed = {}
    for specification in args.alignment:
        graph_id, sample_id, path_text = specification.split("=", 2)
        lengths, aliases = graph_data[graph_id]
        observed[(graph_id, sample_id)] = alignment_evidence(
            Path(path_text), events, lengths, aliases,
            int(rescue_config.get("anchor_bp", 1000)),
            v2_novel_nodes if graph_id == v2_id else set(),
        )

    rows = []
    max_ambiguous = float(rescue_config.get("max_ambiguous_fraction", 0.5))
    for event_id, event in events.items():
        samples = [(event["sample_id"], "rescue_source")]
        samples += [(sample, "independent_carrier") for sample in event.get("independent_carriers", [])]
        samples += [(sample, "noncarrier_control") for sample in event.get("noncarrier_controls", [])]
        independent_replicated = False
        sample_results = []
        for sample_id, evidence_role in samples:
            v1 = observed.get((v1_id, sample_id), {}).get(event_id, {})
            v2 = observed.get((v2_id, sample_id), {}).get(event_id, {})
            v2_mapped = int(v2.get("mapped_records", 0))
            ambiguous_fraction = int(v2.get("secondary_records", 0)) / v2_mapped if v2_mapped else 0
            captured = int(v2.get("both_anchor_records", 0)) > 0 and int(v2.get("alternate_traversal_records", 0)) > 0
            if captured and evidence_role == "independent_carrier":
                independent_replicated = True
            if ambiguous_fraction > max_ambiguous:
                status = "ambiguous_repetitive"
            elif captured:
                status = "captured_in_V2"
            else:
                status = "not_captured"
            sample_results.append((sample_id, evidence_role, status, v1, v2, ambiguous_fraction))
        for sample_id, evidence_role, status, v1, v2, ambiguous_fraction in sample_results:
            if status == "captured_in_V2":
                interpretation = "independently_replicated" if independent_replicated else "source_only_evidence"
            else:
                interpretation = status
            rows.append({
                "rescue_id": event_id,
                "source_sample_id": event["sample_id"],
                "evidence_sample_id": sample_id,
                "evidence_role": evidence_role,
                "svtype": event["svtype"],
                "svlen": event["svlen"],
                "rescue_summary_status": event["status"],
                "chromosome": event.get("chromosome", ""),
                "position_0": event.get("position_0", ""),
                "v1_both_anchor_records": v1.get("both_anchor_records", 0),
                "v2_both_anchor_records": v2.get("both_anchor_records", 0),
                "v1_alternate_traversal_records": v1.get("alternate_traversal_records", 0),
                "v2_alternate_traversal_records": v2.get("alternate_traversal_records", 0),
                "v1_median_max_indel_bp": median(v1.get("max_indel_bp", [])),
                "v2_median_max_indel_bp": median(v2.get("max_indel_bp", [])),
                "v1_median_clipped_bp": median(v1.get("clipped_bp", [])),
                "v2_median_clipped_bp": median(v2.get("clipped_bp", [])),
                "v2_secondary_fraction": ambiguous_fraction,
                "capture_status": status,
                "replication_status": interpretation,
                "independent_carrier_tested": "true" if event.get("independent_carriers") else "false",
                "v2_provenance_warning": "experimental_graph_overall_cohort_accounting_QC_FAIL",
            })

    fields = list(rows[0]) if rows else ["rescue_id", "capture_status"]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.output) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(args.output)


if __name__ == "__main__":
    main()
