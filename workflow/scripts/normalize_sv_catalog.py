#!/usr/bin/env python3
"""Normalize heterogeneous SV tables and cluster comparable anchored events."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from collections import defaultdict
from difflib import SequenceMatcher
from pathlib import Path

import yaml


FIELDS = [
    "stable_event_id", "chm13_chrom", "chm13_start_0", "chm13_end_0",
    "grch38_chrom", "grch38_start_0", "grch38_end_0", "graph_node",
    "graph_offset_0", "svtype", "svlen", "inserted_sequence",
    "inserted_sequence_sha256", "source_graph", "source_sample",
    "discovery_method", "breakpoint_confidence", "sequence_confidence",
    "graph_membership_status", "projection_status", "source_event_id",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    normalize = sub.add_parser("normalize")
    normalize.add_argument("--config", required=True, type=Path)
    normalize.add_argument("--output", required=True, type=Path)
    cluster = sub.add_parser("cluster")
    cluster.add_argument("--catalog", required=True, type=Path)
    cluster.add_argument("--breakpoint-distance", type=int, required=True)
    cluster.add_argument("--length-similarity", type=float, required=True)
    cluster.add_argument("--sequence-similarity", type=float, required=True)
    cluster.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def stable_id(parts: list[str]) -> str:
    return "SV_" + hashlib.sha256("\0".join(parts).encode()).hexdigest()[:16]


def literal_sequence(value: str) -> str:
    value = (value or "").upper()
    return value if value and not value.startswith("<") and set(value) <= {"A", "C", "G", "T", "N"} else ""


def normalize_provisional(row: dict[str, str], source: dict) -> dict[str, object]:
    coordinate = row.get("coordinate_system", "")
    chrom = row.get("chromosome", "")
    start = row.get("position_0", "")
    end = row.get("end_0", "")
    inserted = literal_sequence(row.get("alternate_allele", "")) if row.get("svtype") == "INS" else ""
    source_id = row.get("event_id", "")
    event_id = source_id or stable_id([coordinate, chrom, start, end, row.get("svtype", ""), row.get("svlen", ""), inserted])
    is_chm13 = coordinate.upper() == "CHM13"
    is_grch38 = coordinate.upper() in {"GRCH38", "HG38"}
    return {
        "stable_event_id": event_id,
        "chm13_chrom": chrom if is_chm13 else "",
        "chm13_start_0": start if is_chm13 else "",
        "chm13_end_0": end if is_chm13 else "",
        "grch38_chrom": chrom if is_grch38 else "",
        "grch38_start_0": start if is_grch38 else "",
        "grch38_end_0": end if is_grch38 else "",
        "graph_node": row.get("graph_segment", ""),
        "graph_offset_0": row.get("segment_offset_0", ""),
        "svtype": row.get("svtype", ""),
        "svlen": row.get("svlen", ""),
        "inserted_sequence": inserted,
        "inserted_sequence_sha256": hashlib.sha256(inserted.encode()).hexdigest() if inserted else "",
        "source_graph": source.get("source_graph", ""),
        "source_sample": row.get("carrier_samples", ""),
        "discovery_method": row.get("discovery_methods", source.get("discovery_method", "")),
        "breakpoint_confidence": row.get("confidence", "unresolved"),
        "sequence_confidence": "sequence_resolved" if inserted else "sequence_unavailable_or_not_insertion",
        "graph_membership_status": row.get("graph_representation_status", "unknown"),
        "projection_status": "projected" if is_chm13 or is_grch38 else "unprojected_preserved",
        "source_event_id": source_id,
    }


def parse_info(value: str) -> dict[str, str]:
    result = {}
    for item in value.split(";"):
        key, _, val = item.partition("=")
        result[key] = val
    return result


def normalize_vcf(fields: list[str], source: dict, ordinal: int) -> list[dict[str, object]]:
    chrom, pos, record_id, ref, alt_field, _, filt, info_text = fields[:8]
    info = parse_info(info_text)
    rows = []
    for alt_index, alt in enumerate(alt_field.split(","), start=1):
        svtype = info.get("SVTYPE", "INS" if len(alt) > len(ref) else "DEL")
        svlen_text = info.get("SVLEN", str(abs(len(alt) - len(ref)))).split(",")[0]
        end_0 = int(info.get("END", int(pos) + max(len(ref), 1)))
        inserted = literal_sequence(alt[len(ref):] if svtype == "INS" and alt.startswith(ref) else alt if svtype == "INS" else "")
        source_id = record_id if record_id != "." else f"record_{ordinal}_ALT{alt_index}"
        coordinate = str(source.get("coordinate_system", "unprojected")).upper()
        event_id = stable_id([coordinate, chrom, str(int(pos) - 1), str(end_0), svtype, svlen_text, inserted])
        rows.append({
            "stable_event_id": event_id,
            "chm13_chrom": chrom if coordinate == "CHM13" else "",
            "chm13_start_0": int(pos) - 1 if coordinate == "CHM13" else "",
            "chm13_end_0": end_0 if coordinate == "CHM13" else "",
            "grch38_chrom": chrom if coordinate in {"GRCH38", "HG38"} else "",
            "grch38_start_0": int(pos) - 1 if coordinate in {"GRCH38", "HG38"} else "",
            "grch38_end_0": end_0 if coordinate in {"GRCH38", "HG38"} else "",
            "graph_node": "", "graph_offset_0": "", "svtype": svtype,
            "svlen": svlen_text, "inserted_sequence": inserted,
            "inserted_sequence_sha256": hashlib.sha256(inserted.encode()).hexdigest() if inserted else "",
            "source_graph": source.get("source_graph", ""),
            "source_sample": source.get("source_sample", ""),
            "discovery_method": source.get("discovery_method", "VCF"),
            "breakpoint_confidence": "PASS" if filt == "PASS" else filt,
            "sequence_confidence": "sequence_resolved" if inserted else "sequence_unavailable_or_symbolic",
            "graph_membership_status": source.get("graph_membership_status", "unknown"),
            "projection_status": "projected" if coordinate in {"CHM13", "GRCH38", "HG38"} else "unprojected_preserved",
            "source_event_id": source_id,
        })
    return rows


def normalize(args: argparse.Namespace) -> None:
    with args.config.open() as handle:
        config = yaml.safe_load(handle)
    inputs = config.get("sv_catalog", {}).get("inputs", [])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.output) + ".tmp")
    seen = set()
    with temporary.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=FIELDS, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for source in inputs:
            path = Path(source["path"])
            if not path.is_file():
                raise FileNotFoundError(f"SV catalog input is missing: {path}")
            fmt = source.get("format", "provisional_catalog")
            with open_text(path) as handle:
                if fmt == "provisional_catalog":
                    for row in csv.DictReader(handle, delimiter="\t"):
                        normalized_rows = [normalize_provisional(row, source)]
                        for normalized in normalized_rows:
                            identifier = str(normalized["stable_event_id"])
                            if identifier in seen:
                                identifier = stable_id([identifier, source.get("source_graph", ""), normalized["source_event_id"]])
                                normalized["stable_event_id"] = identifier
                            seen.add(identifier)
                            writer.writerow(normalized)
                elif fmt == "vcf":
                    ordinal = 0
                    for line in handle:
                        if line.startswith("#"):
                            continue
                        ordinal += 1
                        for normalized in normalize_vcf(line.rstrip("\n").split("\t"), source, ordinal):
                            identifier = str(normalized["stable_event_id"])
                            if identifier in seen:
                                identifier = stable_id([identifier, source.get("source_graph", ""), normalized["source_event_id"]])
                                normalized["stable_event_id"] = identifier
                            seen.add(identifier)
                            writer.writerow(normalized)
                else:
                    raise ValueError(f"unsupported SV input format: {fmt}")
    temporary.replace(args.output)


def length_similarity(left: int, right: int) -> float:
    left, right = abs(left), abs(right)
    return min(left, right) / max(left, right) if max(left, right) else 1.0


def cluster(args: argparse.Namespace) -> None:
    with args.catalog.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    anchored = []
    unresolved = []
    for row in rows:
        try:
            anchored.append((row["chm13_chrom"], row["svtype"], int(row["chm13_start_0"]), int(row["svlen"]), row))
        except (KeyError, ValueError):
            unresolved.append(row)
    anchored.sort(key=lambda value: (value[0], value[1], value[2]))
    parent = list(range(len(anchored)))
    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index
    def union(left: int, right: int) -> None:
        a, b = find(left), find(right)
        if a != b:
            parent[b] = a
    active_by_group: dict[tuple[str, str], list[int]] = defaultdict(list)
    for index, (chrom, svtype, start, length, row) in enumerate(anchored):
        group = (chrom, svtype)
        active = [other for other in active_by_group[group] if start - anchored[other][2] <= args.breakpoint_distance]
        active_by_group[group] = active
        for other in active:
            _, _, other_start, other_length, other_row = anchored[other]
            if abs(start - other_start) > args.breakpoint_distance:
                continue
            if length_similarity(length, other_length) < args.length_similarity:
                continue
            left_sequence, right_sequence = row.get("inserted_sequence", ""), other_row.get("inserted_sequence", "")
            if left_sequence and right_sequence and SequenceMatcher(None, left_sequence, right_sequence).ratio() < args.sequence_similarity:
                continue
            union(index, other)
        active_by_group[group].append(index)
    roots = {}
    output = []
    for index, (_, _, _, _, row) in enumerate(anchored):
        root = find(index)
        if root not in roots:
            roots[root] = "SVC_" + hashlib.sha256(anchored[root][4]["stable_event_id"].encode()).hexdigest()[:12]
        output.append({
            "stable_event_id": row["stable_event_id"],
            "cluster_id": roots[root],
            "cluster_status": "clustered_CHM13_anchors",
        })
    for row in unresolved:
        output.append({
            "stable_event_id": row["stable_event_id"],
            "cluster_id": "",
            "cluster_status": "unresolved_or_unprojected_preserved",
        })
    output.sort(key=lambda row: row["stable_event_id"])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.output) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["stable_event_id", "cluster_id", "cluster_status"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    temporary.replace(args.output)


def main() -> None:
    args = parse_args()
    if args.command == "normalize":
        normalize(args)
    else:
        cluster(args)


if __name__ == "__main__":
    main()
