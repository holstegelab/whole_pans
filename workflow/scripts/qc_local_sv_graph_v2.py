#!/usr/bin/env python3
"""QC an incrementally rescued Minigraph rGFA and its cohort SV coverage.

The graph command compares V1 and V2 using stable rGFA segment attributes,
not numeric segment IDs, because Minigraph may renumber existing segments.
The report command combines structural checks with an accounting audit of the
existing cohort discovery evidence. It performs no assembly/read remapping or
new biological validation.
"""

import argparse
import csv
import gzip
import hashlib
import sys
from collections import Counter, defaultdict
from pathlib import Path


CANONICAL_CHROMS = {str(number) for number in range(1, 23)} | {
    "X",
    "Y",
    "M",
    "MT",
}
INTEGRITY_FIELDS = ["check", "expected", "observed", "status", "details"]
DELTA_SEGMENT_FIELDS = [
    "segment_id",
    "length",
    "sn",
    "so",
    "sr",
    "sequence_sha256",
    "reference_component",
    "anchor_status",
]
DELTA_LINK_FIELDS = [
    "left_segment",
    "left_orientation",
    "right_segment",
    "right_orientation",
    "overlap",
    "left_sn",
    "left_so",
    "left_sr",
    "right_sn",
    "right_so",
    "right_sr",
    "reference_component",
    "anchor_status",
]
BUBBLE_DELTA_FIELDS = [
    "stable_source",
    "start_0",
    "end_0",
    "reference_span_bp",
    "bubble_length_bp",
    "column_4",
    "column_5",
    "column_6",
    "column_7",
    "column_9",
    "column_10",
    "column_11",
    "reference_allele",
    "alternate_sha256",
]
RESCUE_QC_FIELDS = [
    "rescue_id",
    "catalog_event_id",
    "plan_status",
    "rescue_status",
    "combined_fasta_status",
    "chromosome",
    "position_0",
    "svtype",
    "svlen",
    "source_allele_spanning_read_count",
    "spanning_read_count",
    "source_allele_fraction",
    "selected_read_count",
    "local_contig_count",
    "selected_alignment_identity",
    "selected_target_coverage",
    "selected_event_discordance_bp",
    "evidence_warnings",
    "nearby_delta_bubble_count",
    "existing_high_confidence_graph_residual_evidence_count",
    "rest_cohort_high_confidence_graph_residual_evidence_count",
    "graph_status",
    "qc_status",
]


def allow_large_tsv_fields():
    """Raise csv's field limit for long SV alleles and catalog evidence."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def iter_tsv(path):
    with open_text(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def read_tsv(path):
    return list(iter_tsv(path))


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def int_value(value, default=0):
    if isinstance(value, bytes):
        value = value.decode("ascii")
    try:
        return int(value)
    except (TypeError, ValueError, OverflowError):
        try:
            return int(float(value))
        except (TypeError, ValueError, OverflowError):
            return default


def float_value(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def true_value(value):
    return str(value).lower() in {"1", "true", "yes"}


def normalize_chromosome(value):
    result = str(value or "")
    for prefix in ("CHM13.", "GRCh38.", "hg38."):
        if result.startswith(prefix):
            result = result[len(prefix) :]
            break
    if result and not result.startswith("chr") and result in CANONICAL_CHROMS:
        result = "chr" + result
    return result


def metric_table(path):
    rows = read_tsv(path)
    if not rows:
        return {}
    key_field = "metric" if "metric" in rows[0] else "key"
    return {row[key_field]: row["value"] for row in rows}


def parse_tags(fields):
    tags = {}
    for field in fields:
        parts = field.split(b":", 2)
        if len(parts) == 3:
            tags[parts[0].decode()] = parts[2].decode()
    return tags


def segment_record(fields):
    if len(fields) < 3:
        raise ValueError("Malformed GFA segment record")
    tags = parse_tags(fields[3:])
    sequence = fields[2]
    if sequence == b"*":
        length = int_value(tags.get("LN"), -1)
        sequence_hash = "*"
    else:
        length = len(sequence)
        sequence_hash = hashlib.sha256(sequence.upper()).hexdigest()
    stable_values = (
        tags.get("SN", ""),
        tags.get("SO", ""),
        tags.get("SR", ""),
        str(length),
        sequence_hash,
    )
    stable_key = hashlib.sha256("\t".join(stable_values).encode()).hexdigest()
    return {
        "segment_id": fields[1].decode(),
        "length": length,
        "sn": tags.get("SN", ""),
        "so": tags.get("SO", ""),
        "sr": tags.get("SR", ""),
        "ln": tags.get("LN", ""),
        "sequence_sha256": sequence_hash,
        "stable_key": stable_key,
    }


def link_signature(fields, segment_keys):
    if len(fields) < 6:
        raise ValueError("Malformed GFA link record")
    left = fields[1].decode()
    right = fields[3].decode()
    if left not in segment_keys or right not in segment_keys:
        return None
    values = (
        segment_keys[left],
        fields[2].decode(),
        segment_keys[right],
        fields[4].decode(),
    ) + tuple(field.decode() for field in fields[5:])
    return hashlib.sha256("\t".join(values).encode()).hexdigest()


class UnionFind:
    def __init__(self):
        self.parent = []
        self.size = []

    def add(self):
        index = len(self.parent)
        self.parent.append(index)
        self.size.append(1)
        return index

    def find(self, value):
        root = value
        while self.parent[root] != root:
            root = self.parent[root]
        while self.parent[value] != value:
            parent = self.parent[value]
            self.parent[value] = root
            value = parent
        return root

    def union(self, left, right):
        left_root = self.find(left)
        right_root = self.find(right)
        if left_root == right_root:
            return
        if self.size[left_root] < self.size[right_root]:
            left_root, right_root = right_root, left_root
        self.parent[right_root] = left_root
        self.size[left_root] += self.size[right_root]


def read_source_graph(path):
    file_hash = hashlib.sha256()
    segment_counts = Counter()
    segment_keys = {}
    link_counts = Counter()
    stable_source_bp = Counter()
    metrics = Counter()

    with Path(path).open("rb") as handle:
        for line_number, line in enumerate(handle, start=1):
            file_hash.update(line)
            fields = line.rstrip(b"\r\n").split(b"\t")
            if not fields:
                continue
            record_type = fields[0]
            if record_type == b"S":
                record = segment_record(fields)
                if record["segment_id"] in segment_keys:
                    raise ValueError(
                        f"Duplicate segment ID in frozen graph at line {line_number}: "
                        f"{record['segment_id']}"
                    )
                segment_keys[record["segment_id"]] = record["stable_key"]
                segment_counts[record["stable_key"]] += 1
                if record["sn"] and record["sr"]:
                    stable_source_bp[(record["sn"], record["sr"])] += record["length"]
                metrics["segments"] += 1
                metrics["segment_bp"] += record["length"]
            elif record_type == b"L":
                signature = link_signature(fields, segment_keys)
                if signature is None:
                    raise ValueError(
                        f"Frozen graph link refers to an unavailable segment at line "
                        f"{line_number}"
                    )
                link_counts[signature] += 1
                metrics["links"] += 1
            elif record_type in {b"P", b"W"}:
                metrics["paths"] += 1
    return {
        "sha256": file_hash.hexdigest(),
        "segment_counts": segment_counts,
        "link_counts": link_counts,
        "stable_source_bp": stable_source_bp,
        "metrics": metrics,
    }


def anchor_label(anchors):
    if len(anchors) == 1:
        return next(iter(anchors)), "ANCHORED"
    if not anchors:
        return "", "UNANCHORED"
    return ";".join(sorted(anchors)), "CROSS_REFERENCE_COMPONENT"


def read_rescued_graph(path, source):
    file_hash = hashlib.sha256()
    remaining_segments = source["segment_counts"].copy()
    remaining_links = source["link_counts"].copy()
    metrics = Counter()
    segment_keys = {}
    segment_records = {}
    segment_indices = {}
    reference_sources = []
    new_segment_rows = []
    new_segment_indices = []
    new_link_rows = []
    pending_edges = []
    duplicate_segments = 0
    missing_tags = Counter()
    ln_mismatches = 0
    stable_name_ranks = defaultdict(set)
    stable_source_bp = Counter()
    union_find = UnionFind()

    with Path(path).open("rb") as handle:
        for line_number, line in enumerate(handle, start=1):
            file_hash.update(line)
            fields = line.rstrip(b"\r\n").split(b"\t")
            if not fields:
                continue
            record_type = fields[0]
            if record_type == b"S":
                record = segment_record(fields)
                segment_id = record["segment_id"]
                if segment_id in segment_keys:
                    duplicate_segments += 1
                    continue
                index = union_find.add()
                segment_keys[segment_id] = record["stable_key"]
                segment_records[segment_id] = record
                segment_indices[segment_id] = index
                reference_sources.append("")
                for tag in ("SN", "SO", "SR", "LN"):
                    if record[tag.lower()] == "":
                        missing_tags[tag] += 1
                if record["ln"] and int_value(record["ln"], -1) != record["length"]:
                    ln_mismatches += 1
                if record["sn"] and record["sr"]:
                    stable_name_ranks[record["sn"]].add(record["sr"])
                    stable_source_bp[(record["sn"], record["sr"])] += record["length"]
                if record["sr"] == "0":
                    reference_sources[index] = normalize_chromosome(record["sn"])
                if remaining_segments[record["stable_key"]] > 0:
                    remaining_segments[record["stable_key"]] -= 1
                else:
                    new_segment_indices.append(index)
                    new_segment_rows.append(dict(record))
                metrics["segments"] += 1
                metrics["segment_bp"] += record["length"]
            elif record_type == b"L":
                metrics["links"] += 1
                left = fields[1].decode() if len(fields) > 1 else ""
                right = fields[3].decode() if len(fields) > 3 else ""
                signature = link_signature(fields, segment_keys)
                if signature is None:
                    pending_edges.append((left, right, line_number))
                    is_new = True
                else:
                    is_new = remaining_links[signature] <= 0
                    if not is_new:
                        remaining_links[signature] -= 1
                if left in segment_indices and right in segment_indices:
                    union_find.union(segment_indices[left], segment_indices[right])
                if is_new:
                    left_record = segment_records.get(left, {})
                    right_record = segment_records.get(right, {})
                    new_link_rows.append(
                        {
                            "left_segment": left,
                            "left_orientation": fields[2].decode() if len(fields) > 2 else "",
                            "right_segment": right,
                            "right_orientation": fields[4].decode() if len(fields) > 4 else "",
                            "overlap": fields[5].decode() if len(fields) > 5 else "",
                            "left_sn": left_record.get("sn", ""),
                            "left_so": left_record.get("so", ""),
                            "left_sr": left_record.get("sr", ""),
                            "right_sn": right_record.get("sn", ""),
                            "right_so": right_record.get("so", ""),
                            "right_sr": right_record.get("sr", ""),
                            "_left_index": segment_indices.get(left),
                        }
                    )
            elif record_type in {b"P", b"W"}:
                metrics["paths"] += 1

    dangling_endpoints = set()
    for left, right, _line_number in pending_edges:
        if left not in segment_indices:
            dangling_endpoints.add(left)
        if right not in segment_indices:
            dangling_endpoints.add(right)
        if left in segment_indices and right in segment_indices:
            union_find.union(segment_indices[left], segment_indices[right])

    component_anchors = defaultdict(set)
    for index, source_name in enumerate(reference_sources):
        if source_name:
            component_anchors[union_find.find(index)].add(source_name)

    for row, index in zip(new_segment_rows, new_segment_indices):
        anchors = component_anchors.get(union_find.find(index), set())
        row["reference_component"], row["anchor_status"] = anchor_label(anchors)
    for row in new_link_rows:
        index = row.pop("_left_index", None)
        anchors = component_anchors.get(union_find.find(index), set()) if index is not None else set()
        row["reference_component"], row["anchor_status"] = anchor_label(anchors)

    missing_source_segments = sum(remaining_segments.values())
    missing_source_links = sum(remaining_links.values())
    rank_conflicts = sum(len(ranks) > 1 for ranks in stable_name_ranks.values())
    cross_reference_components = sum(
        len(anchors) > 1 for anchors in component_anchors.values()
    )
    component_count = len({union_find.find(index) for index in range(len(union_find.parent))})
    missing_stable_source_bp = sum(
        max(0, length - stable_source_bp.get(source_key, 0))
        for source_key, length in source["stable_source_bp"].items()
    )
    missing_stable_sources = sum(
        stable_source_bp.get(source_key, 0) < length
        for source_key, length in source["stable_source_bp"].items()
    )
    return {
        "sha256": file_hash.hexdigest(),
        "metrics": metrics,
        "delta_segments": new_segment_rows,
        "delta_links": new_link_rows,
        "missing_source_segments": missing_source_segments,
        "missing_source_links": missing_source_links,
        "missing_stable_source_bp": missing_stable_source_bp,
        "missing_stable_sources": missing_stable_sources,
        "duplicate_segments": duplicate_segments,
        "dangling_endpoints": dangling_endpoints,
        "missing_tags": missing_tags,
        "ln_mismatches": ln_mismatches,
        "rank_conflicts": rank_conflicts,
        "cross_reference_components": cross_reference_components,
        "components": component_count,
        "unanchored_delta_segments": sum(
            row["anchor_status"] != "ANCHORED" for row in new_segment_rows
        ),
        "unanchored_delta_links": sum(
            row["anchor_status"] != "ANCHORED" for row in new_link_rows
        ),
    }


def build_log_rgfa_warnings(path):
    count = 0
    with Path(path).open(errors="replace") as handle:
        for line in handle:
            if "inconsistent rGFA" in line or "associated with different ranks" in line:
                count += 1
    return count


def run_graph_qc(args):
    source = read_source_graph(args.frozen_graph)
    final = read_rescued_graph(args.rescued_graph, source)
    summary = metric_table(args.graph_summary)
    inventory = metric_table(args.frozen_inventory)
    provenance = metric_table(args.provenance)
    checks = []

    def add(check, expected, observed, status, details=""):
        checks.append(
            {
                "check": check,
                "expected": expected,
                "observed": observed,
                "status": status,
                "details": details,
            }
        )

    add(
        "frozen_graph_sha256_matches_inventory",
        inventory.get("sha256", ""),
        source["sha256"],
        "PASS" if inventory.get("sha256") == source["sha256"] else "FAIL",
    )
    add(
        "rescued_graph_sha256_matches_provenance",
        provenance.get("final_graph_sha256", ""),
        final["sha256"],
        "PASS" if provenance.get("final_graph_sha256") == final["sha256"] else "FAIL",
    )
    for metric in ("segments", "links", "paths", "segment_bp"):
        observed_source = source["metrics"][metric]
        observed_final = final["metrics"][metric]
        expected_source = int_value(summary.get(f"source_{metric}"), -1)
        expected_final = int_value(summary.get(f"final_{metric}"), -1)
        add(
            f"source_{metric}_matches_summary",
            expected_source,
            observed_source,
            "PASS" if expected_source == observed_source else "FAIL",
        )
        add(
            f"final_{metric}_matches_summary",
            expected_final,
            observed_final,
            "PASS" if expected_final == observed_final else "FAIL",
        )
    add(
        "all_frozen_stable_segment_records_preserved",
        0,
        final["missing_source_segments"],
        "PASS" if final["missing_source_segments"] == 0 else "WARN",
        "Exact stable records may change when Minigraph splits existing segments.",
    )
    add(
        "all_frozen_stable_links_preserved",
        0,
        final["missing_source_links"],
        "PASS" if final["missing_source_links"] == 0 else "WARN",
        "Exact links may change when their endpoint segments are split.",
    )
    add(
        "frozen_stable_source_bp_preserved",
        0,
        final["missing_stable_source_bp"],
        "PASS" if final["missing_stable_source_bp"] == 0 else "FAIL",
        f"stable sources with missing bp: {final['missing_stable_sources']}",
    )
    expected_delta_segments = int_value(summary.get("delta_segments"), -1)
    expected_delta_links = int_value(summary.get("delta_links"), -1)
    add(
        "delta_segment_count_matches_summary",
        expected_delta_segments,
        final["metrics"]["segments"] - source["metrics"]["segments"],
        "PASS"
        if expected_delta_segments
        == final["metrics"]["segments"] - source["metrics"]["segments"]
        else "FAIL",
        f"new or re-segmented stable records listed: {len(final['delta_segments'])}",
    )
    add(
        "delta_link_count_matches_summary",
        expected_delta_links,
        final["metrics"]["links"] - source["metrics"]["links"],
        "PASS"
        if expected_delta_links
        == final["metrics"]["links"] - source["metrics"]["links"]
        else "FAIL",
        f"new or reconnected stable links listed: {len(final['delta_links'])}",
    )
    for name, observed in (
        ("duplicate_segment_ids", final["duplicate_segments"]),
        ("dangling_link_endpoints", len(final["dangling_endpoints"])),
        ("segment_ln_mismatches", final["ln_mismatches"]),
        ("stable_name_rank_conflicts", final["rank_conflicts"]),
        ("cross_reference_components", final["cross_reference_components"]),
        ("unanchored_delta_segments", final["unanchored_delta_segments"]),
        ("unanchored_delta_links", final["unanchored_delta_links"]),
    ):
        add(name, 0, observed, "PASS" if observed == 0 else "FAIL")
    for tag in ("SN", "SO", "SR", "LN"):
        observed = final["missing_tags"][tag]
        add(
            f"segments_missing_{tag}_tag",
            0,
            observed,
            "PASS" if observed == 0 else "FAIL",
        )
    rgfa_warnings = build_log_rgfa_warnings(args.build_log)
    add(
        "build_log_rgfa_consistency_warnings",
        0,
        rgfa_warnings,
        "PASS" if rgfa_warnings == 0 else "FAIL",
    )

    metrics = [
        {"metric": "frozen_graph_sha256", "value": source["sha256"]},
        {"metric": "rescued_graph_sha256", "value": final["sha256"]},
        {"metric": "connected_components", "value": final["components"]},
        {
            "metric": "net_delta_segments",
            "value": final["metrics"]["segments"] - source["metrics"]["segments"],
        },
        {
            "metric": "net_delta_links",
            "value": final["metrics"]["links"] - source["metrics"]["links"],
        },
        {"metric": "new_or_resegmented_records", "value": len(final["delta_segments"])},
        {"metric": "new_or_reconnected_links", "value": len(final["delta_links"])},
        {
            "metric": "integrity_failures",
            "value": sum(row["status"] == "FAIL" for row in checks),
        },
    ]
    write_tsv(args.integrity, checks, INTEGRITY_FIELDS)
    write_tsv(args.delta_segments, final["delta_segments"], DELTA_SEGMENT_FIELDS)
    write_tsv(args.delta_links, final["delta_links"], DELTA_LINK_FIELDS)
    write_tsv(args.metrics, metrics, ["metric", "value"])


def bubble_signature(fields):
    if len(fields) < 14:
        raise ValueError(f"Expected at least 14 gfatools bubble fields, found {len(fields)}")
    alternate_hash = hashlib.sha256(fields[13].upper()).hexdigest()
    # Column 12 contains numeric node IDs, which Minigraph may renumber.
    stable_fields = fields[:11] + fields[12:13] + [alternate_hash.encode()]
    return hashlib.sha256(b"\t".join(stable_fields)).hexdigest(), alternate_hash


def bubble_delta(frozen_path, rescued_path):
    frozen = Counter()
    with Path(frozen_path).open("rb") as handle:
        for line in handle:
            fields = line.rstrip(b"\r\n").split(b"\t")
            if fields and fields[0]:
                signature, _alternate_hash = bubble_signature(fields)
                frozen[signature] += 1

    delta = []
    with Path(rescued_path).open("rb") as handle:
        for line in handle:
            fields = line.rstrip(b"\r\n").split(b"\t")
            if not fields or not fields[0]:
                continue
            signature, alternate_hash = bubble_signature(fields)
            if frozen[signature] > 0:
                frozen[signature] -= 1
                continue
            # Required bubble coordinates must never silently become zero.
            start = int(fields[1])
            end = int(fields[2])
            length = int(fields[7])
            if start < 0 or end < start or length < 0:
                raise ValueError(f"Invalid bubble interval: {fields[:8]!r}")
            delta.append(
                {
                    "stable_source": fields[0].decode(),
                    "start_0": start,
                    "end_0": end,
                    "reference_span_bp": max(0, end - start),
                    "bubble_length_bp": length,
                    "column_4": fields[3].decode(),
                    "column_5": fields[4].decode(),
                    "column_6": fields[5].decode(),
                    "column_7": fields[6].decode(),
                    "column_9": fields[8].decode(),
                    "column_10": fields[9].decode(),
                    "column_11": fields[10].decode(),
                    "reference_allele": fields[12].decode(),
                    "alternate_sha256": alternate_hash,
                }
            )
    return delta, sum(frozen.values())


def length_compatible(left, right, similarity):
    left = abs(int_value(left))
    right = abs(int_value(right))
    if not left or not right:
        return True
    return min(left, right) / max(left, right) >= similarity


def bubble_matches_plan(plan, bubble, breakpoint_distance, length_similarity):
    chromosome = normalize_chromosome(bubble["stable_source"])
    if chromosome != normalize_chromosome(plan.get("chromosome")):
        return False
    position = int_value(plan.get("position_0"))
    distance = min(
        abs(position - int_value(bubble["start_0"])),
        abs(position - int_value(bubble["end_0"])),
    )
    bubble_size = max(
        int_value(bubble["reference_span_bp"]),
        int_value(bubble["bubble_length_bp"]),
    )
    return distance <= breakpoint_distance and length_compatible(
        bubble_size, plan.get("svlen"), length_similarity
    )


def evidence_warnings(row, args):
    warnings = []
    allele_reads = int_value(row.get("source_allele_spanning_read_count"))
    spanning_reads = int_value(row.get("spanning_read_count"))
    allele_fraction = allele_reads / spanning_reads if spanning_reads else 0.0
    if allele_reads < args.min_allele_reads_warn:
        warnings.append(f"ALLELE_READS<{args.min_allele_reads_warn}")
    if allele_fraction < args.min_allele_fraction_warn:
        warnings.append(f"ALLELE_FRACTION<{args.min_allele_fraction_warn}")
    if int_value(row.get("selected_read_count")) > args.max_selected_reads_warn:
        warnings.append(f"SELECTED_READS>{args.max_selected_reads_warn}")
    if int_value(row.get("local_contig_count")) > args.max_local_contigs_warn:
        warnings.append(f"LOCAL_CONTIGS>{args.max_local_contigs_warn}")
    if float_value(row.get("selected_target_coverage")) < args.min_target_coverage_warn:
        warnings.append(f"TARGET_COVERAGE<{args.min_target_coverage_warn}")
    return warnings, allele_fraction


def run_report(args):
    plan_rows = read_tsv(args.plan)
    plans = {row["rescue_id"]: row for row in plan_rows}
    rescue_summary = {row["rescue_id"]: row for row in read_tsv(args.rescue_summary)}
    graph_inputs = read_tsv(args.graph_inputs)
    included_rescue_ids = {row["rescue_id"] for row in graph_inputs}
    included_catalog_ids = {
        row["catalog_event_id"] for row in graph_inputs if row.get("catalog_event_id")
    }
    rescue_source_assemblies = {
        plans[rescue_id]["assembly_id"]
        for rescue_id in included_rescue_ids
        if rescue_id in plans and plans[rescue_id].get("assembly_id")
    }

    catalog_rows = read_tsv(args.catalog)
    catalog_ids = {row["event_id"] for row in catalog_rows}
    evidence_to_catalog = {}
    conflicting_evidence_ids = set()
    for row in catalog_rows:
        catalog_event_id = row["event_id"]
        for evidence_id in filter(None, row.get("evidence_ids", "").split(";")):
            previous = evidence_to_catalog.setdefault(evidence_id, catalog_event_id)
            if previous != catalog_event_id:
                conflicting_evidence_ids.add(evidence_id)

    manifest_rows = read_tsv(args.assembly_manifest)
    non_graph_assemblies = {
        row["assembly_id"]
        for row in manifest_rows
        if not true_value(row.get("graph_member"))
    }
    rest_cohort_assemblies = non_graph_assemblies - rescue_source_assemblies

    catalog_hc_evidence = Counter()
    catalog_rest_hc_evidence = Counter()
    all_hc_event_ids = set()
    rest_hc_event_ids = set()
    unaccounted_rest_event_ids = set()
    unaccounted_rest_evidence_ids = set()
    unmapped_rest_evidence_ids = set()
    rest_affected_assemblies = set()
    rest_affected_samples = set()
    evidence_assemblies_missing_from_manifest = set()
    all_hc_evidence_count = 0
    rest_hc_evidence_count = 0
    accounted_rest_hc_evidence_count = 0

    for evidence in iter_tsv(args.catalog_evidence):
        if (
            evidence.get("discovery_method") != "graph_residual"
            or evidence.get("confidence_tier") != "HIGH_CONFIDENCE"
            or true_value(evidence.get("graph_member"))
            or abs(int_value(evidence.get("svlen"))) < args.min_sv_size
        ):
            continue
        evidence_id = evidence.get("evidence_id", "")
        catalog_event_id = evidence_to_catalog.get(evidence_id, "")
        assembly_id = evidence.get("assembly_id", "")
        all_hc_evidence_count += 1
        if catalog_event_id:
            all_hc_event_ids.add(catalog_event_id)
            catalog_hc_evidence[catalog_event_id] += 1
        if assembly_id not in non_graph_assemblies:
            evidence_assemblies_missing_from_manifest.add(assembly_id)
            continue
        if assembly_id not in rest_cohort_assemblies:
            continue
        rest_hc_evidence_count += 1
        rest_affected_assemblies.add(assembly_id)
        rest_affected_samples.add(evidence.get("sample_id", ""))
        if catalog_event_id:
            rest_hc_event_ids.add(catalog_event_id)
            catalog_rest_hc_evidence[catalog_event_id] += 1
        else:
            unmapped_rest_evidence_ids.add(evidence_id)
        if catalog_event_id in included_catalog_ids:
            accounted_rest_hc_evidence_count += 1
        else:
            unaccounted_rest_evidence_ids.add(evidence_id)
            if catalog_event_id:
                unaccounted_rest_event_ids.add(catalog_event_id)

    bubbles, removed_bubbles = bubble_delta(args.frozen_bubbles, args.rescued_bubbles)
    write_tsv(args.bubble_delta, bubbles, BUBBLE_DELTA_FIELDS)

    rescue_qc = []
    for rescue_id in sorted(plans):
        plan = plans[rescue_id]
        summary = rescue_summary.get(rescue_id, {})
        if summary.get("status") == "PASS":
            warnings, allele_fraction = evidence_warnings(summary, args)
            allele_fraction_value = f"{allele_fraction:.8f}"
        else:
            warnings = []
            allele_fraction_value = ""
        nearby_bubbles = sum(
            bubble_matches_plan(
                plan, bubble, args.breakpoint_distance, args.length_similarity
            )
            for bubble in bubbles
        )
        catalog_event_id = plan.get("catalog_event_id", "")
        high_evidence = catalog_hc_evidence[catalog_event_id]
        rest_high_evidence = catalog_rest_hc_evidence[catalog_event_id]
        if plan.get("plan_status") != "READY":
            graph_status = "NOT_ATTEMPTED"
            qc_status = "INFO"
        elif rescue_id not in included_rescue_ids:
            graph_status = "RESCUE_NOT_INCLUDED"
            qc_status = "FAIL"
        elif nearby_bubbles:
            graph_status = "INCLUDED_WITH_NEARBY_DELTA_BUBBLE"
            qc_status = "WARN" if warnings else "PASS"
        else:
            graph_status = "INCLUDED_WITHOUT_UNIQUE_DELTA_BUBBLE"
            qc_status = "WARN"
        rescue_qc.append(
            {
                "rescue_id": rescue_id,
                "catalog_event_id": plan.get("catalog_event_id", ""),
                "plan_status": plan.get("plan_status", ""),
                "rescue_status": summary.get("status", ""),
                "combined_fasta_status": summary.get("combined_fasta_status", ""),
                "chromosome": plan.get("chromosome", ""),
                "position_0": plan.get("position_0", ""),
                "svtype": plan.get("svtype", ""),
                "svlen": plan.get("svlen", ""),
                "source_allele_spanning_read_count": summary.get(
                    "source_allele_spanning_read_count", ""
                ),
                "spanning_read_count": summary.get("spanning_read_count", ""),
                "source_allele_fraction": allele_fraction_value,
                "selected_read_count": summary.get("selected_read_count", ""),
                "local_contig_count": summary.get("local_contig_count", ""),
                "selected_alignment_identity": summary.get(
                    "selected_alignment_identity", ""
                ),
                "selected_target_coverage": summary.get(
                    "selected_target_coverage", ""
                ),
                "selected_event_discordance_bp": summary.get(
                    "selected_event_discordance_bp", ""
                ),
                "evidence_warnings": ";".join(warnings),
                "nearby_delta_bubble_count": nearby_bubbles,
                "existing_high_confidence_graph_residual_evidence_count": high_evidence,
                "rest_cohort_high_confidence_graph_residual_evidence_count": rest_high_evidence,
                "graph_status": graph_status,
                "qc_status": qc_status,
            }
        )
    write_tsv(args.rescue_qc, rescue_qc, RESCUE_QC_FIELDS)

    integrity = read_tsv(args.integrity)
    integrity_failures = sum(row.get("status") == "FAIL" for row in integrity)
    rescue_status_counts = Counter(row["qc_status"] for row in rescue_qc)
    invalid_graph_input_ids = included_catalog_ids - catalog_ids
    invalid_graph_input_rescue_ids = included_rescue_ids - set(plans)
    if (
        unaccounted_rest_evidence_ids
        or unmapped_rest_evidence_ids
        or conflicting_evidence_ids
        or invalid_graph_input_ids
        or invalid_graph_input_rescue_ids
        or evidence_assemblies_missing_from_manifest
    ):
        completeness_status = "FAIL"
    else:
        completeness_status = "PASS"

    if (
        integrity_failures
        or rescue_status_counts["FAIL"]
        or completeness_status == "FAIL"
    ):
        overall_status = "FAIL"
    elif (
        completeness_status == "WARN"
        or rescue_status_counts["WARN"]
        or any(row.get("status") == "WARN" for row in integrity)
    ):
        overall_status = "WARN"
    else:
        overall_status = "PASS"

    cohort_qc = [
        {"metric": "catalog_events", "value": len(catalog_ids), "status": "INFO"},
        {"metric": "included_rescue_catalog_events", "value": len(included_catalog_ids), "status": "INFO"},
        {"metric": "non_graph_cohort_assemblies", "value": len(non_graph_assemblies), "status": "INFO"},
        {"metric": "rest_cohort_assemblies", "value": len(rest_cohort_assemblies), "status": "INFO"},
        {"metric": "existing_non_graph_high_confidence_graph_residual_evidence", "value": all_hc_evidence_count, "status": "INFO"},
        {"metric": "existing_non_graph_high_confidence_catalog_events", "value": len(all_hc_event_ids), "status": "INFO"},
        {"metric": "rest_cohort_high_confidence_graph_residual_evidence", "value": rest_hc_evidence_count, "status": "INFO"},
        {"metric": "rest_cohort_high_confidence_catalog_events", "value": len(rest_hc_event_ids), "status": "INFO"},
        {"metric": "rest_cohort_accounted_high_confidence_evidence", "value": accounted_rest_hc_evidence_count, "status": "INFO"},
        {"metric": "rest_cohort_unaccounted_high_confidence_evidence", "value": len(unaccounted_rest_evidence_ids), "status": completeness_status},
        {"metric": "rest_cohort_unaccounted_catalog_events", "value": len(unaccounted_rest_event_ids), "status": completeness_status},
        {"metric": "rest_cohort_unaccounted_catalog_event_ids", "value": ";".join(sorted(unaccounted_rest_event_ids)), "status": completeness_status},
        {"metric": "rest_cohort_unmapped_evidence_ids", "value": len(unmapped_rest_evidence_ids), "status": "PASS" if not unmapped_rest_evidence_ids else "FAIL"},
        {"metric": "rest_cohort_unmapped_evidence_id_list", "value": ";".join(sorted(unmapped_rest_evidence_ids)), "status": "PASS" if not unmapped_rest_evidence_ids else "FAIL"},
        {"metric": "catalog_evidence_mapping_conflicts", "value": len(conflicting_evidence_ids), "status": "PASS" if not conflicting_evidence_ids else "FAIL"},
        {"metric": "invalid_graph_input_catalog_ids", "value": len(invalid_graph_input_ids), "status": "PASS" if not invalid_graph_input_ids else "FAIL"},
        {"metric": "invalid_graph_input_rescue_ids", "value": len(invalid_graph_input_rescue_ids), "status": "PASS" if not invalid_graph_input_rescue_ids else "FAIL"},
        {"metric": "evidence_assemblies_missing_from_manifest", "value": len(evidence_assemblies_missing_from_manifest), "status": "PASS" if not evidence_assemblies_missing_from_manifest else "FAIL"},
        {"metric": "evidence_assemblies_missing_from_manifest_list", "value": ";".join(sorted(evidence_assemblies_missing_from_manifest)), "status": "PASS" if not evidence_assemblies_missing_from_manifest else "FAIL"},
        {"metric": "rest_cohort_assemblies_with_existing_hc_evidence", "value": len(rest_affected_assemblies), "status": "INFO"},
        {"metric": "rest_cohort_samples_with_existing_hc_evidence", "value": len(rest_affected_samples), "status": "INFO"},
        {"metric": "new_or_modified_bubbles", "value": len(bubbles), "status": "INFO"},
        {"metric": "removed_or_modified_frozen_bubbles", "value": removed_bubbles, "status": "PASS" if removed_bubbles == 0 else "WARN"},
        {"metric": "structural_integrity_failures", "value": integrity_failures, "status": "PASS" if integrity_failures == 0 else "FAIL"},
        {"metric": "rescue_qc_failures", "value": rescue_status_counts["FAIL"], "status": "PASS" if rescue_status_counts["FAIL"] == 0 else "FAIL"},
        {"metric": "overall_qc_status", "value": overall_status, "status": overall_status},
    ]
    write_tsv(args.cohort_qc, cohort_qc, ["metric", "value", "status"])

    if completeness_status == "PASS":
        completeness_statement = (
            f"All previously detected high-confidence graph-residual SV evidence "
            f">= {args.min_sv_size} bp in the rest of the non-graph cohort is "
            "accounted for by catalog events with rescue inputs included in V2."
        )
    else:
        completeness_statement = (
            "The accounting claim is NOT supported: "
            f"{len(unaccounted_rest_evidence_ids)} existing high-confidence "
            f"graph-residual evidence rows from {len(unaccounted_rest_event_ids)} "
            "mapped catalog events are not represented by an included V2 rescue "
            f"input; unmapped evidence rows={len(unmapped_rest_evidence_ids)}."
        )

    lines = [
        "# Local SV graph V2 QC report",
        "",
        f"- Overall QC status: **{overall_status}**.",
        f"- Structural integrity failures: {integrity_failures}.",
        f"- Rescue-level PASS/WARN/FAIL/INFO: {rescue_status_counts['PASS']}/"
        f"{rescue_status_counts['WARN']}/{rescue_status_counts['FAIL']}/"
        f"{rescue_status_counts['INFO']}.",
        f"- New or modified normalized bubbles: {len(bubbles)}.",
        "",
        "## Existing cohort-discovery accounting",
        "",
        f"- {completeness_statement}",
        f"- Non-graph cohort assemblies: {len(non_graph_assemblies)}; "
        f"rest-of-cohort after excluding included-rescue sources: "
        f"{len(rest_cohort_assemblies)}.",
        f"- Existing high-confidence graph-residual evidence in the rest cohort: "
        f"{rest_hc_evidence_count}; accounted: {accounted_rest_hc_evidence_count}; "
        f"unaccounted: {len(unaccounted_rest_evidence_ids)}.",
        f"- Unaccounted catalog events: {len(unaccounted_rest_event_ids)}.",
        f"- Evidence/catalog/manifest consistency anomalies: unmapped evidence="
        f"{len(unmapped_rest_evidence_ids)}, mapping conflicts="
        f"{len(conflicting_evidence_ids)}, missing assemblies="
        f"{len(evidence_assemblies_missing_from_manifest)}.",
        f"- Existing-evidence accounting status: {completeness_status}.",
        "",
        "This is an accounting check over the existing V1 discovery catalog, "
        "not a new V2 cohort screen. It can show that no previously detected "
        "high-confidence graph-residual SV was left without an included rescue, "
        "but it cannot prove that V2 absorbs those alleles or that no additional "
        "SV exists. Those conclusions require remapping or independent validation.",
        "",
        "## Scope",
        "",
        "- Included: GFA integrity, stable V1-to-V2 preservation, graph delta, "
        "normalized bubble delta, rescue evidence QC, and accounting of existing "
        "cohort graph-residual evidence.",
        "- Not included: assembly remapping, sequencing-read remapping, independent "
        "biological validation, or downstream genotyping performance.",
        "",
        "## Tables",
        "",
        "- `graph_integrity.tsv`",
        "- `delta_segments.tsv`",
        "- `delta_links.tsv`",
        "- `bubble_delta.tsv`",
        "- `rescued_sv_graph_qc.tsv`",
        "- `cohort_completeness.tsv`",
    ]
    report = Path(args.report)
    report.parent.mkdir(parents=True, exist_ok=True)
    report.write_text("\n".join(lines) + "\n")


def add_graph_parser(subparsers):
    parser = subparsers.add_parser("graph", help="Check graph integrity and V1-to-V2 delta")
    parser.add_argument("--frozen-graph", required=True)
    parser.add_argument("--rescued-graph", required=True)
    parser.add_argument("--graph-summary", required=True)
    parser.add_argument("--frozen-inventory", required=True)
    parser.add_argument("--provenance", required=True)
    parser.add_argument("--build-log", required=True)
    parser.add_argument("--integrity", required=True)
    parser.add_argument("--delta-segments", required=True)
    parser.add_argument("--delta-links", required=True)
    parser.add_argument("--metrics", required=True)


def add_report_parser(subparsers):
    parser = subparsers.add_parser(
        "report", help="Combine graph QC and existing cohort-evidence accounting"
    )
    for name in (
        "plan",
        "rescue_summary",
        "graph_inputs",
        "catalog",
        "catalog_evidence",
        "assembly_manifest",
        "integrity",
        "frozen_bubbles",
        "rescued_bubbles",
        "bubble_delta",
        "rescue_qc",
        "cohort_qc",
        "report",
    ):
        parser.add_argument(f"--{name.replace('_', '-')}", required=True)
    parser.add_argument("--min-sv-size", type=int, default=50)
    parser.add_argument("--breakpoint-distance", type=int, default=500)
    parser.add_argument("--length-similarity", type=float, default=0.70)
    parser.add_argument("--min-allele-reads-warn", type=int, default=5)
    parser.add_argument("--min-allele-fraction-warn", type=float, default=0.10)
    parser.add_argument("--max-selected-reads-warn", type=int, default=100000)
    parser.add_argument("--max-local-contigs-warn", type=int, default=10000)
    parser.add_argument("--min-target-coverage-warn", type=float, default=0.80)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_graph_parser(subparsers)
    add_report_parser(subparsers)
    return parser.parse_args()


def main():
    allow_large_tsv_fields()
    args = parse_args()
    if args.command == "graph":
        run_graph_qc(args)
    elif args.command == "report":
        run_report(args)
    else:
        raise AssertionError(args.command)


if __name__ == "__main__":
    main()
