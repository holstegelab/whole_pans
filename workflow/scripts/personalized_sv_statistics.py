#!/usr/bin/env python3
"""Calculate mapping and targeted SV statistics for direct graph mapping."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import random
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path


MAPPING_FIELDS = [
    "sample_id",
    "sample_group",
    "batch",
    "target",
    "source_graph",
    "graph_checksum",
    "mapper_version",
    "input_reads",
    "input_read_pairs",
    "input_bases",
    "primary_mapped_reads",
    "mapped_read_fraction",
    "mapq_0_fraction",
    "mapq_ge_5_fraction",
    "mapq_ge_20_fraction",
    "mapq_ge_40_fraction",
    "mean_mapq",
    "median_mapq",
    "mean_alignment_identity",
    "fraction_read_bases_aligned",
    "clipped_bases_per_read",
    "secondary_multimapping_fraction",
    "proper_pair_fraction",
    "common_backbone_mapped_fraction",
    "common_backbone_breadth_fraction",
    "common_backbone_mean_depth",
    "mapping_mode",
    "mapping_graph_size_bytes",
    "mapping_index_size_bytes",
    "index_wall_time_seconds",
    "index_peak_rss_mb",
    "mapping_wall_time_seconds",
    "mapping_cpu_time_seconds",
    "mapping_peak_rss_mb",
    "gaf_size_bytes",
]

EVIDENCE_FIELDS = [
    "sample_id",
    "sample_group",
    "batch",
    "target",
    "source_graph",
    "graph_checksum",
    "mapper_version",
    "sv_id",
    "call_class",
    "target_present",
    "target_callable_status",
    "evidence_available",
    "counting_unit",
    "min_mapq",
    "min_junction_fragments",
    "min_reference_fragments",
    "left_alt_junction_reads",
    "right_alt_junction_reads",
    "minimum_alt_junction_reads",
    "reference_path_reads",
    "insertion_sequence_reads",
    "locus_mapq_qualified_fragments",
    "genome_mapq_qualified_alignments",
    "conflicting_fragments",
    "alt_junction_balance",
    "ambiguous_or_conflicting",
    "note",
]

COMPARE_METRICS = [
    "mapped_read_fraction",
    "mapq_ge_5_fraction",
    "mapq_ge_20_fraction",
    "mapq_ge_40_fraction",
    "mean_alignment_identity",
    "fraction_read_bases_aligned",
    "secondary_multimapping_fraction",
    "proper_pair_fraction",
    "common_backbone_mapped_fraction",
    "common_backbone_breadth_fraction",
    "common_backbone_mean_depth",
    "mapping_wall_time_seconds",
    "mapping_peak_rss_mb",
]

CONTRASTS = [
    ("local_v1_minus_hprc", "local_v1", "hprc"),
    ("local_v2_minus_hprc", "local_v2", "hprc"),
    ("hprc_plus_local_minus_hprc", "hprc_plus_local", "hprc"),
    ("local_v2_minus_local_v1", "local_v2", "local_v1"),
]

DECISIVE_CALLS = {"ALT_SUPPORTED", "HET_LIKE", "REF_SUPPORTED"}
ALT_CALLS = {"ALT_SUPPORTED", "HET_LIKE"}
FRACTION_METRICS = {metric for metric in COMPARE_METRICS if "fraction" in metric} | {"mean_alignment_identity"}

# Graph catalog fields can contain long path or segment lists.  The csv module's
# 128 KiB default is too small for valid records in whole-genome catalogs.
csv.field_size_limit(sys.maxsize)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    mapping = subparsers.add_parser("mapping")
    mapping.add_argument("--gaf", required=True, type=Path)
    mapping.add_argument("--fastq-stats", required=True, type=Path)
    mapping.add_argument("--coverage", required=True, type=Path)
    mapping.add_argument("--flagstat", required=True, type=Path)
    mapping.add_argument("--mapping-graph", required=True, type=Path)
    mapping.add_argument("--index", action="append", default=[])
    mapping.add_argument("--inventory", required=True, type=Path)
    mapping.add_argument("--benchmark", required=True, type=Path)
    mapping.add_argument("--index-benchmark", required=True, type=Path)
    mapping.add_argument("--sample", required=True)
    mapping.add_argument("--sample-group", required=True)
    mapping.add_argument("--batch", required=True)
    mapping.add_argument("--target", required=True)
    mapping.add_argument("--mapper-version", required=True)
    mapping.add_argument("--output", required=True, type=Path)

    evidence = subparsers.add_parser("evidence")
    evidence.add_argument("--gaf", required=True, type=Path)
    evidence.add_argument("--target-table", required=True, type=Path)
    evidence.add_argument("--inventory", required=True, type=Path)
    evidence.add_argument("--sample", required=True)
    evidence.add_argument("--sample-group", required=True)
    evidence.add_argument("--batch", required=True)
    evidence.add_argument("--target", required=True)
    evidence.add_argument("--mapper-version", required=True)
    evidence.add_argument("--min-mapq", type=int, default=5)
    evidence.add_argument("--min-junction-reads", type=int, default=3)
    evidence.add_argument("--min-reference-reads", type=int, default=3)
    evidence.add_argument("--output", required=True, type=Path)

    catalog = subparsers.add_parser("catalog")
    catalog.add_argument("--bubble", action="append", default=[])
    catalog.add_argument("--reference-sample", default="CHM13")
    catalog.add_argument("--output", required=True, type=Path)

    merge = subparsers.add_parser("merge")
    merge.add_argument("--input", action="append", default=[])
    merge.add_argument("--output", required=True, type=Path)

    compare = subparsers.add_parser("compare")
    compare.add_argument("--metrics", required=True, type=Path)
    compare.add_argument("--output", required=True, type=Path)
    compare.add_argument("--permutations", type=int, default=10000)
    compare.add_argument("--seed", type=int, default=20260902)

    report = subparsers.add_parser("report")
    report.add_argument("--samples", required=True, type=Path)
    report.add_argument("--inventory", required=True, type=Path)
    report.add_argument("--membership", required=True, type=Path)
    report.add_argument("--catalog", required=True, type=Path)
    report.add_argument("--metrics", required=True, type=Path)
    report.add_argument("--evidence", required=True, type=Path)
    report.add_argument("--comparison", required=True, type=Path)
    report.add_argument("--targets", required=True, type=Path)
    report.add_argument("--validation", required=True, type=Path)
    report.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def number(value: str | None, default: float = 0.0) -> float:
    try:
        return float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return default


def ratio(numerator: float, denominator: float) -> float | str:
    return numerator / denominator if denominator else ""


def parse_fastq_stats(path: Path) -> tuple[int, int]:
    rows = read_tsv(path)
    reads = int(sum(number(row.get("num_seqs")) for row in rows))
    bases = int(sum(number(row.get("sum_len")) for row in rows))
    return reads, bases


def gaf_path_tokens(path: str) -> list[str]:
    if not path.startswith((">", "<")):
        return []
    return [direction + node for direction, node in re.findall(r"([<>])([^<>]+)", path)]


def gaf_tags(fields: list[str]) -> dict[str, str]:
    result = {}
    for field in fields[12:]:
        parts = field.split(":", 2)
        if len(parts) == 3:
            result[parts[0]] = parts[2]
    return result


def normalized_read_name(name: str) -> str:
    return re.sub(r"/(?:1|2)$", "", name.split()[0])


def histogram_median(histogram: Counter) -> float | str:
    count = sum(histogram.values())
    if not count:
        return ""
    ranks = ((count - 1) // 2, count // 2)
    values = []
    cumulative = 0
    for value, frequency in sorted(histogram.items()):
        values.extend(value for rank in ranks if cumulative <= rank < cumulative + frequency)
        cumulative += frequency
    return statistics.mean(values)


def mapping_from_gaf(path: Path, input_reads: int, input_bases: int) -> dict[str, object]:
    primary = secondary = aligned_bases = clipped_bases = 0
    identity_sum = mapq_sum = 0.0
    histogram = Counter()
    with open_text(path) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            try:
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                matches = int(fields[9])
                block_length = int(fields[10])
                mapq = int(fields[11])
            except ValueError:
                continue
            tags = gaf_tags(fields)
            if tags.get("tp") == "S":
                secondary += 1
                continue
            if fields[5] == "*" or block_length == 0:
                continue
            primary += 1
            if 0 <= mapq < 255:
                histogram[mapq] += 1
                mapq_sum += mapq
            identity_sum += matches / block_length
            aligned_bases += max(query_end - query_start, 0)
            clipped_bases += max(query_start, 0) + max(query_length - query_end, 0)
    return {
        "primary_mapped_reads": primary,
        "mapped_read_fraction": ratio(primary, input_reads),
        "mapq_0_fraction": ratio(histogram[0], input_reads),
        "mapq_ge_5_fraction": ratio(sum(count for value, count in histogram.items() if value >= 5), input_reads),
        "mapq_ge_20_fraction": ratio(sum(count for value, count in histogram.items() if value >= 20), input_reads),
        "mapq_ge_40_fraction": ratio(sum(count for value, count in histogram.items() if value >= 40), input_reads),
        "mean_mapq": ratio(mapq_sum, sum(histogram.values())),
        "median_mapq": histogram_median(histogram),
        "mean_alignment_identity": ratio(identity_sum, primary),
        "fraction_read_bases_aligned": ratio(aligned_bases, input_bases),
        "clipped_bases_per_read": ratio(clipped_bases, primary),
        "secondary_multimapping_fraction": ratio(secondary, input_reads),
    }


def parse_flagstat(path: Path) -> tuple[float | str, float | str]:
    mapped = ""
    proper = ""
    for line in path.read_text().splitlines():
        match = re.search(r"\(([0-9.]+)%", line)
        if not match:
            continue
        if re.match(r"^\d+ \+ \d+ mapped \(", line):
            mapped = float(match.group(1)) / 100
        elif "properly paired" in line:
            proper = float(match.group(1)) / 100
    return mapped, proper


def parse_coverage(path: Path) -> tuple[float | str, float | str]:
    length = covered = depth_sum = 0.0
    with path.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            region_length = number(row.get("endpos")) - number(row.get("startpos")) + 1
            if region_length <= 0:
                continue
            length += region_length
            covered += number(row.get("covbases"))
            depth_sum += number(row.get("meandepth")) * region_length
    return ratio(covered, length), ratio(depth_sum, length)


def benchmark_values(path: Path) -> tuple[float | str, float | str, float | str]:
    rows = read_tsv(path)
    if not rows:
        return "", "", ""
    row = rows[0]
    wall = row.get("s") or row.get("runtime") or ""
    rss = row.get("max_rss") or row.get("max_rss_mb") or ""
    cpu = row.get("cpu_time") or ""
    return wall, rss, cpu


def inventory_row(path: Path, target: str) -> dict[str, str]:
    rows = read_tsv(path)
    matches = [row for row in rows if row.get("target") == target]
    if len(matches) != 1:
        raise ValueError(f"expected one inventory row for {target} in {path}")
    return matches[0]


def mapping_stats(args: argparse.Namespace) -> None:
    input_reads, input_bases = parse_fastq_stats(args.fastq_stats)
    inventory = inventory_row(args.inventory, args.target)
    breadth, depth = parse_coverage(args.coverage)
    backbone_mapped, proper_pair = parse_flagstat(args.flagstat)
    wall, rss, cpu = benchmark_values(args.benchmark)
    index_wall, index_rss, _ = benchmark_values(args.index_benchmark)
    row = {
        "sample_id": args.sample,
        "sample_group": args.sample_group,
        "batch": args.batch,
        "target": args.target,
        "source_graph": inventory["source_graph"],
        "graph_checksum": inventory["graph_checksum"],
        "mapper_version": args.mapper_version,
        "input_reads": input_reads,
        "input_read_pairs": input_reads // 2,
        "input_bases": input_bases,
        **mapping_from_gaf(args.gaf, input_reads, input_bases),
        "proper_pair_fraction": proper_pair,
        "common_backbone_mapped_fraction": backbone_mapped,
        "common_backbone_breadth_fraction": breadth,
        "common_backbone_mean_depth": depth,
        "mapping_mode": "direct_full_graph",
        "mapping_graph_size_bytes": args.mapping_graph.stat().st_size,
        "mapping_index_size_bytes": sum(Path(path).stat().st_size for path in args.index),
        "index_wall_time_seconds": index_wall,
        "index_peak_rss_mb": index_rss,
        "mapping_wall_time_seconds": wall,
        "mapping_cpu_time_seconds": cpu,
        "mapping_peak_rss_mb": rss,
        "gaf_size_bytes": args.gaf.stat().st_size,
    }
    write_tsv(args.output, [row], MAPPING_FIELDS)


def parse_edges(value: str) -> set[tuple[str, str]]:
    edges = set()
    for item in value.split(";"):
        tokens = [token.strip() for token in item.split(",") if token.strip()]
        if tokens and (len(tokens) != 2 or any(not re.fullmatch(r"[<>][1-9]\d*", token) for token in tokens)):
            raise ValueError(f"Expected an oriented edge in vg node IDs, got: {item}")
        if len(tokens) == 2:
            edge = (tokens[0], tokens[1])
            edges.add(edge)
            reverse = tuple(
                ("<" if token.startswith(">") else ">") + token[1:]
                for token in reversed(edge)
            )
            edges.add(reverse)
    return edges


def evidence_stats(args: argparse.Namespace) -> None:
    targets = read_tsv(args.target_table)
    if len(targets) != 1:
        raise ValueError(f"expected one target row in {args.target_table}")
    target = targets[0]
    if target["graph_id"] != args.target:
        raise ValueError("TMEM target table belongs to a different graph")
    if not 0 <= args.min_mapq < 255 or min(args.min_junction_reads, args.min_reference_reads) < 1:
        raise ValueError("MAPQ must be in 0..254 and fragment thresholds must be positive")
    callable_status = target.get("callable_status", "not_callable")
    inventory = inventory_row(args.inventory, args.target)
    row = {
        "sample_id": args.sample, "sample_group": args.sample_group, "batch": args.batch,
        "target": args.target, "source_graph": inventory["source_graph"],
        "graph_checksum": inventory["graph_checksum"], "mapper_version": args.mapper_version,
        "sv_id": "TMEM106B_approximately_19kb_insertion",
        "target_present": target.get("target_present", "unknown"),
        "target_callable_status": callable_status,
        "evidence_available": "false", "counting_unit": "read_name_collapsed_fragments",
        "min_mapq": args.min_mapq, "min_junction_fragments": args.min_junction_reads,
        "min_reference_fragments": args.min_reference_reads,
        "ambiguous_or_conflicting": target.get("ambiguity", ""),
    }
    if not callable_status.startswith("callable"):
        row.update(call_class="NO_DECISION", note=target.get("failure_reason") or callable_status)
        # Blank counts mean unavailable evidence. Do not scan a genome to manufacture zeros.
        write_tsv(args.output, [row], EVIDENCE_FIELDS)
        return
    if target.get("edge_namespace") != "vg_node_ids":
        raise ValueError("TMEM edges must be translated using the GBZ that produced the read GAF")
    left_edges = parse_edges(target.get("left_junction", ""))
    right_edges = parse_edges(target.get("right_junction", ""))
    reference_edges = parse_edges(target.get("reference_defining_edges", ""))
    insertion_edges = parse_edges(target.get("insertion_edges", ""))
    insertion_nodes = set(filter(None, target.get("insertion_nodes", "").split(",")))
    if not left_edges or not right_edges or not reference_edges:
        raise ValueError("Callable target is missing allele-defining edges")
    if any(not re.fullmatch(r"[1-9]\d*", node) for node in insertion_nodes):
        raise ValueError("Insertion nodes must use numeric vg IDs")
    locus_nodes = insertion_nodes | {
        token[1:] for edge in left_edges | right_edges | reference_edges | insertion_edges for token in edge
    }
    left_reads: set[str] = set()
    right_reads: set[str] = set()
    reference_reads: set[str] = set()
    insertion_reads: set[str] = set()
    locus_fragments: set[str] = set()
    qualified_alignments = 0
    namespace_checked = False
    with open_text(args.gaf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[5] == "*":
                continue
            try:
                mapq = int(fields[11])
                aligned_length = int(fields[10])
            except ValueError:
                continue
            if aligned_length <= 0 or not args.min_mapq <= mapq < 255 or gaf_tags(fields).get("tp") == "S":
                continue
            qualified_alignments += 1
            tokens = gaf_path_tokens(fields[5])
            if not tokens:
                raise ValueError("Read GAF must contain vg node traversals, not named stable paths")
            if not namespace_checked:
                if any(not re.fullmatch(r"[<>][1-9]\d*", token) for token in tokens):
                    raise ValueError("Read GAF does not use numeric vg node IDs")
                namespace_checked = True
            nodes = {token[1:] for token in tokens}
            if not nodes & locus_nodes:
                continue
            name = normalized_read_name(fields[0])
            locus_fragments.add(name)
            edges = set(zip(tokens, tokens[1:]))
            if edges & left_edges:
                left_reads.add(name)
            if edges & right_edges:
                right_reads.add(name)
            if edges & reference_edges:
                reference_reads.add(name)
            if edges & insertion_edges or nodes & insertion_nodes:
                insertion_reads.add(name)

    conflicting = (left_reads | right_reads) & reference_reads
    left_reads.difference_update(conflicting)
    right_reads.difference_update(conflicting)
    reference_reads.difference_update(conflicting)
    left = len(left_reads)
    right = len(right_reads)
    reference = len(reference_reads)
    minimum_alt = min(left, right)
    both_alt = minimum_alt >= args.min_junction_reads
    ref_supported = reference >= args.min_reference_reads
    if len(conflicting) >= min(args.min_junction_reads, args.min_reference_reads):
        call = "CONFLICTING"
        note = "multiple fragments support incompatible REF and ALT edges; excluded from allele counts"
    elif both_alt and ref_supported:
        call = "HET_LIKE"
        note = "both ALT junctions and REF-specific edges have independent fragment support; not a formal genotype"
    elif both_alt:
        call = "ALT_SUPPORTED"
        note = "both ALT junctions meet the fragment threshold; zygosity is undetermined"
    elif ref_supported and left < args.min_junction_reads and right < args.min_junction_reads:
        call = "REF_SUPPORTED"
        note = "REF-specific edges meet the threshold; this alone does not exclude an ALT allele"
    elif max(left, right) >= args.min_junction_reads:
        call = "CONFLICTING"
        note = "only one ALT junction meets the configured threshold"
    else:
        call = "NO_DECISION"
        note = ("insufficient MAPQ-qualified junction and reference evidence" if locus_fragments else
                "no MAPQ-qualified fragments on target nodes; check locus coverage and graph provenance")

    row.update({
        "call_class": call,
        "evidence_available": "true",
        "left_alt_junction_reads": left,
        "right_alt_junction_reads": right,
        "minimum_alt_junction_reads": minimum_alt,
        "reference_path_reads": reference,
        "insertion_sequence_reads": len(insertion_reads),
        "locus_mapq_qualified_fragments": len(locus_fragments),
        "genome_mapq_qualified_alignments": qualified_alignments,
        "conflicting_fragments": len(conflicting),
        "alt_junction_balance": ratio(minimum_alt, max(left, right)),
        "note": note,
    })
    write_tsv(args.output, [row], EVIDENCE_FIELDS)


def canonical_reference_path(path: str, sample: str = "CHM13") -> str:
    for prefix in (f"{sample}#0#", f"{sample}."):
        if path.startswith(prefix):
            return path[len(prefix):]
    return path


def build_catalog(args: argparse.Namespace) -> None:
    fields = [
        "target", "sv_id", "reference_path", "start_0based", "end_0based",
        "reference_length", "sv_type", "size_difference", "size_bin", "segments", "graph_paths",
        "bidirectional", "minimum_allele_length", "maximum_allele_length",
        "segment_ids", "canonical_reference_path", "coordinate_key", "shared_graph_count",
        "shared_graphs", "representation_class", "match_basis", "allele_equivalence",
    ]
    rows = []
    for specification in args.bubble:
        target, filename = specification.split("=", 1)
        with Path(filename).open() as handle:
            for line_number, line in enumerate(handle, start=1):
                values = line.rstrip("\n").split("\t")
                if len(values) < 12:
                    continue
                minimum = int(values[6])
                maximum = int(values[7])
                reference_length = int(values[2]) - int(values[1])
                difference = max(abs(minimum - reference_length), abs(maximum - reference_length))
                bidirectional = values[5] == "1"
                if difference < 50 and not bidirectional:
                    continue
                if bidirectional:
                    sv_type = "INV_OR_COMPLEX"
                elif maximum >= reference_length + 50 and minimum >= reference_length - 49:
                    sv_type = "INS"
                elif minimum <= reference_length - 50 and maximum <= reference_length + 49:
                    sv_type = "DEL"
                else:
                    sv_type = "MULTIALLELIC_OR_COMPLEX"
                reference_path = canonical_reference_path(values[0], args.reference_sample)
                coordinate_key = f"{reference_path}:{values[1]}-{values[2]}:{sv_type}"
                if difference < 50:
                    size_bin = "complex_under_50bp_length_difference"
                elif difference < 100:
                    size_bin = "50-99bp"
                elif difference < 1000:
                    size_bin = "100-999bp"
                elif difference < 10000:
                    size_bin = "1-9kb"
                else:
                    size_bin = "10kb_or_larger"
                rows.append({
                    "target": target,
                    "sv_id": f"{target}:{values[0]}:{values[1]}-{values[2]}:{line_number}",
                    "reference_path": values[0],
                    "canonical_reference_path": reference_path,
                    "start_0based": values[1],
                    "end_0based": values[2],
                    "reference_length": reference_length,
                    "sv_type": sv_type,
                    "size_difference": difference,
                    "size_bin": size_bin,
                    "segments": values[3],
                    "graph_paths": values[4],
                    "bidirectional": str(bidirectional).lower(),
                    "minimum_allele_length": minimum,
                    "maximum_allele_length": maximum,
                    "segment_ids": values[11],
                    "coordinate_key": coordinate_key,
                    "match_basis": "exact_normalized_coordinates_and_sv_type",
                    "allele_equivalence": "not_assessed",
                })
    graph_counts = defaultdict(set)
    for row in rows:
        graph_counts[row["coordinate_key"]].add(row["target"])
    for row in rows:
        targets = graph_counts[row["coordinate_key"]]
        count = len(targets)
        row["shared_graph_count"] = count
        row["shared_graphs"] = ",".join(sorted(targets))
        if count == 4:
            representation = "shared_all_four"
        elif count == 1 and row["target"] == "hprc":
            representation = "hprc_only"
        elif count == 1 and row["target"] in {"local_v1", "local_v2"}:
            representation = "local_only"
        elif count == 1:
            representation = "combined_only"
        elif "hprc" in targets:
            representation = "shared_with_hprc"
        else:
            representation = "shared_without_hprc"
        row["representation_class"] = representation
    write_tsv(args.output, rows, fields)


def merge_tables(args: argparse.Namespace) -> None:
    if not args.input:
        raise ValueError("at least one --input is required")
    all_rows = []
    fields: list[str] | None = None
    for filename in args.input:
        with Path(filename).open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            current_fields = list(reader.fieldnames or [])
            if fields is None:
                fields = current_fields
            elif current_fields != fields:
                raise ValueError(f"table header differs: {filename}")
            all_rows.extend(reader)
    write_tsv(args.output, all_rows, fields or [])


def sign_flip_pvalue(values: list[float], permutations: int, rng: random.Random) -> float | str:
    if not values:
        return ""
    observed = abs(statistics.mean(values))
    if len(values) <= 16:
        total = 2 ** len(values)
        extreme = 0
        for mask in range(total):
            test = [value if mask & (1 << index) else -value for index, value in enumerate(values)]
            extreme += abs(statistics.mean(test)) >= observed
        return extreme / total
    extreme = 0
    for _ in range(permutations):
        test = [value if rng.random() < 0.5 else -value for value in values]
        extreme += abs(statistics.mean(test)) >= observed
    return (extreme + 1) / (permutations + 1)


def paired_comparison(args: argparse.Namespace) -> None:
    rows = read_tsv(args.metrics)
    by_sample = defaultdict(dict)
    groups = {}
    batches = {}
    input_bases = {}
    input_reads = {}
    for row in rows:
        sample = row["sample_id"]
        if row["target"] in by_sample[sample]:
            raise ValueError(f"Duplicate mapping result: {sample}/{row['target']}")
        if sample in groups and (
            groups[sample] != row["sample_group"] or batches[sample] != row["batch"]
            or input_bases[sample] != number(row.get("input_bases"))
            or input_reads[sample] != number(row.get("input_reads"))
        ):
            raise ValueError(f"Sample metadata or read input differs between graphs: {sample}")
        by_sample[row["sample_id"]][row["target"]] = row
        groups[row["sample_id"]] = row["sample_group"]
        batches[row["sample_id"]] = row["batch"]
        input_bases[row["sample_id"]] = number(row.get("input_bases"))
        input_reads[sample] = number(row.get("input_reads"))
    expected_targets = {target for _, left, right in CONTRASTS for target in (left, right)}
    if not rows:
        raise ValueError("No mapping results to compare")
    for sample, results in by_sample.items():
        if set(results) != expected_targets:
            raise ValueError(f"Incomplete four-graph comparison for {sample}: {sorted(results)}")
    ranked_samples = sorted(input_bases, key=lambda sample: (input_bases[sample], sample))
    depth_bins = {
        sample: ("low" if rank * 3 < len(ranked_samples) else
                 "medium" if rank * 3 < 2 * len(ranked_samples) else "high")
        for rank, sample in enumerate(ranked_samples)
    }
    output_rows = []
    rng = random.Random(args.seed)
    for contrast, left, right in CONTRASTS:
        for metric in COMPARE_METRICS:
            sample_differences = []
            for sample in sorted(by_sample):
                if left not in by_sample[sample] or right not in by_sample[sample]:
                    continue
                left_value = number(by_sample[sample][left].get(metric), math.nan)
                right_value = number(by_sample[sample][right].get(metric), math.nan)
                if not math.isfinite(left_value) or not math.isfinite(right_value):
                    continue
                difference = left_value - right_value
                sample_differences.append((sample, difference))
                output_rows.append({
                    "level": "sample",
                    "sample_id": sample,
                    "sample_group": groups[sample],
                    "batch": batches[sample],
                    "read_depth_bin": depth_bins[sample],
                    "stratum_type": "",
                    "stratum": "",
                    "contrast": contrast,
                    "left_target": left,
                    "right_target": right,
                    "metric": metric,
                    "left_value": left_value,
                    "right_value": right_value,
                    "difference": difference,
                    "difference_percentage_points": 100 * difference if metric in FRACTION_METRICS else "",
                    "read_count_difference": round(difference * input_reads[sample])
                    if metric == "mapped_read_fraction" or metric.startswith("mapq_ge_") else "",
                    "n_pairs": "",
                    "mean_difference": "",
                    "median_difference": "",
                    "sign_flip_pvalue": "",
                })

            strata = [("all", "all", [value for _, value in sample_differences])]
            for stratum_type, labels in [
                ("batch", batches),
                ("sample_group", groups),
                ("read_depth_bin", depth_bins),
            ]:
                for stratum in sorted(set(labels.values())):
                    values = [value for sample, value in sample_differences if labels[sample] == stratum]
                    strata.append((stratum_type, stratum, values))

            for stratum_type, stratum, differences in strata:
                output_rows.append({
                    "level": "summary",
                    "sample_id": "",
                    "sample_group": stratum if stratum_type == "sample_group" else "all",
                    "batch": stratum if stratum_type == "batch" else "all",
                    "read_depth_bin": stratum if stratum_type == "read_depth_bin" else "all",
                    "stratum_type": stratum_type,
                    "stratum": stratum,
                    "contrast": contrast,
                    "left_target": left,
                    "right_target": right,
                    "metric": metric,
                    "left_value": "",
                    "right_value": "",
                    "difference": "",
                    "n_pairs": len(differences),
                    "mean_difference": statistics.mean(differences) if differences else "",
                    "median_difference": statistics.median(differences) if differences else "",
                    "mean_difference_percentage_points": 100 * statistics.mean(differences)
                    if differences and metric in FRACTION_METRICS else "",
                    "min_difference": min(differences) if differences else "",
                    "max_difference": max(differences) if differences else "",
                    "n_positive": sum(value > 0 for value in differences),
                    "n_negative": sum(value < 0 for value in differences),
                    "n_tied": sum(value == 0 for value in differences),
                    "sign_flip_pvalue": sign_flip_pvalue(differences, args.permutations, rng),
                })
    # Adjust the overall mapping comparisons together; strata remain exploratory.
    family = sorted(
        [row for row in output_rows if row.get("stratum_type") == "all" and row["sign_flip_pvalue"] != ""],
        key=lambda row: row["sign_flip_pvalue"],
    )
    adjusted = 1.0
    for rank in range(len(family), 0, -1):
        row = family[rank - 1]
        adjusted = min(adjusted, row["sign_flip_pvalue"] * len(family) / rank)
        row["pvalue_bh_overall"] = adjusted
    fields = [
        "level", "sample_id", "sample_group", "batch", "read_depth_bin",
        "stratum_type", "stratum", "contrast", "left_target", "right_target",
        "metric", "left_value", "right_value", "difference",
        "difference_percentage_points", "read_count_difference",
        "n_pairs", "mean_difference", "median_difference", "mean_difference_percentage_points",
        "min_difference", "max_difference", "n_positive", "n_negative", "n_tied",
        "sign_flip_pvalue", "pvalue_bh_overall",
    ]
    write_tsv(args.output, output_rows, fields)


def mean_by_target(rows: list[dict[str, str]], field: str) -> dict[str, float]:
    values = defaultdict(list)
    for row in rows:
        value = number(row.get(field), math.nan)
        if math.isfinite(value):
            values[row["target"]].append(value)
    return {target: statistics.mean(target_values) for target, target_values in values.items()}


def render_report(args: argparse.Namespace) -> None:
    samples = read_tsv(args.samples)
    inventory = read_tsv(args.inventory)
    membership = read_tsv(args.membership)
    catalog = read_tsv(args.catalog)
    metrics = read_tsv(args.metrics)
    evidence = read_tsv(args.evidence)
    comparisons = read_tsv(args.comparison)
    graph_targets = read_tsv(args.targets)
    validation = read_tsv(args.validation)
    targets = sorted({row["target"] for row in inventory})
    sample_ids = {row["sample_id"] for row in samples}
    if not sample_ids or len(sample_ids) != len(samples) or len(targets) != 4 or len(inventory) != 4:
        raise ValueError("Report requires distinct pilot samples and four distinct graph inventories")
    expected = {(sample, target) for sample in sample_ids for target in targets}
    sample_metadata = {row["sample_id"]: row for row in samples}
    checksums = {row["target"]: row["graph_checksum"] for row in inventory}
    for label, rows in (("mapping", metrics), ("evidence", evidence), ("membership", membership)):
        keys = [(row["sample_id"], row["target"]) for row in rows]
        if len(keys) != len(expected) or set(keys) != expected:
            raise ValueError(f"Incomplete or duplicated {label} sample/graph results")
        for row in rows:
            if row["sample_group"] != sample_metadata[row["sample_id"]]["sample_group"]:
                raise ValueError(f"Sample group mismatch in {label}: {row['sample_id']}")
            if label != "membership" and row["batch"] != sample_metadata[row["sample_id"]]["batch"]:
                raise ValueError(f"Sample batch mismatch in {label}: {row['sample_id']}")
            if label != "membership" and row["graph_checksum"] != checksums[row["target"]]:
                raise ValueError(f"Graph checksum mismatch in {label}: {row['target']}")
    if (len(validation) != len(samples)
            or {row["sample_id"] for row in validation} != sample_ids
            or any(row["runtime_validation"] != "PASS" for row in validation)):
        raise ValueError("TMEM prior-group validation must contain exactly the selected samples, all PASS")
    if len(graph_targets) != 4 or {row["graph_id"] for row in graph_targets} != set(targets):
        raise ValueError("Report requires one TMEM target audit per graph")
    if any(row["is_graph_member"] != "false" or row["membership_status"] != "known" for row in membership):
        raise ValueError("Cannot compare held-out samples with positive or unknown graph membership")
    thresholds = sorted({(row["min_mapq"], row["min_junction_fragments"], row["min_reference_fragments"]) for row in evidence})
    if len(thresholds) != 1:
        raise ValueError("Cannot compare evidence tables generated with different support thresholds")

    def percent(value) -> str:
        value = number(value, math.nan)
        return f"{100 * value:.3f}%" if math.isfinite(value) else "N/A"

    def display(value, digits=3) -> str:
        value = number(value, math.nan)
        return f"{value:.{digits}f}" if math.isfinite(value) else "N/A"

    mapped = mean_by_target(metrics, "mapped_read_fraction")
    mapq20 = mean_by_target(metrics, "mapq_ge_20_fraction")
    breadth = mean_by_target(metrics, "common_backbone_breadth_fraction")
    identity = mean_by_target(metrics, "mean_alignment_identity")
    callable_counts = Counter(row["target"] for row in evidence if row["call_class"] in DECISIVE_CALLS)
    catalog_counts = Counter(row["target"] for row in catalog)
    best_mapping = ", ".join(
        target for target in targets if target in mapq20 and max(mapq20.values()) - mapq20[target] <= 0.00001
    ) or "not available"
    lines = [
        "# Direct SV graph comparison", "",
        f"Pilot samples: **{len(samples)}**. Source graphs: **{len(inventory)}**.", "",
        "Each sample was mapped directly to the complete source graph. Pathless rGFA inputs were indexed "
        "with the greedy path cover generated by `vg autoindex`; no sample-specific haplotype selection was used.", "",
        "## Main result", "",
        f"- Highest mean MAPQ >=20 fraction: **{best_mapping}** (ties within 0.001 percentage points).",
        f"- TMEM sample/graph pairs with usable allele evidence: **{sum(callable_counts.values())}/{len(evidence)}**. "
        "`CONFLICTING` and `NO_DECISION` are excluded from usable decisions.",
        "- Higher mapping confidence is not a measurement of mapping or genotyping accuracy.",
        "- General graph-wide SV genotyping is not inferred from bubble inventory alone; "
        "the read-evidence table in this pilot is restricted to the approximately 19-kb TMEM insertion.",
        "", "## Per-target summary", "",
        "Unweighted means across the same samples. Mapped/MAPQ fractions use all input reads as denominator.", "",
        "| Target | Mapped | MAPQ >=20 | Alignment identity | CHM13 breadth | TMEM usable / total | Catalogue entries |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for target in targets:
        lines.append(
            f"| {target} | {percent(mapped.get(target))} | {percent(mapq20.get(target))} | "
            f"{percent(identity.get(target))} | {percent(breadth.get(target))} | "
            f"{callable_counts[target]}/{len(samples)} | {catalog_counts[target]} |"
        )
    lines.extend([
        "", "## Paired mapping effects", "",
        "Differences are left minus right, in percentage points (pp). The sample is the paired unit. "
        "Positive/negative counts describe direction, not truth-set accuracy. BH adjustment covers all overall "
        "mapping contrasts and metrics; batch, prior-group and input-depth strata remain exploratory.", "",
        "| Contrast | Metric | Pairs | Mean pp | Median pp | Min pp | Max pp | Positive / negative / tied | Sign-flip p | BH p |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in comparisons:
        if row.get("level") != "summary" or row.get("stratum_type") != "all":
            continue
        if row["metric"] not in {"mapped_read_fraction", "mapq_ge_20_fraction", "common_backbone_breadth_fraction"}:
            continue
        values = [display(100 * number(row.get(field), math.nan), 5)
                  for field in ("mean_difference", "median_difference", "min_difference", "max_difference")]
        lines.append(
            f"| {row['contrast']} | {row['metric']} | {row['n_pairs']} | " + " | ".join(values)
            + f" | {row['n_positive']}/{row['n_negative']}/{row['n_tied']} | "
            + f"{display(row.get('sign_flip_pvalue'), 6)} | {display(row.get('pvalue_bh_overall'), 6)} |"
        )
    lines.extend([
        "", "## TMEM graph callability", "",
        "A failed query audit does not establish allele absence. Traversals are source GFA segments; "
        "allele-defining edges and insertion nodes are translated to the existing mapping GBZ node IDs.", "",
        "| Graph | Status | REF coverage / identity | ALT coverage / identity | Edge namespace | Failure reason |",
        "|---|---|---|---|---|---|",
    ])
    for row in graph_targets:
        lines.append(
            f"| {row['graph_id']} | {row['callable_status']} | "
            f"{percent(row['reference_query_coverage'])} / {percent(row['reference_query_identity'])} | "
            f"{percent(row['alternate_query_coverage'])} / {percent(row['alternate_query_identity'])} | "
            f"{row['edge_namespace']} | {row['failure_reason'] or 'none'} |"
        )
    if not sum(callable_counts.values()):
        lines.extend(["", "**TMEM comparison unavailable: no usable sample-level allele decisions. "
                      "No graph can be ranked for TMEM discrimination.**"])
    lines.extend([
        "", "## TMEM evidence by prior group", "",
        "ALT-like combines `ALT_SUPPORTED` and `HET_LIKE`. Rates exclude unresolved/conflicting results; "
        "the explicit counts show the resulting denominator. These selected prior groups are not truth "
        "genotypes or disease/control groups, so these are not sensitivity, specificity or population frequencies.", "",
        "| Graph | Prior group | Total | ALT-like | REF-supported | Conflicting | No decision | ALT-like / usable |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ])
    for target in targets:
        for group in sorted({row["sample_group"] for row in samples}):
            rows = [row for row in evidence if row["target"] == target and row["sample_group"] == group]
            counts = Counter(row["call_class"] for row in rows)
            alt = sum(counts[call] for call in ALT_CALLS)
            usable = alt + counts["REF_SUPPORTED"]
            lines.append(
                f"| {target} | {group} | {len(rows)} | {alt} | {counts['REF_SUPPORTED']} | "
                f"{counts['CONFLICTING']} | {counts['NO_DECISION']} | {percent(ratio(alt, usable))} |"
            )
    evidence_by_pair = {(row["sample_id"], row["target"]): row for row in evidence}
    lines.extend([
        "", "## Paired TMEM agreement", "",
        "Comparison is restricted to samples with usable evidence in both graphs. Agreement compares "
        "ALT-like versus REF-supported evidence, not zygosity; agreement is not accuracy.", "",
        "| Contrast | Jointly usable | Both ALT-like | Both REF-supported | Left ALT / right REF | Left REF / right ALT | Excluded |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ])
    for contrast, left, right in CONTRASTS:
        counts = Counter()
        for sample in sorted(sample_ids):
            lcall = evidence_by_pair[sample, left]["call_class"]
            rcall = evidence_by_pair[sample, right]["call_class"]
            if lcall in DECISIVE_CALLS and rcall in DECISIVE_CALLS:
                counts[lcall in ALT_CALLS, rcall in ALT_CALLS] += 1
        usable = sum(counts.values())
        lines.append(
            f"| {contrast} | {usable} | {counts[True, True]} | {counts[False, False]} | "
            f"{counts[True, False]} | {counts[False, True]} | {len(samples) - usable} |"
        )
    lines.extend([
        "", "## Per-sample TMEM support", "",
        "Cells show evidence class and left ALT / right ALT / REF fragment counts. "
        "Blank evidence is shown as N/A. Mate names are collapsed; PCR duplicates are not removed. "
        "Fragments supporting incompatible REF and ALT edges are excluded from allele counts and reported separately in the TSV.", "",
        "| Sample | " + " | ".join(targets) + " |",
        "|---|" + "---|" * len(targets),
    ])
    for sample in sorted(sample_ids):
        cells = []
        for target in targets:
            row = evidence_by_pair[sample, target]
            support = "/".join(row.get(field) or "N/A" for field in
                               ("left_alt_junction_reads", "right_alt_junction_reads", "reference_path_reads"))
            cells.append(f"{row['call_class']} ({support})")
        lines.append(f"| {sample} | " + " | ".join(cells) + " |")
    shared = {row["coordinate_key"] for row in catalog if row["shared_graph_count"] == "4"}
    lines.extend([
        "", "## Catalogue and provenance", "",
        f"- Normalized coordinate/type keys present in all four graphs: **{len(shared)}**.",
        "- Chromosome aliases are harmonized only for the configured common reference sample. "
        "Exact coordinates/type do not establish allele-sequence identity; unmatched keys do not establish novelty.",
        "- Known graph-member sample/target pairs: **0**. Relatedness to donors was not evaluated.",
        f"- Evidence thresholds (minimum MAPQ, fragments per ALT junction, REF fragments): `{thresholds}`.",
        f"- Per-sample mapping differences, strata and overall BH adjustments: `{args.comparison}`.",
        f"- Full TMEM fragment evidence: `{args.evidence}`. Graph target audits: `{args.targets}`.",
        "", "## Interpretation notes", "",
        "- Junction support, not insertion-sequence coverage alone, defines ALT support.",
        "- `NO_DECISION` is not evidence for absence.",
        "- `ALT_SUPPORTED` does not distinguish heterozygous from homozygous insertion carriage. "
        "`HET_LIKE` is supporting evidence, not a likelihood-based genotype.",
        "- `REF_SUPPORTED` establishes reference-allele evidence, not the absence of the insertion allele.",
        "- Bubble counts measure representation, not mapping improvement.",
        "- Greedy path-cover threads are indexing aids, not observed biological haplotypes.",
        "- The rGFA bubble catalog distinguishes insertions, deletions, bidirectional/complex events, and size bins; "
        "it cannot reliably distinguish tandem duplications from other insertion-like bubbles without extra annotation.",
        "- Read-depth strata are rank-based tertiles of input bases across the pilot samples.",
        "- Zero secondary alignments/clipping in Giraffe GAF are not independent proof of unique or perfect mapping. "
        "The query-span aligned fraction may be identical to the mapped-read fraction.",
        "- Mapping benchmarks include Giraffe, compression, surjection, BAM sorting and indexing; "
        "they are not isolated mapper runtimes.",
        "- Results can depend on graph construction, path-cover structure, read depth, and repetitive sequence.",
    ])
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(str(args.output) + ".tmp")
    temporary.write_text("\n".join(lines) + "\n")
    temporary.replace(args.output)


def main() -> None:
    args = parse_args()
    commands = {
        "mapping": mapping_stats,
        "evidence": evidence_stats,
        "catalog": build_catalog,
        "merge": merge_tables,
        "compare": paired_comparison,
        "report": render_report,
    }
    commands[args.command](args)


if __name__ == "__main__":
    main()
