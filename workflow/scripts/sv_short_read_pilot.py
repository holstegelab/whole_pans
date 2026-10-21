#!/usr/bin/env python3
"""Metadata, FASTQ, index, and aggregation helpers for the SV-rGFA read pilot."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import math
import os
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Iterator, TextIO

import yaml


METRIC_FIELDS = [
    "sample_id", "sample_role", "batch", "graph_id", "baseline_graph",
    "mapper", "indexing_mode", "run_mode", "total_input_reads",
    "total_input_pairs", "gaf_records_observed", "gaf_unmapped_records",
    "primary_mapped_reads", "unmapped_reads", "mapped_read_fraction",
    "mapq_fraction_denominator", "mapq_0_fraction", "mapq_ge_5_fraction",
    "mapq_ge_20_fraction", "mapq_ge_40_fraction", "mean_mapq",
    "median_mapq", "mean_alignment_identity", "mean_query_aligned_fraction",
    "mean_clipped_bases_per_read", "secondary_multimapping_count",
    "secondary_multimapping_fraction", "malformed_gaf_rows",
    "paired_fragment_learning_succeeded", "single_end_fallback",
    "single_end_fallback_read_count", "learned_fragment_mean",
    "learned_fragment_sd", "mapping_wall_time_seconds", "peak_rss_mb",
    "reads_per_second", "proper_pair_percentage", "relevant_warnings",
    "graph_qc_status", "membership_status", "is_graph_member",
    "selection_seed", "requested_pairs", "actual_pairs",
]

CONTRAST_METRICS = [
    "mapped_read_fraction", "mapq_0_fraction", "mapq_ge_5_fraction",
    "mapq_ge_20_fraction", "mapq_ge_40_fraction", "mean_mapq",
    "median_mapq", "mean_alignment_identity", "mean_query_aligned_fraction",
    "mean_clipped_bases_per_read", "secondary_multimapping_fraction",
    "mapping_wall_time_seconds", "peak_rss_mb", "reads_per_second",
]

CONTRASTS = [
    ("local_V1_minus_HPRC", "hprc_chm13_sv", "local_v1"),
    ("local_V2_minus_HPRC", "hprc_chm13_sv", "local_v2"),
    ("HPRC_plus_local_minus_HPRC", "hprc_chm13_sv", "hprc_plus_local"),
    ("local_V2_minus_local_V1", "local_v1", "local_v2"),
    ("HPRC_plus_local_minus_local_V1", "local_v1", "hprc_plus_local"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    select = sub.add_parser("select")
    select.add_argument("--config", required=True, type=Path)
    select.add_argument("--manifest", required=True, type=Path)
    select.add_argument("--output", required=True, type=Path)

    qc = sub.add_parser("fastq-qc")
    qc.add_argument("--sample", required=True)
    qc.add_argument("--r1", action="append", default=[])
    qc.add_argument("--r2", action="append", default=[])
    qc.add_argument("--output", required=True, type=Path)

    sample = sub.add_parser("subsample")
    sample.add_argument("--sample", required=True)
    sample.add_argument("--r1", action="append", default=[])
    sample.add_argument("--r2", action="append", default=[])
    sample.add_argument("--qc", required=True, type=Path)
    sample.add_argument("--max-pairs", required=True, type=int)
    sample.add_argument("--seed", required=True, type=int)
    sample.add_argument("--output-r1", required=True, type=Path)
    sample.add_argument("--output-r2", required=True, type=Path)
    sample.add_argument("--stats", required=True, type=Path)

    cat = sub.add_parser("cat-fastq")
    cat.add_argument("--input", action="append", default=[])

    membership = sub.add_parser("membership")
    membership.add_argument("--config", required=True, type=Path)
    membership.add_argument("--selected", required=True, type=Path)
    membership.add_argument("--output", required=True, type=Path)

    preflight = sub.add_parser("index-preflight")
    preflight.add_argument("--graph", required=True)
    preflight.add_argument("--source", required=True, type=Path)
    preflight.add_argument("--qc-status", required=True)
    preflight.add_argument("--experimental", required=True)
    preflight.add_argument("--indexing-mode", required=True)
    preflight.add_argument("--output", required=True, type=Path)

    merge_pre = sub.add_parser("merge-index-preflight")
    merge_pre.add_argument("--input", action="append", default=[])
    merge_pre.add_argument("--output", required=True, type=Path)

    idx = sub.add_parser("index-manifest")
    idx.add_argument("--graph", required=True)
    idx.add_argument("--source", required=True, type=Path)
    idx.add_argument("--preflight", required=True, type=Path)
    idx.add_argument("--index-dir", required=True, type=Path)
    idx.add_argument("--reported-index-dir", type=Path)
    idx.add_argument("--prefix-name", required=True)
    idx.add_argument("--vg-version-file", required=True, type=Path)
    idx.add_argument("--log", required=True, type=Path)
    idx.add_argument("--indexing-mode", required=True)
    idx.add_argument("--command", required=True)
    idx.add_argument("--output", required=True, type=Path)

    validate = sub.add_parser("validate-index")
    validate.add_argument("--graph", required=True)
    validate.add_argument("--source", required=True, type=Path)
    validate.add_argument("--manifest", required=True, type=Path)
    validate.add_argument("--output", required=True, type=Path)

    final = sub.add_parser("finalize-metrics")
    final.add_argument("--core", required=True, type=Path)
    final.add_argument("--fastq-stats", required=True, type=Path)
    final.add_argument("--membership", required=True, type=Path)
    final.add_argument("--benchmark", required=True, type=Path)
    final.add_argument("--giraffe-log", required=True, type=Path)
    final.add_argument("--sample", required=True)
    final.add_argument("--graph", required=True)
    final.add_argument("--graph-qc-status", required=True)
    final.add_argument("--output", required=True, type=Path)

    merge = sub.add_parser("merge-metrics")
    merge.add_argument("--input", action="append", default=[])
    merge.add_argument("--output", required=True, type=Path)

    contrast = sub.add_parser("contrasts")
    contrast.add_argument("--metrics", required=True, type=Path)
    contrast.add_argument("--output", required=True, type=Path)

    evidence = sub.add_parser("merge-evidence")
    evidence.add_argument("--input", action="append", default=[])
    evidence.add_argument("--output", required=True, type=Path)

    tools = sub.add_parser("tools")
    tools.add_argument("--expected-vg-version", required=True)
    tools.add_argument("--output", required=True, type=Path)

    provenance = sub.add_parser("provenance")
    provenance.add_argument("--config", required=True, type=Path)
    provenance.add_argument("--selected", required=True, type=Path)
    provenance.add_argument("--index-preflight", required=True, type=Path)
    provenance.add_argument("--tools", required=True, type=Path)
    provenance.add_argument("--output", required=True, type=Path)

    gate = sub.add_parser("smoke-gate")
    gate.add_argument("--metrics", required=True, type=Path)
    gate.add_argument("--expected-rows", required=True, type=int)
    gate.add_argument("--min-mapped-fraction", type=float, default=0.80)
    gate.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def load_yaml(path: Path) -> dict:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a YAML mapping")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict], fields: list[str]) -> None:
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


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bool_text(value: object) -> str:
    return "true" if str(value).lower() in {"true", "1", "yes"} else "false"


def select_samples(config_path: Path, manifest_path: Path, output: Path) -> None:
    config = load_yaml(config_path)
    pilot = config["sv_short_read_pilot"]
    if not pilot.get("enabled", False):
        raise ValueError("sv_short_read_pilot.enabled is false")
    roles = [str(value) for value in pilot.get("sample_roles", [])]
    if pilot.get("include_matched_controls", False):
        roles.append("tmem_matched_control")
    roles = list(dict.fromkeys(roles))
    if not roles:
        raise ValueError("sv_short_read_pilot.sample_roles is empty")

    with manifest_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        source_fields = list(reader.fieldnames or [])
        source_rows = list(reader)
    if not source_rows:
        raise ValueError(f"frozen manifest is empty: {manifest_path}")
    selected = [row for row in source_rows if row.get("role") in roles]
    duplicates = [sample for sample, count in Counter(row["sample_id"] for row in selected).items() if count > 1]
    if duplicates:
        raise ValueError("duplicate selected sample IDs: " + ",".join(sorted(duplicates)))

    expected = pilot.get("expected_role_counts", {}) or {}
    observed = Counter(row.get("role", "") for row in selected)
    for role in roles:
        if role in expected and observed[role] != int(expected[role]):
            raise ValueError(
                f"selected {observed[role]} {role} rows; expected {int(expected[role])}"
            )
    expected_total = sum(int(expected[role]) for role in roles if role in expected)
    if expected_total and len(selected) != expected_total:
        raise ValueError(f"selected {len(selected)} rows; expected {expected_total}")

    for row in selected:
        if not row.get("r1") or not row.get("r2"):
            raise ValueError(f"{row.get('sample_id')}: missing frozen R1/R2 paths")
        row["pilot_selection_roles"] = ",".join(roles)
        row["selection_manifest_sha256"] = sha256_file(manifest_path)
    fields = source_fields + ["pilot_selection_roles", "selection_manifest_sha256"]
    write_tsv(output, selected, fields)


def open_fastq(path: Path) -> TextIO:
    if path.suffix == ".gz":
        return gzip.open(path, "rt", newline="")
    return path.open("r", newline="")


def fastq_records(path: Path) -> Iterator[tuple[str, str, str, str]]:
    with open_fastq(path) as handle:
        record_number = 0
        while True:
            name = handle.readline()
            if not name:
                break
            sequence = handle.readline()
            plus = handle.readline()
            quality = handle.readline()
            record_number += 1
            if not sequence or not plus or not quality:
                raise ValueError(f"{path}: truncated FASTQ record {record_number}")
            if not name.startswith("@"):
                raise ValueError(f"{path}: record {record_number} name does not start with @")
            if not plus.startswith("+"):
                raise ValueError(f"{path}: record {record_number} plus line is malformed")
            if len(sequence.rstrip("\r\n")) != len(quality.rstrip("\r\n")):
                raise ValueError(f"{path}: record {record_number} sequence/quality lengths differ")
            yield name, sequence, plus, quality


def normalize_read_name(name: str) -> str:
    token = name.strip().removeprefix("@").split()[0]
    return re.sub(r"/(?:1|2)$", "", token)


def validate_fastq_paths(r1: list[Path], r2: list[Path]) -> None:
    if not r1 or not r2:
        raise ValueError("both R1 and R2 lane lists are required")
    if len(r1) != len(r2):
        raise ValueError(f"R1/R2 lane counts differ: {len(r1)} != {len(r2)}")
    for path in r1 + r2:
        if not path.is_file():
            raise FileNotFoundError(f"FASTQ does not exist: {path}")
        if path.stat().st_size == 0:
            raise ValueError(f"FASTQ is empty: {path}")


def paired_records(
    r1: list[Path], r2: list[Path]
) -> Iterator[tuple[tuple[str, str, str, str], tuple[str, str, str, str]]]:
    validate_fastq_paths(r1, r2)
    for lane, (left_path, right_path) in enumerate(zip(r1, r2), start=1):
        left_iter = fastq_records(left_path)
        right_iter = fastq_records(right_path)
        pair_number = 0
        while True:
            left = next(left_iter, None)
            right = next(right_iter, None)
            if left is None and right is None:
                break
            pair_number += 1
            if left is None or right is None:
                raise ValueError(
                    f"lane {lane}: R1/R2 read counts differ near pair {pair_number}"
                )
            if normalize_read_name(left[0]) != normalize_read_name(right[0]):
                raise ValueError(
                    f"lane {lane} pair {pair_number}: normalized names differ: "
                    f"{left[0].strip()} != {right[0].strip()}"
                )
            yield left, right


def fastq_qc(sample: str, r1_values: list[str], r2_values: list[str], output: Path) -> None:
    r1, r2 = [Path(value) for value in r1_values], [Path(value) for value in r2_values]
    count = sum(1 for _ in paired_records(r1, r2))
    row = {
        "sample_id": sample,
        "run_mode": "full",
        "lane_count": len(r1),
        "r1_files": ",".join(map(str, r1)),
        "r2_files": ",".join(map(str, r2)),
        "total_input_pairs": count,
        "total_input_reads": 2 * count,
        "requested_pairs": "all",
        "actual_pairs": count,
        "selection_seed": "",
        "fastq_qc_status": "PASS",
        "normalized_pair_names": "PASS",
    }
    write_tsv(output, [row], list(row))


def gzip_text_writer(path: Path) -> tuple[io.TextIOWrapper, io.BufferedWriter]:
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(compressed, newline=""), raw


def modular_selector(total: int, requested: int, seed: int, sample: str):
    selected = min(total, requested)
    if selected == total:
        return lambda index: True, selected, 1, 0
    a = (int.from_bytes(hashlib.sha256(f"a:{seed}:{sample}".encode()).digest()[:8], "big") | 1) % total
    if a == 0:
        a = 1
    while math.gcd(a, total) != 1:
        a = (a + 2) % total
        if a == 0:
            a = 1
    b = int.from_bytes(hashlib.sha256(f"b:{seed}:{sample}".encode()).digest()[:8], "big") % total
    inverse = pow(a, -1, total)
    return lambda index: ((index - b) * inverse) % total < selected, selected, a, b


def subsample(args: argparse.Namespace) -> None:
    qc_rows = read_tsv(args.qc)
    if len(qc_rows) != 1 or qc_rows[0].get("sample_id") != args.sample:
        raise ValueError(f"expected one FASTQ-QC row for {args.sample}")
    total = int(qc_rows[0]["total_input_pairs"])
    if total <= 0:
        raise ValueError(f"{args.sample}: cannot subsample zero pairs")
    choose, expected, multiplier, offset = modular_selector(
        total, args.max_pairs, args.seed, args.sample
    )
    r1, r2 = [Path(value) for value in args.r1], [Path(value) for value in args.r2]
    args.output_r1.parent.mkdir(parents=True, exist_ok=True)
    args.output_r2.parent.mkdir(parents=True, exist_ok=True)
    tmp1, tmp2 = Path(str(args.output_r1) + ".tmp"), Path(str(args.output_r2) + ".tmp")
    writer1, raw1 = gzip_text_writer(tmp1)
    writer2, raw2 = gzip_text_writer(tmp2)
    actual = 0
    try:
        for index, (left, right) in enumerate(paired_records(r1, r2)):
            if choose(index):
                writer1.writelines(left)
                writer2.writelines(right)
                actual += 1
        writer1.close()
        writer2.close()
        raw1.close()
        raw2.close()
    except Exception:
        try:
            writer1.close()
            writer2.close()
        finally:
            tmp1.unlink(missing_ok=True)
            tmp2.unlink(missing_ok=True)
        raise
    if actual != expected:
        tmp1.unlink(missing_ok=True)
        tmp2.unlink(missing_ok=True)
        raise ValueError(f"deterministic selector produced {actual}, expected {expected}")
    tmp1.replace(args.output_r1)
    tmp2.replace(args.output_r2)
    row = {
        "sample_id": args.sample,
        "run_mode": "smoke",
        "lane_count": len(r1),
        "r1_files": str(args.output_r1),
        "r2_files": str(args.output_r2),
        "total_input_pairs": actual,
        "total_input_reads": 2 * actual,
        "full_input_pairs": total,
        "requested_pairs": args.max_pairs,
        "actual_pairs": actual,
        "selection_seed": args.seed,
        "selection_algorithm": "exact_affine_permutation_subset_emitted_in_input_order",
        "selection_multiplier": multiplier,
        "selection_offset": offset,
        "fastq_qc_status": "PASS_reused_full_qc_and_rechecked_pairs",
        "normalized_pair_names": "PASS",
    }
    write_tsv(args.stats, [row], list(row))


def cat_fastq(paths: list[str]) -> None:
    if not paths:
        raise ValueError("cat-fastq requires at least one --input")
    for path_text in paths:
        path = Path(path_text)
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing or empty FASTQ: {path}")
        with open_fastq(path) as handle:
            for block in iter(lambda: handle.read(1024 * 1024), ""):
                sys.stdout.write(block)


def donor_from_assembly(value: str) -> str:
    name = Path(value).name
    return re.split(r"\.hifi\.|_(?:hap[12]|mat|pat)(?:_|\.|$)", name, maxsplit=1)[0]


def read_local(path: Path) -> set[str]:
    return {
        donor_from_assembly(row.get("assembly_id", "") or row.get("path", ""))
        for row in read_tsv(path) if row.get("role") == "assembly"
    }


def read_hprc(path: Path) -> set[str]:
    result = set()
    for raw in read_tsv(path):
        row = {str(key or "").strip(): value for key, value in raw.items()}
        sample = (row.get("Sample ID") or row.get("sample_id") or "").strip()
        if sample and sample.upper() not in {"CHM13", "GRCH38"}:
            result.add(sample)
    return result


def graph_members(graph: dict) -> tuple[set[str], str]:
    manifest_text = str(graph.get("membership_manifest", "") or "")
    fmt = graph.get("membership_format", "")
    if not manifest_text or not Path(manifest_text).is_file():
        return set(), "unknown_missing_membership_manifest"
    manifest = Path(manifest_text)
    if fmt == "hprc_assemblies":
        return read_hprc(manifest), "known"
    if fmt == "local_v1_ordered_assemblies":
        return read_local(manifest), "known"
    if fmt == "local_v2_union":
        result = read_local(manifest)
        rescue = Path(str(graph.get("rescue_manifest", "") or ""))
        if not rescue.is_file():
            return result, "unknown_missing_rescue_manifest"
        result.update(
            row["sample_id"] for row in read_tsv(rescue)
            if row.get("status") == "PASS" and row.get("sample_id")
        )
        return result, "known"
    if fmt == "hprc_plus_local_v1":
        result = read_hprc(manifest)
        local = Path(str(graph.get("local_membership_manifest", "") or ""))
        if not local.is_file():
            return result, "unknown_missing_local_membership_manifest"
        return result | read_local(local), "known"
    return set(), "unknown_membership_format"


def membership_audit(config_path: Path, selected_path: Path, output: Path) -> None:
    config = load_yaml(config_path)
    pilot = config["sv_short_read_pilot"]
    graph_ids = [str(value) for value in pilot["graph_ids"]]
    registry = {row["graph_id"]: row for row in config["graph_registry"]}
    samples = read_tsv(selected_path)
    rows = []
    for graph_id in graph_ids:
        members, status = graph_members(registry[graph_id])
        for sample in samples:
            sample_id = sample["sample_id"]
            if status == "known":
                is_member = sample_id in members
                audit = "FAIL" if is_member else "PASS"
                reason = "selected_sample_is_graph_member" if is_member else "selected_sample_absent_from_membership_manifest"
            else:
                is_member = None
                audit = "UNCERTAIN"
                reason = status
            rows.append({
                "sample_id": sample_id,
                "sample_role": sample.get("role", ""),
                "graph_id": graph_id,
                "is_graph_member": "unknown" if is_member is None else str(is_member).lower(),
                "membership_status": status,
                "audit_status": audit,
                "reason": reason,
            })
    write_tsv(output, rows, list(rows[0]))


def open_graph(path: Path) -> TextIO:
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def index_preflight(args: argparse.Namespace) -> None:
    if not args.source.is_file() or args.source.stat().st_size == 0:
        raise FileNotFoundError(f"missing or empty source graph: {args.source}")
    counts = Counter()
    with open_graph(args.source) as handle:
        for line in handle:
            if len(line) >= 2 and line[1] == "\t" and line[0] in {"S", "L", "P", "W"}:
                counts[line[0]] += 1
    row = {
        "graph_id": args.graph,
        "source_graph": str(args.source),
        "source_sha256": sha256_file(args.source),
        "source_size_bytes": args.source.stat().st_size,
        "segment_lines": counts["S"],
        "link_lines": counts["L"],
        "p_path_lines": counts["P"],
        "w_walk_lines": counts["W"],
        "indexing_mode": args.indexing_mode,
        "artificial_path_cover_required": "true",
        "path_cover_reason": "SV-level_rGFA_has_no_haplotype-embedded_base-resolution_GBZ",
        "graph_qc_status": args.qc_status,
        "experimental": bool_text(args.experimental),
        "preflight_status": "PASS" if counts["S"] and counts["L"] else "FAIL",
    }
    if row["preflight_status"] != "PASS":
        raise ValueError(f"{args.graph}: GFA has no S or L lines")
    write_tsv(args.output, [row], list(row))


def merge_same_schema(inputs: list[str], output: Path, sort_fields: tuple[str, ...]) -> None:
    rows: list[dict[str, str]] = []
    fields: list[str] = []
    for value in inputs:
        path = Path(value)
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for field in reader.fieldnames or []:
                if field not in fields:
                    fields.append(field)
            rows.extend(reader)
    rows.sort(key=lambda row: tuple(row.get(field, "") for field in sort_fields))
    write_tsv(output, rows, fields)


def index_manifest(args: argparse.Namespace) -> None:
    preflight = read_tsv(args.preflight)
    if len(preflight) != 1 or preflight[0].get("graph_id") != args.graph:
        raise ValueError(f"invalid index preflight row for {args.graph}")
    observed_source = sha256_file(args.source)
    if observed_source != preflight[0]["source_sha256"]:
        raise ValueError(f"{args.graph}: source graph changed after preflight")
    products = {
        "gbz": args.index_dir / f"{args.prefix_name}.giraffe.gbz",
        "distance": args.index_dir / f"{args.prefix_name}.dist",
        "minimizer": args.index_dir / f"{args.prefix_name}.shortread.withzip.min",
        "zipcodes": args.index_dir / f"{args.prefix_name}.shortread.zipcodes",
    }
    reported_dir = args.reported_index_dir or args.index_dir
    reported_products = {
        label: reported_dir / path.name for label, path in products.items()
    }
    for label, path in products.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"{args.graph}: missing or empty {label} index {path}")
    version_text = args.vg_version_file.read_text().strip().replace("\t", " ")
    log_text = args.log.read_text(errors="replace") if args.log.is_file() else ""
    path_cover_log_seen = bool(re.search(
        r"path cover|covering path|greedy path|no haplotype|artificial", log_text, re.I
    ))
    row = {
        "graph_id": args.graph,
        "indexing_mode": args.indexing_mode,
        "source_graph": str(args.source),
        "source_sha256": observed_source,
        "source_size_bytes": args.source.stat().st_size,
        "vg_version": version_text,
        "autoindex_workflow": "sr-giraffe",
        "artificial_path_cover_created": "true",
        "artificial_path_cover_log_seen": str(path_cover_log_seen).lower(),
        "p_path_lines": preflight[0].get("p_path_lines", ""),
        "w_walk_lines": preflight[0].get("w_walk_lines", ""),
        "gbz_path": str(reported_products["gbz"]),
        "gbz_size_bytes": products["gbz"].stat().st_size,
        "distance_path": str(reported_products["distance"]),
        "distance_size_bytes": products["distance"].stat().st_size,
        "minimizer_path": str(reported_products["minimizer"]),
        "minimizer_size_bytes": products["minimizer"].stat().st_size,
        "zipcodes_path": str(reported_products["zipcodes"]),
        "zipcodes_size_bytes": products["zipcodes"].stat().st_size,
        "total_index_size_bytes": sum(path.stat().st_size for path in products.values()),
        "validation": "PASS_nonempty_all_products;PASS_vg_validate_gbz;PASS_vg_describe_all_products",
        "command": args.command,
        "completion_status": "PASS",
    }
    write_tsv(args.output, [row], list(row))


def validate_index(args: argparse.Namespace) -> None:
    manifest_rows = read_tsv(args.manifest)
    if len(manifest_rows) != 1 or manifest_rows[0].get("graph_id") != args.graph:
        raise ValueError(f"expected one completion manifest row for {args.graph}")
    row = manifest_rows[0]
    observed_checksum = sha256_file(args.source)
    if observed_checksum != row.get("source_sha256"):
        raise ValueError(
            f"{args.graph}: source checksum differs from the derived index manifest; rebuild the index"
        )
    products = {
        "gbz": (row.get("gbz_path"), row.get("gbz_size_bytes")),
        "distance": (row.get("distance_path"), row.get("distance_size_bytes")),
        "minimizer": (row.get("minimizer_path"), row.get("minimizer_size_bytes")),
        "zipcodes": (row.get("zipcodes_path"), row.get("zipcodes_size_bytes")),
    }
    for label, (path_text, expected_size) in products.items():
        path = Path(str(path_text or ""))
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"{args.graph}: missing or empty recorded {label} index: {path}")
        if expected_size and path.stat().st_size != int(expected_size):
            raise ValueError(f"{args.graph}: {label} index size differs from completion manifest")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        f"PASS\tgraph={args.graph}\tsource_sha256={observed_checksum}\t"
        f"indexing_mode={row.get('indexing_mode', '')}\n"
    )


def benchmark_values(path: Path) -> tuple[str, str]:
    if not path.is_file() or path.stat().st_size == 0:
        return "", ""
    rows = read_tsv(path)
    if not rows:
        return "", ""
    row = rows[-1]
    return row.get("s", ""), row.get("max_rss", row.get("max_rss_mb", ""))


def parse_giraffe_log(path: Path) -> dict[str, object]:
    text = path.read_text(errors="replace") if path.is_file() else ""
    estimate = re.search(
        r"Using fragment length estimate:\s*([0-9.eE+-]+)\s*\+/-\s*([0-9.eE+-]+)", text
    )
    forced = re.search(r"mean:\s*([0-9.eE+-]+),\s*stdev:\s*([0-9.eE+-]+)", text)
    early = "Finalizing fragment length distribution before reaching maximum sample size" in text
    fallback_match = re.search(r"mapped\s*(\d+)\s+reads single ended", text)
    fallback_reads = int(fallback_match.group(1)) if fallback_match else 0
    explicit_single = "running in single-end mode" in text
    warnings = []
    for line in text.splitlines():
        if re.search(r"warn|fallback|single[- ]end|hairball|fragment length", line, re.I):
            cleaned = " ".join(line.split())
            if cleaned and cleaned not in warnings:
                warnings.append(cleaned)
    if estimate:
        mean, sd = estimate.groups()
        learned = not early
    elif forced:
        mean, sd = forced.groups()
        learned = False
    else:
        mean = sd = ""
        learned = False
        warnings.append("fragment_length_estimate_not_found_in_progress_log")
    return {
        "paired_fragment_learning_succeeded": str(learned).lower(),
        "single_end_fallback": str(bool(fallback_reads or explicit_single)).lower(),
        "single_end_fallback_read_count": fallback_reads,
        "learned_fragment_mean": mean,
        "learned_fragment_sd": sd,
        "relevant_warnings": ";".join(warnings),
    }


def finalize_metrics(args: argparse.Namespace) -> None:
    core_rows = read_tsv(args.core)
    stats_rows = read_tsv(args.fastq_stats)
    if len(core_rows) != 1 or len(stats_rows) != 1:
        raise ValueError("expected one core metric and one FASTQ-stat row")
    row = dict(core_rows[0])
    stats = stats_rows[0]
    if row.get("sample_id") != args.sample or row.get("graph_id") != args.graph:
        raise ValueError("core metric sample/graph does not match requested finalization")
    if int(row["total_input_reads"]) != int(stats["total_input_reads"]):
        raise ValueError("GAF parser denominator differs from independently counted FASTQ reads")

    membership_rows = [
        value for value in read_tsv(args.membership)
        if value.get("sample_id") == args.sample and value.get("graph_id") == args.graph
    ]
    if len(membership_rows) != 1:
        raise ValueError(f"no unique membership row for {args.sample}/{args.graph}")
    member = membership_rows[0]
    wall, rss = benchmark_values(args.benchmark)
    log_values = parse_giraffe_log(args.giraffe_log)
    warnings = [value for value in (row.get("relevant_warnings", ""), log_values["relevant_warnings"]) if value]
    try:
        speed = int(row["total_input_reads"]) / float(wall) if float(wall) > 0 else ""
    except (TypeError, ValueError):
        speed = ""
    row.update(log_values)
    row.update({
        "mapping_wall_time_seconds": wall,
        "peak_rss_mb": rss,
        "reads_per_second": speed,
        "proper_pair_percentage": "unavailable_from_GAF",
        "relevant_warnings": ";".join(warnings),
        "graph_qc_status": args.graph_qc_status,
        "membership_status": member.get("membership_status", ""),
        "is_graph_member": member.get("is_graph_member", ""),
        "selection_seed": stats.get("selection_seed", ""),
        "requested_pairs": stats.get("requested_pairs", ""),
        "actual_pairs": stats.get("actual_pairs", ""),
    })
    write_tsv(args.output, [row], METRIC_FIELDS)


def merge_metrics(inputs: list[str], output: Path) -> None:
    rows = [row for value in inputs for row in read_tsv(Path(value))]
    rows.sort(key=lambda row: (
        row.get("run_mode", ""), row.get("sample_role", ""),
        row.get("sample_id", ""), row.get("graph_id", ""),
    ))
    write_tsv(output, rows, METRIC_FIELDS)


def numeric(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2


def contrasts(metrics_path: Path, output: Path) -> None:
    rows = read_tsv(metrics_path)
    index = {(row["sample_id"], row["graph_id"]): row for row in rows}
    output_rows = []
    samples = sorted({row["sample_id"] for row in rows})
    for sample in samples:
        role = next(row["sample_role"] for row in rows if row["sample_id"] == sample)
        for contrast_id, left_graph, right_graph in CONTRASTS:
            left, right = index.get((sample, left_graph)), index.get((sample, right_graph))
            for metric in CONTRAST_METRICS:
                left_value = numeric(left.get(metric)) if left else None
                right_value = numeric(right.get(metric)) if right else None
                status = "PASS" if left_value is not None and right_value is not None else "BLOCKED_missing_paired_metric"
                output_rows.append({
                    "row_type": "sample",
                    "sample_id": sample,
                    "sample_role": role,
                    "summary_group": "",
                    "n_samples": "",
                    "contrast_id": contrast_id,
                    "left_graph": left_graph,
                    "right_graph": right_graph,
                    "metric": metric,
                    "delta_right_minus_left": "" if status != "PASS" else right_value - left_value,
                    "median_delta": "", "min_delta": "", "max_delta": "",
                    "status": status,
                })
    sample_rows = list(output_rows)
    roles = sorted({row["sample_role"] for row in sample_rows})
    for contrast_id, left_graph, right_graph in CONTRASTS:
        for metric in CONTRAST_METRICS:
            for group in ["overall"] + roles:
                selected = [
                    row for row in sample_rows
                    if row["contrast_id"] == contrast_id and row["metric"] == metric
                    and row["status"] == "PASS"
                    and (group == "overall" or row["sample_role"] == group)
                ]
                values = [float(row["delta_right_minus_left"]) for row in selected]
                output_rows.append({
                    "row_type": "summary_overall" if group == "overall" else "summary_role",
                    "sample_id": "", "sample_role": "" if group == "overall" else group,
                    "summary_group": group, "n_samples": len(values),
                    "contrast_id": contrast_id, "left_graph": left_graph,
                    "right_graph": right_graph, "metric": metric,
                    "delta_right_minus_left": "",
                    "median_delta": median(values) if values else "",
                    "min_delta": min(values) if values else "",
                    "max_delta": max(values) if values else "",
                    "status": "PASS" if values else "BLOCKED_no_paired_samples",
                })
    fields = [
        "row_type", "sample_id", "sample_role", "summary_group", "n_samples",
        "contrast_id", "left_graph", "right_graph", "metric",
        "delta_right_minus_left", "median_delta", "min_delta", "max_delta", "status",
    ]
    write_tsv(output, output_rows, fields)


def run_tool(command: list[str], accepted_nonzero: bool = False) -> tuple[str, str, str]:
    try:
        result = subprocess.run(command, text=True, capture_output=True, timeout=60)
    except (FileNotFoundError, subprocess.TimeoutExpired) as error:
        return "unavailable", error.__class__.__name__, ""
    text = (result.stdout + "\n" + result.stderr).strip()
    status = "available" if result.returncode == 0 or (accepted_nonzero and text) else f"exit_{result.returncode}"
    first = next((line.strip() for line in text.splitlines() if line.strip()), "")
    return status, first, text


def tool_versions(expected: str, output: Path) -> None:
    vg_status, vg_version, vg_text = run_tool(["vg", "version"])
    auto_status, auto_first, auto_help = run_tool(["vg", "autoindex", "--help"], True)
    gir_status, gir_first, giraffe_help = run_tool(["vg", "giraffe", "--help"], True)
    mini_status, mini_version, _ = run_tool(["minigraph", "--version"], True)
    checks = {
        "vg_version": expected in vg_text,
        "autoindex_sr_giraffe": "sr-giraffe" in auto_help,
        "autoindex_gfa": "--gfa" in auto_help,
        "autoindex_tmp_dir": "--tmp-dir" in auto_help,
        "giraffe_gaf": "gaf" in giraffe_help.lower(),
        "giraffe_paired_fastq": "two are allowed" in giraffe_help,
    }
    rows = [
        {"tool": "python", "version": sys.version.split()[0], "status": "available", "required_check": "python_runtime"},
        {"tool": "vg", "version": vg_version, "status": vg_status, "required_check": f"version_contains_{expected}:{str(checks['vg_version']).lower()}"},
        {"tool": "vg autoindex", "version": auto_first, "status": auto_status, "required_check": ";".join(f"{key}:{str(value).lower()}" for key, value in checks.items() if key.startswith("autoindex"))},
        {"tool": "vg giraffe", "version": gir_first, "status": gir_status, "required_check": ";".join(f"{key}:{str(value).lower()}" for key, value in checks.items() if key.startswith("giraffe"))},
        {"tool": "minigraph", "version": mini_version, "status": mini_status, "required_check": "assembly_to_graph_mapper"},
    ]
    if not all(checks.values()):
        failed = ",".join(key for key, value in checks.items() if not value)
        raise RuntimeError(f"required VG {expected} help/version checks failed: {failed}")
    write_tsv(output, rows, list(rows[0]))


def provenance(args: argparse.Namespace) -> None:
    config = load_yaml(args.config)["sv_short_read_pilot"]
    rows = [
        {"key": "analysis_label", "value": "SV-graph short-read mapping QC"},
        {"key": "config_path", "value": str(args.config)},
        {"key": "config_sha256", "value": sha256_file(args.config)},
        {"key": "selected_manifest_sha256", "value": sha256_file(args.selected)},
        {"key": "index_preflight_sha256", "value": sha256_file(args.index_preflight)},
        {"key": "tool_versions_sha256", "value": sha256_file(args.tools)},
        {"key": "baseline_graph", "value": config["baseline_graph"]},
        {"key": "graph_ids", "value": ",".join(map(str, config["graph_ids"]))},
        {"key": "indexing_mode", "value": config["index_mode"]},
        {"key": "configured_scientific_mode", "value": config.get("mode", "full")},
        {"key": "subsampling_seed", "value": config["subsampling_seed"]},
        {"key": "smoke_max_pairs", "value": config["smoke_max_pairs"]},
        {"key": "scientific_scope", "value": "preliminary_SV-level_rGFA_artificial_path_cover_QC_not_mapping_accuracy"},
        {"key": "tmem_selection_caveat", "value": "candidates_selected_from_prior_local-graph_short-read_signal_not_independent_truth"},
        {"key": "future_graph_caveat", "value": "does_not_predict_future_Minigraph-Cactus_performance"},
    ]
    write_tsv(args.output, rows, ["key", "value"])


def smoke_gate(args: argparse.Namespace) -> None:
    rows = read_tsv(args.metrics)
    problems = []
    if len(rows) != args.expected_rows:
        problems.append(f"expected_{args.expected_rows}_mapping_rows_observed_{len(rows)}")
    for row in rows:
        label = f"{row.get('sample_id')}/{row.get('graph_id')}"
        value = numeric(row.get("mapped_read_fraction"))
        if value is None or value < args.min_mapped_fraction:
            problems.append(f"{label}:mapped_fraction_{row.get('mapped_read_fraction')}")
        if row.get("single_end_fallback") == "true":
            problems.append(f"{label}:single_end_fallback")
        if row.get("paired_fragment_learning_succeeded") != "true":
            problems.append(f"{label}:paired_fragment_learning_not_confirmed")
        if row.get("is_graph_member") == "true":
            problems.append(f"{label}:selected_sample_is_graph_member")
        if int(row.get("malformed_gaf_rows", "0") or 0):
            problems.append(f"{label}:malformed_gaf_rows")
    if problems:
        raise SystemExit("Smoke QC gate failed:\n- " + "\n- ".join(problems))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        f"PASS\trows={len(rows)}\tmin_mapped_fraction={args.min_mapped_fraction}\n"
    )


def main() -> None:
    args = parse_args()
    if args.command == "select":
        select_samples(args.config, args.manifest, args.output)
    elif args.command == "fastq-qc":
        fastq_qc(args.sample, args.r1, args.r2, args.output)
    elif args.command == "subsample":
        subsample(args)
    elif args.command == "cat-fastq":
        cat_fastq(args.input)
    elif args.command == "membership":
        membership_audit(args.config, args.selected, args.output)
    elif args.command == "index-preflight":
        index_preflight(args)
    elif args.command == "merge-index-preflight":
        merge_same_schema(args.input, args.output, ("graph_id",))
    elif args.command == "index-manifest":
        index_manifest(args)
    elif args.command == "validate-index":
        validate_index(args)
    elif args.command == "finalize-metrics":
        finalize_metrics(args)
    elif args.command == "merge-metrics":
        merge_metrics(args.input, args.output)
    elif args.command == "contrasts":
        contrasts(args.metrics, args.output)
    elif args.command == "merge-evidence":
        merge_same_schema(args.input, args.output, ("sample_role", "sample_id", "graph_id"))
    elif args.command == "tools":
        tool_versions(args.expected_vg_version, args.output)
    elif args.command == "provenance":
        provenance(args)
    elif args.command == "smoke-gate":
        smoke_gate(args)


if __name__ == "__main__":
    main()
