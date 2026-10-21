#!/usr/bin/env python3
"""Build graph/sample registries and gate expensive benchmark stages."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import re
import subprocess
from pathlib import Path

import yaml


ROLES = {
    "local_v1_member_control",
    "hprc_member_control",
    "v2_rescue_source_control",
    "held_out_primary",
    "tmem_candidate",
    "tmem_matched_control",
    "tmem_random_blinded",
}
GRAPH_FIELDS = [
    "graph_id", "label", "comparison_scope", "graph_type", "backbone",
    "reference_version", "path", "index_dist", "index_minimizer",
    "index_zipcodes", "membership_manifest", "checksum", "checksum_file",
    "mapping_capabilities", "experimental", "qc_status", "provenance",
]
PREFLIGHT_FIELDS = [
    "graph_id", "available", "severity", "graph_type", "backbone_compatible",
    "representation_check", "long_read_mapper", "short_read_mapper",
    "short_read_compatible", "membership_auditable", "checksum_status",
    "size_bytes", "problems", "qc_status", "experimental",
]
SAMPLE_FIELDS = [
    "sample_id", "data_type", "long_read", "r1", "r2",
    "nominal_coverage", "selected_coverage", "cohort", "population", "sex",
    "batch", "archive_cohort", "stage_from_archive", "truth_assembly_hap1",
    "truth_assembly_hap2", "truth_callset", "callable_regions", "role",
    "include", "inclusion_reason", "exclusion_reason",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    pre = sub.add_parser("preflight")
    pre.add_argument("--config", required=True, type=Path)
    pre.add_argument("--registry", required=True, type=Path)
    pre.add_argument("--preflight", required=True, type=Path)
    pre.add_argument("--report", required=True, type=Path)

    samples = sub.add_parser("samples")
    samples.add_argument("--config", required=True, type=Path)
    samples.add_argument("--tmem-manifest", required=True, type=Path)
    samples.add_argument("--samples", required=True, type=Path)
    samples.add_argument("--membership", required=True, type=Path)
    samples.add_argument("--leakage", required=True, type=Path)

    gate = sub.add_parser("gate")
    gate.add_argument("--preflight", required=True, type=Path)
    gate.add_argument("--leakage", required=True, type=Path)
    gate.add_argument("--mode", choices=("long", "short", "rescue"), required=True)
    gate.add_argument("--graph-id", action="append", default=[])
    gate.add_argument("--output", required=True, type=Path)

    versions = sub.add_parser("versions")
    versions.add_argument("--output", required=True, type=Path)

    provenance = sub.add_parser("provenance")
    provenance.add_argument("--config", required=True, type=Path)
    provenance.add_argument("--preflight", required=True, type=Path)
    provenance.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def load_config(path: Path) -> dict:
    with path.open() as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"{path}: expected a YAML mapping")
    return value


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else path.open()


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
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


def bool_text(value: object) -> str:
    return "true" if bool(value) else "false"


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def expected_checksum(graph: dict) -> tuple[str, str]:
    configured = str(graph.get("checksum", "") or "").strip().split()[0:1]
    if configured:
        return configured[0], "configured"
    checksum_file = str(graph.get("checksum_file", "") or "").strip()
    if checksum_file:
        path = Path(checksum_file)
        if path.is_file() and path.stat().st_size:
            value = path.read_text().strip().split()[0]
            return value, "checksum_file"
        return "", "checksum_file_missing"
    return "", "not_available"


def inspect_gfa(path: Path) -> str:
    saw_segment = False
    saw_rgfa = False
    try:
        with open_text(path) as handle:
            for line_number, line in enumerate(handle, start=1):
                if line.startswith("S\t"):
                    saw_segment = True
                    fields = line.rstrip("\n").split("\t")
                    tags = fields[3:]
                    if any(tag.startswith("SN:Z:") for tag in tags) and any(
                        tag.startswith("SO:i:") for tag in tags
                    ):
                        saw_rgfa = True
                        break
                if line_number >= 5000:
                    break
    except (OSError, UnicodeError) as error:
        return f"unreadable:{error.__class__.__name__}"
    if saw_rgfa:
        return "rgfa_SN_SO_tags_seen"
    if saw_segment:
        return "gfa_segments_seen_rgfa_tags_not_seen"
    return "no_GFA_segments_seen"


def preflight(config: dict) -> tuple[list[dict], list[dict]]:
    graphs = config.get("graph_registry", [])
    if not isinstance(graphs, list):
        raise ValueError("graph_registry must be a YAML list so duplicate IDs can be audited")
    graph_ids = [str(graph.get("graph_id", "")) for graph in graphs]
    duplicates = {graph_id for graph_id in graph_ids if graph_ids.count(graph_id) > 1}
    registry_rows: list[dict] = []
    status_rows: list[dict] = []

    for graph in graphs:
        graph_id = str(graph.get("graph_id", ""))
        path_text = str(graph.get("path", "") or "")
        path = Path(path_text) if path_text else Path("/__missing_graph_path__")
        graph_type = str(graph.get("graph_type", ""))
        capabilities = set(graph.get("mapping_capabilities", []))
        indexes = graph.get("indexes", {}) or {}
        problems: list[str] = []
        warnings: list[str] = []
        size = 0

        if not graph_id:
            problems.append("missing_graph_id")
        if graph_id in duplicates:
            problems.append("duplicate_graph_id")
        if not path_text:
            problems.append("missing_graph_path")
        elif any(
            path.name.endswith(suffix) or f"{suffix}." in path.name
            for suffix in (".tmp", ".part", ".partial", ".incomplete")
        ):
            problems.append("temporary_path_declared_as_final")
        elif not path.is_file():
            problems.append("graph_missing")
            temporary_companions = [
                Path(path_text + suffix)
                for suffix in (".tmp", ".part", ".partial", ".incomplete")
            ]
            if any(candidate.exists() for candidate in temporary_companions):
                problems.append("temporary_incomplete_artifact_present")
        else:
            size = path.stat().st_size
            if size == 0:
                problems.append("graph_empty")

        backbone_ok = True
        if graph.get("comparison_scope") == "headline":
            backbone_ok = str(graph.get("backbone", "")).upper() == "CHM13"
            if not backbone_ok:
                problems.append("headline_backbone_not_CHM13")
            if "CHM13" not in str(graph.get("reference_version", "")).upper():
                problems.append("headline_reference_version_not_explicit_CHM13")

        representation = "not_inspected"
        long_mapper = "unsupported"
        short_mapper = "unsupported"
        short_compatible = False
        if graph_type == "linear_fasta":
            long_mapper = "minigraph;minimap2"
            representation = "single_path_fasta"
            if "short_reads" in capabilities:
                problems.append("linear_short_mapper_not_implemented")
        elif graph_type == "minigraph_rgfa":
            long_mapper = "minigraph"
            if path.is_file() and size:
                representation = inspect_gfa(path)
                if representation.startswith("unreadable") or representation == "no_GFA_segments_seen":
                    problems.append("declared_minigraph_graph_not_readable_as_GFA")
                elif representation == "gfa_segments_seen_rgfa_tags_not_seen":
                    warnings.append("rgfa_tags_not_seen_in_first_5000_lines")
            if "short_reads" in capabilities:
                problems.append("minigraph_rgfa_not_directly_supported_by_giraffe")
        elif graph_type == "gbz":
            representation = "gbz_with_explicit_indexes"
            required_indexes = ("dist", "minimizer", "zipcodes")
            missing_indexes = []
            for key in required_indexes:
                index_text = str(indexes.get(key, "") or "")
                index_path = Path(index_text) if index_text else None
                if not index_path or not index_path.is_file() or index_path.stat().st_size == 0:
                    missing_indexes.append(key)
            if missing_indexes:
                problems.append("missing_or_empty_indexes:" + ",".join(missing_indexes))
            else:
                short_mapper = "vg_giraffe"
                short_compatible = "short_reads" in capabilities
        else:
            problems.append("unsupported_graph_type")

        if "long_reads" in capabilities and long_mapper == "unsupported":
            problems.append("declared_long_read_capability_incompatible")
        if "short_reads" in capabilities and not short_compatible:
            problems.append("declared_short_read_capability_incompatible")

        membership = str(graph.get("membership_manifest", "") or "")
        membership_auditable = bool(
            membership and Path(membership).is_file() and Path(membership).stat().st_size
        )
        if graph.get("comparison_scope") == "headline" and graph_type != "linear_fasta":
            if not membership_auditable:
                problems.append("headline_membership_manifest_missing")
        if graph.get("membership_format") == "hprc_plus_local_v1":
            local_manifest = str(graph.get("local_membership_manifest", "") or "")
            if not local_manifest or not Path(local_manifest).is_file():
                problems.append("local_union_membership_manifest_missing")

        provenance = str(graph.get("provenance", "") or "")
        if provenance and "://" not in provenance and "/" in provenance:
            provenance_path = Path(provenance)
            if not provenance_path.is_file() or provenance_path.stat().st_size == 0:
                problems.append("provenance_artifact_missing_or_empty")
        qc_report = str(graph.get("qc_report", "") or "")
        if qc_report:
            qc_path = Path(qc_report)
            if not qc_path.is_file() or qc_path.stat().st_size == 0:
                problems.append("qc_report_missing_or_empty")

        expected, checksum_source = expected_checksum(graph)
        checksum_status = checksum_source
        if expected and path.is_file() and size:
            if graph.get("verify_checksum", False):
                observed = file_sha256(path)
                checksum_status = "verified" if observed == expected else "mismatch"
                if observed != expected:
                    problems.append("checksum_mismatch")
            else:
                checksum_status = "recorded_not_recomputed"

        registry_rows.append({
            "graph_id": graph_id,
            "label": graph.get("label", ""),
            "comparison_scope": graph.get("comparison_scope", ""),
            "graph_type": graph_type,
            "backbone": graph.get("backbone", ""),
            "reference_version": graph.get("reference_version", ""),
            "path": path_text,
            "index_dist": indexes.get("dist", ""),
            "index_minimizer": indexes.get("minimizer", ""),
            "index_zipcodes": indexes.get("zipcodes", ""),
            "membership_manifest": membership,
            "checksum": expected,
            "checksum_file": graph.get("checksum_file", ""),
            "mapping_capabilities": ";".join(sorted(capabilities)),
            "experimental": bool_text(graph.get("experimental", False)),
            "qc_status": graph.get("qc_status", ""),
            "provenance": graph.get("provenance", ""),
        })
        available = not problems
        all_notes = problems + ["WARN:" + value for value in warnings]
        status_rows.append({
            "graph_id": graph_id,
            "available": bool_text(available),
            "severity": "FAIL" if problems else ("WARN" if warnings else "PASS"),
            "graph_type": graph_type,
            "backbone_compatible": bool_text(backbone_ok),
            "representation_check": representation,
            "long_read_mapper": long_mapper,
            "short_read_mapper": short_mapper,
            "short_read_compatible": bool_text(short_compatible),
            "membership_auditable": bool_text(membership_auditable),
            "checksum_status": checksum_status,
            "size_bytes": size,
            "problems": ";".join(all_notes),
            "qc_status": graph.get("qc_status", ""),
            "experimental": bool_text(graph.get("experimental", False)),
        })
    return registry_rows, status_rows


def write_preflight_report(path: Path, registry: list[dict], status: list[dict]) -> None:
    labels = {row["graph_id"]: row["label"] for row in registry}
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Read graph benchmark preflight", "",
        "Preflight inspects artifacts and mapper compatibility only; it does not validate biological completeness.", "",
        "| Graph | Available | Type | Long mapper | Short mapper | QC/provenance status | Problems |",
        "|---|---:|---|---|---|---|---|",
    ]
    for row in status:
        lines.append(
            f"| {labels.get(row['graph_id'], row['graph_id'])} (`{row['graph_id']}`) "
            f"| {row['available']} | {row['graph_type']} | {row['long_read_mapper']} "
            f"| {row['short_read_mapper']} | {row['qc_status']} | {row['problems'] or 'none'} |"
        )
    lines.extend([
        "", "## Interpretation", "",
        "- Only CHM13-backbone whole-genome targets may enter the headline long-read contrasts.",
        "- Minigraph SV-level rGFA targets use Minigraph for like-for-like long-read mapping.",
        "- `vg giraffe` is enabled only for GBZ targets with explicit dist/minimizer/zipcodes indexes.",
        "- Local V2 remains experimental because cohort-accounting QC is FAIL even when structural integrity passes.",
    ])
    temporary = Path(str(path) + ".tmp")
    temporary.write_text("\n".join(lines) + "\n")
    temporary.replace(path)


def donor_from_assembly(value: str) -> str:
    name = Path(value).name
    return re.split(r"\.hifi\.|_(?:hap[12]|mat|pat)(?:_|\.|$)", name, maxsplit=1)[0]


def read_local_v1(path: Path) -> set[str]:
    donors: set[str] = set()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("role") == "assembly":
                donors.add(donor_from_assembly(row.get("assembly_id", "") or row.get("path", "")))
    return donors


def read_rescue_sources(path: Path) -> set[str]:
    with path.open(newline="") as handle:
        return {
            row["sample_id"] for row in csv.DictReader(handle, delimiter="\t")
            if row.get("status") == "PASS" and row.get("sample_id")
        }


def read_hprc(path: Path) -> set[str]:
    donors: set[str] = set()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            row = {str(key or "").strip(): value for key, value in raw.items()}
            sample = (row.get("Sample ID") or row.get("sample_id") or "").strip()
            if sample and sample.upper() not in {"CHM13", "GRCH38"}:
                donors.add(sample)
    return donors


def graph_members(graph: dict) -> tuple[set[str], str]:
    fmt = graph.get("membership_format", "")
    manifest = str(graph.get("membership_manifest", "") or "")
    if graph.get("graph_type") == "linear_fasta":
        return set(), "known_reference_only"
    if not manifest or not Path(manifest).is_file():
        return set(), "unknown_missing_manifest"
    if fmt == "local_v1_ordered_assemblies":
        return read_local_v1(Path(manifest)), "known"
    if fmt == "local_v2_union":
        members = read_local_v1(Path(manifest))
        rescue = str(graph.get("rescue_manifest", "") or "")
        if not rescue or not Path(rescue).is_file():
            return members, "unknown_missing_rescue_manifest"
        return members | read_rescue_sources(Path(rescue)), "known"
    if fmt == "hprc_assemblies":
        return read_hprc(Path(manifest)), "known"
    if fmt == "hprc_plus_local_v1":
        local = str(graph.get("local_membership_manifest", "") or "")
        if not local or not Path(local).is_file():
            return read_hprc(Path(manifest)), "unknown_missing_local_manifest"
        return read_hprc(Path(manifest)) | read_local_v1(Path(local)), "known"
    return set(), "unknown_membership_format"


def load_tmem_samples(path: Path) -> list[dict]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [{
        "sample_id": row.get("sample_id", ""),
        "data_type": "illumina_pe",
        "long_read": "",
        "r1": row.get("r1", ""),
        "r2": row.get("r2", ""),
        "nominal_coverage": row.get("coverage_proxy", ""),
        "selected_coverage": row.get("coverage_proxy", ""),
        "cohort": "DemenTree",
        "population": "",
        "sex": "",
        "batch": row.get("batch", ""),
        "truth_assembly_hap1": "",
        "truth_assembly_hap2": "",
        "truth_callset": "",
        "callable_regions": "",
        "role": row.get("role", ""),
        "include": "true",
        "inclusion_reason": row.get("selection_reason", ""),
        "exclusion_reason": "",
    } for row in rows]


def prepare_samples(config: dict, tmem_path: Path) -> tuple[list[dict], list[dict], list[dict]]:
    rows = []
    seen: set[str] = set()
    for raw in list(config.get("samples", [])) + load_tmem_samples(tmem_path):
        row = {field: raw.get(field, "") for field in SAMPLE_FIELDS}
        sample = str(row["sample_id"])
        if not sample:
            raise ValueError("sample entry has an empty sample_id")
        if sample in seen:
            raise ValueError(f"duplicate sample_id: {sample}")
        seen.add(sample)
        if row["role"] not in ROLES:
            raise ValueError(f"{sample}: unsupported role {row['role']!r}")
        row["include"] = bool_text(str(row["include"]).lower() not in {"false", "0", "no", ""})
        row["stage_from_archive"] = bool_text(
            str(row["stage_from_archive"]).lower() not in {"false", "0", "no", ""}
        )
        rows.append(row)

    graphs = config.get("graph_registry", [])
    memberships = {}
    membership_status = {}
    for graph in graphs:
        members, status = graph_members(graph)
        graph_id = graph["graph_id"]
        memberships[graph_id] = members
        membership_status[graph_id] = status

    matrix = []
    audit = []
    for row in rows:
        sample = str(row["sample_id"])
        for graph in graphs:
            graph_id = graph["graph_id"]
            status = membership_status[graph_id]
            member = sample in memberships[graph_id]
            matrix.append({
                "sample_id": sample,
                "sample_role": row["role"],
                "graph_id": graph_id,
                "is_member": bool_text(member) if status.startswith("known") else "unknown",
                "membership_status": status,
            })
            if row["role"] == "held_out_primary" and graph.get("comparison_scope") == "headline":
                if status.startswith("unknown"):
                    result = "FAIL"
                    reason = "membership_not_auditable"
                elif member:
                    result = "FAIL"
                    reason = "held_out_sample_is_graph_member"
                else:
                    result = "PASS"
                    reason = "sample_absent_from_graph_membership"
                audit.append({
                    "sample_id": sample,
                    "sample_role": row["role"],
                    "graph_id": graph_id,
                    "audit_status": result,
                    "reason": reason,
                })
    if not audit:
        audit.append({
            "sample_id": "",
            "sample_role": "held_out_primary",
            "graph_id": "",
            "audit_status": "NOT_RUN",
            "reason": "no_held_out_primary_samples_configured_headline_results_unavailable",
        })
    return rows, matrix, audit


def read_rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run_gate(args: argparse.Namespace) -> None:
    preflight_rows = read_rows(args.preflight)
    leakage_rows = read_rows(args.leakage)
    requested = set(args.graph_id)
    problems = []
    for row in preflight_rows:
        graph_id = row["graph_id"]
        relevant = graph_id in requested if requested else (
            (args.mode in {"long", "rescue"} and row["long_read_mapper"] != "unsupported")
            or (args.mode == "short" and row["short_read_mapper"] != "unsupported")
        )
        if relevant and row["available"] != "true":
            problems.append(f"{graph_id}: {row['problems'] or 'unavailable'}")
        if args.mode == "short" and relevant and row["short_read_compatible"] != "true":
            problems.append(f"{graph_id}: not short-read compatible")
    for row in leakage_rows:
        if row.get("audit_status") == "FAIL":
            problems.append(
                f"held-out leakage {row.get('sample_id')} / {row.get('graph_id')}: {row.get('reason')}"
            )
    if problems:
        raise SystemExit(
            "Preflight gate failed before expensive jobs:\n- " + "\n- ".join(problems)
            + f"\nInspect {args.preflight} and {args.leakage}."
        )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(f"PASS\tmode={args.mode}\n")


def tool_versions(output: Path) -> None:
    commands = {
        "python": ["python", "--version"],
        "minigraph": ["minigraph", "--version"],
        "minimap2": ["minimap2", "--version"],
        "seqtk": ["seqtk"],
        "vg": ["vg", "version"],
    }
    rows = []
    for tool, command in commands.items():
        try:
            result = subprocess.run(command, text=True, capture_output=True, timeout=30)
            text = (result.stdout or result.stderr).strip().splitlines()
            version = text[0] if text else f"exit_{result.returncode}"
            status = "available" if result.returncode == 0 or tool == "seqtk" else "error"
        except (FileNotFoundError, subprocess.TimeoutExpired) as error:
            version = error.__class__.__name__
            status = "unavailable"
        rows.append({"tool": tool, "version": version, "status": status})
    write_tsv(output, rows, ["tool", "version", "status"])


def write_provenance(config_path: Path, preflight_path: Path, output: Path) -> None:
    config = load_config(config_path)
    rows = [
        {"key": "config_path", "value": str(config_path)},
        {"key": "config_sha256", "value": file_sha256(config_path)},
        {"key": "selection_seed", "value": config.get("long_reads", {}).get("selection_seed", "")},
        {"key": "tmem_random_seed", "value": config.get("tmem_pilot", {}).get("random_seed", "")},
        {"key": "simulation_seed", "value": config.get("simulation", {}).get("seed", "")},
        {"key": "v2_qc_interpretation", "value": "EXPERIMENTAL: structural integrity passes but overall cohort-accounting QC FAIL"},
        {"key": "preflight_sha256", "value": file_sha256(preflight_path)},
    ]
    for row in read_rows(preflight_path):
        rows.append({"key": f"graph.{row['graph_id']}.available", "value": row["available"]})
        rows.append({"key": f"graph.{row['graph_id']}.qc_status", "value": row["qc_status"]})
    write_tsv(output, rows, ["key", "value"])


def main() -> None:
    args = parse_args()
    if args.command == "preflight":
        config = load_config(args.config)
        registry, status = preflight(config)
        write_tsv(args.registry, registry, GRAPH_FIELDS)
        write_tsv(args.preflight, status, PREFLIGHT_FIELDS)
        write_preflight_report(args.report, registry, status)
    elif args.command == "samples":
        rows, matrix, audit = prepare_samples(load_config(args.config), args.tmem_manifest)
        write_tsv(args.samples, rows, SAMPLE_FIELDS)
        write_tsv(
            args.membership, matrix,
            ["sample_id", "sample_role", "graph_id", "is_member", "membership_status"],
        )
        write_tsv(
            args.leakage, audit,
            ["sample_id", "sample_role", "graph_id", "audit_status", "reason"],
        )
    elif args.command == "gate":
        run_gate(args)
    elif args.command == "versions":
        tool_versions(args.output)
    elif args.command == "provenance":
        write_provenance(args.config, args.preflight, args.output)


if __name__ == "__main__":
    main()
