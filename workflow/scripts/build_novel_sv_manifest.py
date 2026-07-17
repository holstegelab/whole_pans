#!/usr/bin/env python3
"""Build the joined assembly manifest used by novel-SV discovery."""

import argparse
import csv
import gzip
import os
import re
from collections import defaultdict
from pathlib import Path


ASSEMBLY_RE = re.compile(
    r"^(?P<sample>.+?)\.hifi\.hifiasm\.bp\.(?P<haplotype>hap[12])\.p_ctg(?:\..*)?$"
)

MANIFEST_FIELDS = [
    "assembly_id",
    "sample_id",
    "haplotype",
    "cleaned_fasta_path",
    "graph_member",
    "graph_order",
    "graph_role",
    "graph_context",
    "mate_status",
    "rescue_tier",
    "discovery_priority",
    "required_action",
    "original_qc_status",
    "post_qc_status",
    "contig_count",
    "total_length_bp",
    "contig_n50_bp",
    "n_percent",
    "compleasm_complete_percent",
    "compleasm_duplicated_percent",
    "best_query_aligned_percent",
    "best_reference_covered_percent",
    "best_alignment_identity_percent",
    "original_fail_reasons",
    "original_warning_reasons",
    "post_qc_fail_reasons",
    "post_qc_warning_reasons",
    "removed_contigs",
    "removed_bp",
    "review_contigs",
    "review_bp",
    "mash_distance_chm13",
    "mash_distance_graph_cohort",
    "population",
    "phenotype",
    "sequencing_batch",
    "sex",
    "paternal_haplotype",
    "read_access_status",
    "hifi_path_or_accession",
    "ont_path_or_accession",
    "read_platform",
    "read_coverage",
    "read_checksum",
]

FEASIBILITY_FIELDS = [
    "assembly_id",
    "sample_id",
    "haplotype",
    "graph_context",
    "mate_status",
    "rescue_tier",
    "discovery_priority",
    "required_action",
    "post_qc_status",
    "contig_count",
    "total_length_bp",
    "contig_n50_bp",
    "compleasm_complete_percent",
    "compleasm_duplicated_percent",
    "best_query_aligned_percent",
    "best_reference_covered_percent",
    "best_alignment_identity_percent",
    "fail_reasons",
    "warning_reasons",
    "cleaned_fasta_path",
]

PROVISIONAL_HIFI_FIELDS = [
    "sample_id",
    "assembly_ids",
    "graph_membership",
    "missing_mate_status",
    "rescue_tiers",
    "priority_tier",
    "reason_codes",
    "recommendation",
    "intended_use",
    "requested_read_type",
    "read_access_status",
    "path_or_accession",
    "platform",
    "coverage",
    "checksum",
    "population",
    "sequencing_batch",
]


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def delimiter_for(path):
    return "," if str(path).lower().endswith(".csv") else "\t"


def read_rows(path):
    if not path:
        return []
    if not Path(path).exists():
        raise FileNotFoundError(path)
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter_for(path)))


def write_rows(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def strip_fasta_suffix(value):
    name = os.path.basename(str(value))
    for suffix in (".fasta.gz", ".fna.gz", ".fa.gz", ".fasta", ".fna", ".fa"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return name


def normalize_assembly_id(value):
    name = strip_fasta_suffix(value)
    return name[: -len(".clean")] if name.endswith(".clean") else name


def sample_haplotype(assembly_id):
    match = ASSEMBLY_RE.match(normalize_assembly_id(assembly_id))
    if not match:
        raise ValueError(f"Assembly ID does not match expected hifiasm name: {assembly_id}")
    return match.group("sample"), match.group("haplotype")


def index_unique(rows, field="assembly_id"):
    result = {}
    for row in rows:
        key = normalize_assembly_id(row[field])
        if key in result:
            raise ValueError(f"Duplicate {field}: {key}")
        result[key] = row
    return result


def index_samples(rows):
    result = {}
    for row in rows:
        key = row.get("sample_id", row.get("sample", ""))
        if key:
            result[key] = row
    return result


def first(row, *names, default=""):
    for name in names:
        value = row.get(name, "") if row else ""
        if value not in {"", None, "NA", "nan"}:
            return value
    return default


def number(row, *names, default=0.0):
    value = first(row, *names, default="")
    if value == "":
        return default
    try:
        return float(str(value).replace(",", ""))
    except ValueError:
        return default


def normalize_sex(value):
    normalized = str(value).strip().lower()
    if normalized in {"m", "male", "1", "xy"}:
        return "male"
    if normalized in {"f", "female", "2", "xx"}:
        return "female"
    if normalized in {"", "0", "na", "n/a", "unknown", "not_reported"}:
        return "unknown"
    raise ValueError(f"Unrecognized sample sex value: {value!r}")


def normalize_paternal_haplotype(value):
    normalized = str(value).strip().lower()
    if normalized in {"", "na", "n/a", "unknown", "not_reported"}:
        return ""
    if normalized not in {"hap1", "hap2"}:
        raise ValueError(
            f"paternal_haplotype must be hap1 or hap2, got: {value!r}"
        )
    return normalized


def classify_rescue(post, original, previous, graph_member, restore_missing_mate):
    if graph_member:
        return "graph_member"
    status = first(post, "assembly_status", default="FAIL")
    if restore_missing_mate or status == "PASS":
        return "best_rescue"
    if status == "WARN":
        return "reasonable_rescue"
    previous_tier = first(previous, "rescue_tier")
    if previous_tier in {
        "best_rescue",
        "reasonable_rescue",
        "fragmented_rescue",
        "not_recommended",
    }:
        return previous_tier
    n50 = number(post, "contig_n50_bp", default=number(original, "contig_n50_bp"))
    length = number(post, "total_length_bp", default=number(original, "total_length_bp"))
    if n50 >= 5_000_000 and 2_700_000_000 <= length <= 3_300_000_000:
        return "fragmented_rescue"
    return "not_recommended"


def discovery_priority(graph_member, restore_missing_mate, tier):
    if restore_missing_mate:
        return 1
    if graph_member:
        return 2
    return {
        "best_rescue": 3,
        "reasonable_rescue": 4,
        "fragmented_rescue": 5,
        "not_recommended": 6,
    }.get(tier, 7)


def build(args):
    cleaned_rows = read_rows(args.cleaned_manifest)
    original = index_unique(read_rows(args.original_qc))
    post = index_unique(read_rows(args.post_qc))
    contamination = index_unique(read_rows(args.contamination))
    previous = index_unique(read_rows(args.previous_feasibility)) if args.previous_feasibility else {}
    metadata = index_samples(read_rows(args.sample_metadata)) if args.sample_metadata else {}
    reads = index_samples(read_rows(args.read_manifest)) if args.read_manifest else {}

    graph = {}
    for row in read_rows(args.graph_assemblies):
        if row.get("role", "assembly") != "assembly":
            continue
        assembly_id = normalize_assembly_id(row.get("assembly_id", row.get("path", "")))
        graph[assembly_id] = row
    graph_samples = defaultdict(set)
    for assembly_id in graph:
        sample, haplotype = sample_haplotype(assembly_id)
        graph_samples[sample].add(haplotype)

    manifest = []
    seen = set()
    for cleaned in cleaned_rows:
        assembly_id = normalize_assembly_id(cleaned["assembly_id"])
        if assembly_id in seen:
            raise ValueError(f"Duplicate cleaned assembly: {assembly_id}")
        seen.add(assembly_id)
        sample, haplotype = sample_haplotype(assembly_id)
        original_row = original.get(assembly_id, {})
        post_row = post.get(assembly_id, {})
        contamination_row = contamination.get(assembly_id, {})
        previous_row = previous.get(assembly_id, {})
        graph_row = graph.get(assembly_id, {})
        sample_metadata = metadata.get(sample, {})
        read_row = reads.get(sample, {})
        is_graph_member = assembly_id in graph
        mate_haplotype = "hap2" if haplotype == "hap1" else "hap1"
        restore = not is_graph_member and mate_haplotype in graph_samples.get(sample, set())
        tier = classify_rescue(post_row, original_row, previous_row, is_graph_member, restore)
        if is_graph_member:
            context = "current_graph_member"
            mate_status = "IN_GRAPH" if mate_haplotype in graph_samples[sample] else "MISSING_GRAPH_MATE"
            action = "retain_in_frozen_graph_control_set"
        elif restore:
            context = "restore_missing_mate_for_graph_sample"
            mate_status = "WARN"
            action = "screen_then_validate_for_incremental_test"
        else:
            context = first(previous_row, "graph_context", default="outside_graph")
            mate_status = first(previous_row, "mate_status", default="OUTSIDE_GRAPH")
            if first(post_row, "assembly_status") in {"PASS", "WARN"}:
                action = "reconsider_after_current_post_qc_then_incremental_test"
            else:
                action = first(
                    previous_row,
                    "required_action",
                    default="do_not_add_without_new_assembly_or_manual_review",
                )

        cleaned_path = first(cleaned, "path", "cleaned_path", "cleaned_fasta_path")
        row = {
            "assembly_id": assembly_id,
            "sample_id": sample,
            "haplotype": haplotype,
            "cleaned_fasta_path": cleaned_path,
            "graph_member": str(is_graph_member).lower(),
            "graph_order": first(graph_row, "order"),
            "graph_role": first(graph_row, "role", default="excluded_assembly"),
            "graph_context": context,
            "mate_status": mate_status,
            "rescue_tier": tier,
            "discovery_priority": discovery_priority(is_graph_member, restore, tier),
            "required_action": action,
            "original_qc_status": first(original_row, "assembly_status"),
            "post_qc_status": first(post_row, "assembly_status"),
            "contig_count": first(post_row, "contig_count", default=first(original_row, "contig_count")),
            "total_length_bp": first(post_row, "total_length_bp", default=first(original_row, "total_length_bp")),
            "contig_n50_bp": first(post_row, "contig_n50_bp", default=first(original_row, "contig_n50_bp")),
            "n_percent": first(post_row, "n_percent", default=first(original_row, "n_percent")),
            "compleasm_complete_percent": first(original_row, "compleasm_complete_percent"),
            "compleasm_duplicated_percent": first(original_row, "compleasm_duplicated_percent"),
            "best_query_aligned_percent": first(original_row, "best_query_aligned_percent"),
            "best_reference_covered_percent": first(original_row, "best_reference_covered_percent"),
            "best_alignment_identity_percent": first(original_row, "best_alignment_identity_percent"),
            "original_fail_reasons": first(original_row, "fail_reasons"),
            "original_warning_reasons": first(original_row, "warning_reasons"),
            "post_qc_fail_reasons": first(post_row, "fail_reasons"),
            "post_qc_warning_reasons": first(post_row, "warning_reasons"),
            "removed_contigs": first(contamination_row, "removed_contigs", default="0"),
            "removed_bp": first(contamination_row, "removed_whole_contig_bp", default="0"),
            "review_contigs": first(contamination_row, "review_contigs", default="0"),
            "review_bp": first(contamination_row, "review_bp", default="0"),
            "mash_distance_chm13": first(graph_row, "mash_distance", default=first(previous_row, "mash_distance_chm13")),
            "mash_distance_graph_cohort": first(previous_row, "mash_distance_graph_cohort"),
            "population": first(sample_metadata, "population", "ancestry"),
            "phenotype": first(sample_metadata, "phenotype", "case_control_status"),
            "sequencing_batch": first(sample_metadata, "sequencing_batch", "batch"),
            "sex": normalize_sex(first(sample_metadata, "sex", "reported_sex")),
            "paternal_haplotype": normalize_paternal_haplotype(
                first(sample_metadata, "paternal_haplotype", "paternal_hap")
            ),
            "read_access_status": first(read_row, "read_access_status", "status", default="not_checked"),
            "hifi_path_or_accession": first(read_row, "hifi_path", "hifi_accession", "path_or_accession"),
            "ont_path_or_accession": first(read_row, "ont_path", "ont_accession"),
            "read_platform": first(read_row, "platform"),
            "read_coverage": first(read_row, "coverage", "estimated_coverage"),
            "read_checksum": first(read_row, "checksum"),
        }
        manifest.append(row)

    manifest.sort(key=lambda row: (int(row["discovery_priority"]), row["sample_id"], row["haplotype"]))
    write_rows(args.output, manifest, MANIFEST_FIELDS)
    excluded = [row for row in manifest if row["graph_member"] == "false"]
    feasibility_rows = [
        {
            "assembly_id": row["assembly_id"],
            "sample_id": row["sample_id"],
            "haplotype": row["haplotype"],
            "graph_context": row["graph_context"],
            "mate_status": row["mate_status"],
            "rescue_tier": row["rescue_tier"],
            "discovery_priority": row["discovery_priority"],
            "required_action": row["required_action"],
            "post_qc_status": row["post_qc_status"],
            "contig_count": row["contig_count"],
            "total_length_bp": row["total_length_bp"],
            "contig_n50_bp": row["contig_n50_bp"],
            "compleasm_complete_percent": row["compleasm_complete_percent"],
            "compleasm_duplicated_percent": row["compleasm_duplicated_percent"],
            "best_query_aligned_percent": row["best_query_aligned_percent"],
            "best_reference_covered_percent": row["best_reference_covered_percent"],
            "best_alignment_identity_percent": row["best_alignment_identity_percent"],
            "fail_reasons": row["post_qc_fail_reasons"],
            "warning_reasons": row["post_qc_warning_reasons"],
            "cleaned_fasta_path": row["cleaned_fasta_path"],
        }
        for row in excluded
    ]
    write_rows(args.refreshed_feasibility, feasibility_rows, FEASIBILITY_FIELDS)

    by_sample = defaultdict(list)
    for row in manifest:
        by_sample[row["sample_id"]].append(row)
    provisional = []
    low_tiers = {"fragmented_rescue", "not_recommended"}
    for sample in sorted(by_sample):
        rows = sorted(by_sample[sample], key=lambda row: row["haplotype"])
        low = [row for row in rows if row["graph_member"] == "false" and row["rescue_tier"] in low_tiers]
        if not low:
            continue
        reasons = ["P1_DISCOVERY_BLIND_SPOT"]
        if len(low) == 2:
            reasons.append("BOTH_HAPLOTYPES_LOW_TIER")
        if any(row["graph_context"] == "restore_missing_mate_for_graph_sample" for row in rows):
            reasons.append("RESTORE_MISSING_GRAPH_MATE")
        first_row = rows[0]
        provisional.append(
            {
                "sample_id": sample,
                "assembly_ids": ";".join(row["assembly_id"] for row in rows),
                "graph_membership": ";".join(f"{row['haplotype']}={row['graph_member']}" for row in rows),
                "missing_mate_status": ";".join(sorted({row["mate_status"] for row in rows})),
                "rescue_tiers": ";".join(f"{row['haplotype']}={row['rescue_tier']}" for row in rows),
                "priority_tier": "P1_DISCOVERY_BLIND_SPOT",
                "reason_codes": ";".join(reasons),
                "recommendation": "Retrieve sample-level HiFi reads for read-based SV discovery; do not interpret assembly no-calls as reference genotypes.",
                "intended_use": "discovery",
                "requested_read_type": "HiFi",
                "read_access_status": first_row["read_access_status"],
                "path_or_accession": first_row["hifi_path_or_accession"],
                "platform": first_row["read_platform"],
                "coverage": first_row["read_coverage"],
                "checksum": first_row["read_checksum"],
                "population": first_row["population"],
                "sequencing_batch": first_row["sequencing_batch"],
            }
        )
    provisional.sort(
        key=lambda row: (
            "BOTH_HAPLOTYPES_LOW_TIER" not in row["reason_codes"],
            row["sample_id"],
        )
    )
    write_rows(args.provisional_hifi, provisional, PROVISIONAL_HIFI_FIELDS)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cleaned-manifest", required=True)
    parser.add_argument("--graph-assemblies", required=True)
    parser.add_argument("--original-qc", required=True)
    parser.add_argument("--post-qc", required=True)
    parser.add_argument("--contamination", required=True)
    parser.add_argument("--previous-feasibility")
    parser.add_argument("--sample-metadata")
    parser.add_argument("--read-manifest")
    parser.add_argument("--output", required=True)
    parser.add_argument("--refreshed-feasibility", required=True)
    parser.add_argument("--provisional-hifi", required=True)
    return parser.parse_args()


if __name__ == "__main__":
    build(parse_args())
