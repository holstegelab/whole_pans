#!/usr/bin/env python3
"""Rank sample-level HiFi retrieval requests after assembly-only discovery."""

import argparse
import bisect
import csv
import gzip
import sys
from collections import Counter, defaultdict
from pathlib import Path


FIELDS = [
    "rank",
    "sample_id",
    "assembly_ids",
    "graph_membership",
    "missing_mate_status",
    "rescue_tiers",
    "callable_fractions",
    "uncallable_haplotype_count",
    "priority_tier",
    "reason_codes",
    "recommendation",
    "candidate_count",
    "candidate_counts_by_method",
    "candidate_counts_by_confidence",
    "candidate_counts_by_svtype",
    "important_candidate_ids",
    "population",
    "sequencing_batch",
    "intended_use",
    "requested_read_type",
    "read_access_status",
    "path_or_accession",
    "platform",
    "coverage",
    "checksum",
]


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def allow_large_tsv_fields():
    """Raise csv's conservative default field limit for aggregated catalogs."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def read_tsv(path):
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def iter_tsv(path):
    with open_text(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def split_values(value):
    return [item for item in str(value).split(";") if item]


def fraction(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return 0.0


def counter_text(counter):
    return ";".join(f"{key}={counter[key]}" for key in sorted(counter))


def summarize_catalog(path, important_limit=25):
    """Collect the small per-sample summary needed from a potentially huge catalog."""
    summaries = defaultdict(
        lambda: {
            "candidate_count": 0,
            "pending_hifi_count": 0,
            "high_confidence_count": 0,
            "difficult_candidate_count": 0,
            "method": Counter(),
            "confidence": Counter(),
            "svtype": Counter(),
            "important_candidates": [],
        }
    )
    difficult_types = {
        "BND",
        "TRA",
        "CTX",
        "INV",
        "DUP",
        "DUP:TANDEM",
        "DUP:INT",
        "CNV",
        "COMPLEX_INDEL",
    }
    for row in iter_tsv(path):
        pending = row.get("validation_status") == "PENDING_HIFI"
        high = row.get("confidence") == "HIGH"
        difficult = row.get("svtype") in difficult_types
        candidate = (
            not pending,
            not high,
            not difficult,
            row["event_id"],
        )
        for sample in split_values(row.get("carrier_samples", "")):
            summary = summaries[sample]
            summary["candidate_count"] += 1
            summary["pending_hifi_count"] += pending
            summary["high_confidence_count"] += high
            summary["difficult_candidate_count"] += difficult
            for method in split_values(row.get("discovery_methods", "")):
                summary["method"][method] += 1
            summary["confidence"][row.get("confidence", "unknown")] += 1
            summary["svtype"][row.get("svtype", "unknown")] += 1
            bisect.insort(summary["important_candidates"], candidate)
            if len(summary["important_candidates"]) > important_limit:
                summary["important_candidates"].pop()
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-manifest", required=True)
    parser.add_argument("--screen-summary", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument(
        "--evidence",
        help="Deprecated; candidate summaries are now derived from the filtered catalog.",
    )
    parser.add_argument("--output", required=True)
    parser.add_argument("--callable-threshold", type=float, default=0.85)
    parser.add_argument("--validation-count", type=int, default=50)
    parser.add_argument("--control-count", type=int, default=20)
    args = parser.parse_args()

    allow_large_tsv_fields()
    assemblies = read_tsv(args.assembly_manifest)
    screen = {row["assembly_id"]: row for row in read_tsv(args.screen_summary)}

    by_sample = defaultdict(list)
    for row in assemblies:
        by_sample[row["sample_id"]].append(row)
    catalog_by_sample = summarize_catalog(args.catalog)
    low_tiers = {"fragmented_rescue", "not_recommended"}
    sample_rows_for_ranking = []
    for sample, sample_rows in sorted(by_sample.items()):
        sample_rows.sort(key=lambda row: row["haplotype"])
        sample_catalog = catalog_by_sample.get(
            sample,
            {
                "candidate_count": 0,
                "pending_hifi_count": 0,
                "high_confidence_count": 0,
                "difficult_candidate_count": 0,
                "method": Counter(),
                "confidence": Counter(),
                "svtype": Counter(),
                "important_candidates": [],
            },
        )
        callability = {}
        for row in sample_rows:
            summary = screen.get(row["assembly_id"], {})
            callability[row["haplotype"]] = max(
                fraction(summary.get("primary_callable_fraction")),
                fraction(summary.get("sensitivity_callable_fraction")),
            )
        uncallable = sum(value < args.callable_threshold for value in callability.values())
        low_tier = any(
            row.get("graph_member") == "false"
            and row.get("rescue_tier") in low_tiers
            for row in sample_rows
        )
        missing_mate = any(
            row.get("graph_context") == "restore_missing_mate_for_graph_sample"
            for row in sample_rows
        )
        priority_candidate_count = (
            sample_catalog["pending_hifi_count"]
            + sample_catalog["high_confidence_count"]
            + sample_catalog["difficult_candidate_count"]
        )
        first = sample_rows[0]
        output_row = {
            "sample_id": sample,
            "assembly_ids": ";".join(row["assembly_id"] for row in sample_rows),
            "graph_membership": ";".join(f"{row['haplotype']}={row['graph_member']}" for row in sample_rows),
            "missing_mate_status": ";".join(sorted({row["mate_status"] for row in sample_rows})),
            "rescue_tiers": ";".join(f"{row['haplotype']}={row['rescue_tier']}" for row in sample_rows),
            "callable_fractions": ";".join(f"{hap}={callability[hap]:.6f}" for hap in sorted(callability)),
            "uncallable_haplotype_count": uncallable,
            "candidate_count": sample_catalog["candidate_count"],
            "candidate_counts_by_method": counter_text(sample_catalog["method"]),
            "candidate_counts_by_confidence": counter_text(sample_catalog["confidence"]),
            "candidate_counts_by_svtype": counter_text(sample_catalog["svtype"]),
            "important_candidate_ids": ";".join(
                candidate[3] for candidate in sample_catalog["important_candidates"]
            ),
            "population": first.get("population", ""),
            "sequencing_batch": first.get("sequencing_batch", ""),
            "intended_use": "",
            "requested_read_type": "HiFi",
            "read_access_status": first.get("read_access_status", "not_checked"),
            "path_or_accession": first.get("hifi_path_or_accession", ""),
            "platform": first.get("read_platform", ""),
            "coverage": first.get("read_coverage", ""),
            "checksum": first.get("read_checksum", ""),
            "__min_callable": min(callability.values()) if callability else 0.0,
            "__blind_spot": bool(low_tier or uncallable or missing_mate),
            "__missing_mate": missing_mate,
            "__graph_member": any(
                row.get("graph_member") == "true" for row in sample_rows
            ),
            "__priority_candidate_count": priority_candidate_count,
            "__pending_count": sample_catalog["pending_hifi_count"],
            "__high_count": sample_catalog["high_confidence_count"],
            "__difficult_count": sample_catalog["difficult_candidate_count"],
        }
        sample_rows_for_ranking.append(output_row)

    mandatory = [row for row in sample_rows_for_ranking if row["__blind_spot"]]
    validation_pool = [
        row
        for row in sample_rows_for_ranking
        if not row["__blind_spot"] and row["__priority_candidate_count"] > 0
    ]
    validation_pool.sort(
        key=lambda row: (
            -row["__pending_count"],
            -row["__high_count"],
            -row["__difficult_count"],
            -int(row["candidate_count"]),
            row["sample_id"],
        )
    )
    selected_validation = validation_pool[: args.validation_count]
    selected_ids = {row["sample_id"] for row in mandatory + selected_validation}
    control_pool = [
        row
        for row in sample_rows_for_ranking
        if row["sample_id"] not in selected_ids
        and row["__graph_member"]
        and int(row["uncallable_haplotype_count"]) == 0
    ]
    control_pool.sort(
        key=lambda row: (
            int(row["candidate_count"]) != 0,
            int(row["candidate_count"]),
            -row["__min_callable"],
            row["sample_id"],
        )
    )
    controls = control_pool[: args.control_count]

    for row in mandatory:
        reasons = ["P1_DISCOVERY_BLIND_SPOT"]
        intended = ["discovery"]
        if row["__priority_candidate_count"]:
            reasons.append("P1_VALIDATE_CANDIDATE")
            intended.append("validation")
        if row["__missing_mate"]:
            reasons.append("RESTORE_MISSING_GRAPH_MATE")
        row.update(
            {
                "priority_tier": "P1_DISCOVERY_BLIND_SPOT",
                "reason_codes": ";".join(reasons),
                "recommendation": "Retrieve HiFi reads for sensitive read-based discovery; assembly no-calls are not negative evidence.",
                "intended_use": ";".join(intended),
            }
        )
    for row in selected_validation:
        row.update(
            {
                "priority_tier": "P1_VALIDATE_CANDIDATE",
                "reason_codes": "P1_VALIDATE_CANDIDATE",
                "recommendation": "Retrieve HiFi reads for breakpoint-spanning validation of prioritized graph-residual candidates.",
                "intended_use": "validation",
            }
        )
    for row in controls:
        row.update(
            {
                "priority_tier": "P3_CONTROL",
                "reason_codes": "P3_CONTROL",
                "recommendation": "Retrieve or reuse HiFi reads as a graph-member assembly/read concordance control.",
                "intended_use": "control",
            }
        )

    ranked = mandatory + selected_validation + controls
    tier_order = {
        "P1_DISCOVERY_BLIND_SPOT": 0,
        "P1_VALIDATE_CANDIDATE": 1,
        "P3_CONTROL": 3,
    }
    ranked.sort(
        key=lambda row: (
            tier_order[row["priority_tier"]],
            -int(row["uncallable_haplotype_count"]),
            row["__min_callable"],
            -int(row["candidate_count"]),
            row["sample_id"],
        )
    )
    for rank, row in enumerate(ranked, start=1):
        row["rank"] = rank
        for field in list(row):
            if field.startswith("__"):
                row.pop(field)
    write_tsv(args.output, ranked)


if __name__ == "__main__":
    main()
