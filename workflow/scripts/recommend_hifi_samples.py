#!/usr/bin/env python3
"""Rank sample-level HiFi retrieval requests after assembly-only discovery."""

import argparse
import csv
import gzip
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assembly-manifest", required=True)
    parser.add_argument("--screen-summary", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--evidence", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--callable-threshold", type=float, default=0.85)
    parser.add_argument("--control-count", type=int, default=20)
    args = parser.parse_args()

    assemblies = read_tsv(args.assembly_manifest)
    screen = {row["assembly_id"]: row for row in read_tsv(args.screen_summary)}
    catalog = read_tsv(args.catalog)

    by_sample = defaultdict(list)
    for row in assemblies:
        by_sample[row["sample_id"]].append(row)
    evidence_counts = defaultdict(
        lambda: {
            "method": Counter(),
            "confidence": Counter(),
            "svtype": Counter(),
        }
    )
    for row in iter_tsv(args.evidence):
        if row.get("sample_id"):
            counts = evidence_counts[row["sample_id"]]
            counts["method"][row.get("discovery_method", "unknown")] += 1
            counts["confidence"][row.get("confidence_tier", "unknown")] += 1
            counts["svtype"][row.get("svtype", "unknown")] += 1
    catalog_by_sample = defaultdict(list)
    for row in catalog:
        for sample in split_values(row.get("carrier_samples", "")):
            catalog_by_sample[sample].append(row)

    low_tiers = {"fragmented_rescue", "not_recommended"}
    difficult_types = {"BND", "TRA", "CTX", "INV", "DUP", "CNV", "COMPLEX_INDEL"}
    ranked = []
    control_pool = []
    for sample, sample_rows in sorted(by_sample.items()):
        sample_rows.sort(key=lambda row: row["haplotype"])
        sample_evidence_counts = evidence_counts[sample]
        sample_catalog = catalog_by_sample.get(sample, [])
        callability = {}
        for row in sample_rows:
            summary = screen.get(row["assembly_id"], {})
            callability[row["haplotype"]] = max(
                fraction(summary.get("primary_callable_fraction")),
                fraction(summary.get("sensitivity_callable_fraction")),
            )
        uncallable = sum(value < args.callable_threshold for value in callability.values())
        low_tier = any(
            row["graph_member"] == "false" and row["rescue_tier"] in low_tiers
            for row in sample_rows
        )
        missing_mate = any(
            row["graph_context"] == "restore_missing_mate_for_graph_sample"
            for row in sample_rows
        )
        validation_candidate = any(
            row.get("validation_status") == "PENDING_HIFI"
            or row.get("independent_sample_count") == "1"
            or row.get("svtype") in difficult_types
            for row in sample_catalog
        )
        uncertain_genotype = bool(sample_catalog) and not validation_candidate

        reasons = []
        intended = []
        if low_tier or uncallable:
            reasons.append("P1_DISCOVERY_BLIND_SPOT")
            intended.append("discovery")
        if validation_candidate:
            reasons.append("P1_VALIDATE_CANDIDATE")
            intended.append("validation")
        if uncertain_genotype:
            reasons.append("P2_GENOTYPE_OR_PHASE")
            intended.extend(["genotyping", "phasing"])
        if missing_mate:
            reasons.append("RESTORE_MISSING_GRAPH_MATE")

        if reasons:
            if "P1_DISCOVERY_BLIND_SPOT" in reasons:
                tier = "P1_DISCOVERY_BLIND_SPOT"
            elif "P1_VALIDATE_CANDIDATE" in reasons:
                tier = "P1_VALIDATE_CANDIDATE"
            else:
                tier = "P2_GENOTYPE_OR_PHASE"
        else:
            tier = "P3_CONTROL"
            intended = ["control"]

        methods = sample_evidence_counts["method"]
        confidence = sample_evidence_counts["confidence"]
        svtypes = sample_evidence_counts["svtype"]
        important = sorted(
            sample_catalog,
            key=lambda row: (
                row.get("validation_status") != "PENDING_HIFI",
                row.get("confidence") != "HIGH",
                row["event_id"],
            ),
        )
        first = sample_rows[0]
        recommendation = {
            "P1_DISCOVERY_BLIND_SPOT": "Retrieve HiFi reads for sensitive read-based discovery; assembly no-calls are not negative evidence.",
            "P1_VALIDATE_CANDIDATE": "Retrieve HiFi reads for breakpoint-spanning validation of linked candidates.",
            "P2_GENOTYPE_OR_PHASE": "Retrieve HiFi reads to resolve genotype, breakpoint, or haplotype assignment.",
            "P3_CONTROL": "Retrieve or reuse HiFi reads as an assembly/read concordance control.",
        }[tier]
        output_row = {
            "sample_id": sample,
            "assembly_ids": ";".join(row["assembly_id"] for row in sample_rows),
            "graph_membership": ";".join(f"{row['haplotype']}={row['graph_member']}" for row in sample_rows),
            "missing_mate_status": ";".join(sorted({row["mate_status"] for row in sample_rows})),
            "rescue_tiers": ";".join(f"{row['haplotype']}={row['rescue_tier']}" for row in sample_rows),
            "callable_fractions": ";".join(f"{hap}={callability[hap]:.6f}" for hap in sorted(callability)),
            "uncallable_haplotype_count": uncallable,
            "priority_tier": tier,
            "reason_codes": ";".join(dict.fromkeys(reasons or ["P3_CONTROL"])),
            "recommendation": recommendation,
            "candidate_count": len(sample_catalog),
            "candidate_counts_by_method": counter_text(methods),
            "candidate_counts_by_confidence": counter_text(confidence),
            "candidate_counts_by_svtype": counter_text(svtypes),
            "important_candidate_ids": ";".join(row["event_id"] for row in important[:25]),
            "population": first.get("population", ""),
            "sequencing_batch": first.get("sequencing_batch", ""),
            "intended_use": ";".join(dict.fromkeys(intended)),
            "requested_read_type": "HiFi",
            "read_access_status": first.get("read_access_status", "not_checked"),
            "path_or_accession": first.get("hifi_path_or_accession", ""),
            "platform": first.get("read_platform", ""),
            "coverage": first.get("read_coverage", ""),
            "checksum": first.get("read_checksum", ""),
            "__min_callable": min(callability.values()) if callability else 0.0,
        }
        if tier == "P3_CONTROL":
            control_pool.append(output_row)
        else:
            ranked.append(output_row)

    # Controls are deliberately small and deterministic.  Prefer graph-member
    # samples and alternate positive-candidate/no-candidate controls.
    control_pool.sort(
        key=lambda row: (
            "=true" not in row["graph_membership"],
            row["candidate_count"] == 0,
            row["sample_id"],
        )
    )
    ranked.extend(control_pool[: args.control_count])
    tier_order = {
        "P1_DISCOVERY_BLIND_SPOT": 0,
        "P1_VALIDATE_CANDIDATE": 1,
        "P2_GENOTYPE_OR_PHASE": 2,
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
        row.pop("__min_callable", None)
    write_tsv(args.output, ranked)


if __name__ == "__main__":
    main()
