#!/usr/bin/env python3
"""Normalize and conservatively cluster assembly-only novel-SV evidence."""

import argparse
import csv
import gzip
import hashlib
import os
import re
import sqlite3
import tempfile
from collections import defaultdict
from pathlib import Path


EVIDENCE_FIELDS = [
    "evidence_id",
    "source_event_id",
    "coordinate_system",
    "chromosome",
    "sex_chromosome_status",
    "position_0",
    "end_0",
    "svtype",
    "svlen",
    "alternate_sequence",
    "discovery_method",
    "caller",
    "sample_id",
    "assembly_id",
    "haplotype",
    "graph_member",
    "assembly_qc_tier",
    "source_path",
    "source_filter",
    "graph_segment",
    "segment_offset_0",
    "confidence_tier",
]

CATALOG_FIELDS = [
    "event_id",
    "coordinate_system",
    "chromosome",
    "sex_chromosome_status",
    "position_0",
    "end_0",
    "svtype",
    "svlen",
    "reference_allele",
    "alternate_allele",
    "graph_segment",
    "segment_offset_0",
    "carrier_assemblies",
    "carrier_haplotypes",
    "carrier_samples",
    "carrier_assembly_count",
    "independent_sample_count",
    "discovery_methods",
    "caller_support",
    "assembly_qc_tiers",
    "graph_representation_status",
    "validation_status",
    "confidence",
    "evidence_ids",
]


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def read_tsv(path):
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def iter_tsv(path):
    with open_text(path, "rt") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def write_tsv(path, rows, fields):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_info(value):
    result = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            result[key] = val
        else:
            result[item] = True
    return result


def int_value(value, default=0):
    try:
        return int(float(str(value).split(",")[0]))
    except (TypeError, ValueError):
        return default


def infer_svtype(ref, alt, info):
    if info.get("SVTYPE"):
        return str(info["SVTYPE"]).upper()
    if "[" in alt or "]" in alt:
        return "BND"
    if alt.startswith("<") and alt.endswith(">"):
        return alt[1:-1].split(":", 1)[0].upper()
    delta = len(alt) - len(ref)
    return "INS" if delta > 0 else "DEL" if delta < 0 else "SUB"


def infer_svlen(ref, alt, info, position_0):
    if info.get("SVLEN") not in {None, ""}:
        return abs(int_value(info["SVLEN"]))
    if info.get("END") not in {None, ""}:
        return abs(int_value(info["END"]) - position_0 - 1)
    if not (alt.startswith("<") or "[" in alt or "]" in alt):
        return abs(len(alt) - len(ref))
    return 0


def alternate_sequence(ref, alt, svtype):
    if svtype == "INS" and not alt.startswith("<") and "[" not in alt and "]" not in alt:
        return alt[len(ref) :] if alt.startswith(ref) else alt
    return ""


def evidence_id(parts):
    return "EVD_" + hashlib.sha1("|".join(map(str, parts)).encode()).hexdigest()[:16]


def sex_chromosome_status(coordinate_system, chromosome):
    if coordinate_system == "GRAPH":
        return "UNRESOLVED_GRAPH_ORIGIN"
    if re.fullmatch(r"(?:chr)?[XY]", str(chromosome), flags=re.IGNORECASE):
        return "KNOWN_SEX_CHROMOSOME"
    return "NOT_KNOWN_SEX_CHROMOSOME"


def parse_vcf(source, min_sv_size, excluded_contig=None):
    source_regex = source.get("exclude_contig_regex", "")
    source_excluded_contig = re.compile(source_regex) if source_regex else excluded_contig
    with open_text(source["path"], "rt") as handle:
        samples = []
        for line_number, line in enumerate(handle, start=1):
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                raise ValueError(f"Malformed VCF row at {source['path']}:{line_number}")
            chrom, pos, variant_id, ref, alts, qual, filters, info_text = fields[:8]
            if source_excluded_contig and source_excluded_contig.fullmatch(chrom):
                continue
            for alt_index, alt in enumerate(alts.split(","), start=1):
                info = parse_info(info_text)
                svtype = infer_svtype(ref, alt, info)
                position_0 = int(pos) - 1
                svlen = infer_svlen(ref, alt, info, position_0)
                if svtype not in {"BND", "TRA", "CTX", "INV"} and svlen < min_sv_size:
                    continue
                end_0 = int_value(info.get("END"), position_0 + max(1, svlen))
                source_id = variant_id if variant_id not in {"", "."} else f"{chrom}:{pos}:{ref}:{alt_index}"
                row = {
                    "source_event_id": source_id,
                    "coordinate_system": source["coordinate_system"],
                    "chromosome": chrom,
                    "sex_chromosome_status": sex_chromosome_status(
                        source["coordinate_system"], chrom
                    ),
                    "position_0": position_0,
                    "end_0": end_0,
                    "svtype": svtype,
                    "svlen": svlen,
                    "alternate_sequence": alternate_sequence(ref, alt, svtype),
                    "discovery_method": f"assembly_{source['caller']}",
                    "caller": source["caller"],
                    "sample_id": source.get("sample_id", ""),
                    "assembly_id": source.get("assembly_id", ""),
                    "haplotype": source.get("haplotype", ""),
                    "source_path": source["path"],
                    "source_filter": filters,
                    "graph_segment": "",
                    "segment_offset_0": "",
                    "confidence_tier": "CALLER_PASS" if filters in {"PASS", "."} else "CALLER_FILTERED_REVIEW",
                }
                row["evidence_id"] = evidence_id(
                    [source["caller"], source["path"], source_id, alt_index]
                )
                yield row


def parse_graph_candidates(path, min_sv_size, excluded_contig=None):
    for source in iter_tsv(path):
        size = int_value(source.get("event_size_bp"))
        if size < min_sv_size:
            continue
        stable_source = source.get("stable_source", "")
        source_rank = source.get("source_rank", "")
        if stable_source and str(source_rank) == "0":
            coordinate_system = "CHM13"
            chromosome = source.get("chromosome", stable_source)
            position = int_value(source.get("stable_position_0"))
        else:
            coordinate_system = "GRAPH"
            chromosome = source.get("graph_segment", "unprojected")
            position = int_value(source.get("segment_offset_0"))
        filter_labels = {
            chromosome,
            source.get("chromosome", ""),
            stable_source,
        }
        if excluded_contig and any(
            excluded_contig.fullmatch(label) for label in filter_labels if label
        ):
            continue
        svtype = source.get("svtype", "")
        svlen = abs(int_value(source.get("svlen"), size))
        row = {
            "source_event_id": source.get("event_id", ""),
            "coordinate_system": coordinate_system,
            "chromosome": chromosome,
            "sex_chromosome_status": sex_chromosome_status(
                coordinate_system, chromosome
            ),
            "position_0": position,
            "end_0": position + (svlen if svtype == "DEL" else 1),
            "svtype": svtype,
            "svlen": svlen,
            "alternate_sequence": "",
            "discovery_method": "graph_residual",
            "caller": "minigraph_cigar",
            "sample_id": source.get("sample_id", ""),
            "assembly_id": source.get("assembly_id", ""),
            "haplotype": source.get("haplotype", ""),
            "source_path": path,
            "source_filter": source.get("filter_reasons", ""),
            "graph_segment": source.get("graph_segment", ""),
            "segment_offset_0": source.get("segment_offset_0", ""),
            "confidence_tier": source.get("confidence_tier", "REVIEW"),
        }
        row["evidence_id"] = evidence_id(["graph", row["source_event_id"]])
        yield row


def length_compatible(left, right, similarity):
    if left == 0 or right == 0:
        return True
    return min(left, right) / max(left, right) >= similarity


def compatible(event, representative, distance, similarity):
    if event["coordinate_system"] != representative["coordinate_system"]:
        return False
    if event["chromosome"] != representative["chromosome"]:
        return False
    if event["svtype"] != representative["svtype"]:
        return False
    if abs(event["position_0"] - representative["position_0"]) > distance:
        return False
    if not length_compatible(event["svlen"], representative["svlen"], similarity):
        return False
    left_sequence = event.get("alternate_sequence", "")
    right_sequence = representative.get("alternate_sequence", "")
    if left_sequence and right_sequence and left_sequence != right_sequence:
        # Sequence-resolved insertions are kept separate unless their exact
        # allele agrees.  This intentionally avoids distance-only overmerging.
        return False
    return True


def cluster_evidence(evidence, distance, similarity):
    clusters = []
    ordered = sorted(
        evidence,
        key=lambda row: (
            row["coordinate_system"],
            row["chromosome"],
            row["position_0"],
            row["svtype"],
            row["svlen"],
            row["evidence_id"],
        ),
    )
    for event in ordered:
        match = None
        for cluster in reversed(clusters):
            representative = cluster[0]
            if representative["coordinate_system"] != event["coordinate_system"] or representative["chromosome"] != event["chromosome"]:
                continue
            if event["position_0"] - representative["position_0"] > distance:
                break
            if compatible(event, representative, distance, similarity):
                match = cluster
                break
        if match is None:
            clusters.append([event])
        else:
            match.append(event)
    return clusters


def cluster_sorted_evidence(evidence, distance, similarity):
    """Cluster coordinate-sorted evidence with bounded in-memory state."""
    active = []
    active_key = None
    for event in evidence:
        key = (event["coordinate_system"], event["chromosome"])
        if key != active_key:
            for cluster in active:
                yield cluster
            active = []
            active_key = key

        keep = []
        for cluster in active:
            if event["position_0"] - cluster[0]["position_0"] > distance:
                yield cluster
            else:
                keep.append(cluster)
        active = keep

        match = None
        for cluster in reversed(active):
            if compatible(event, cluster[0], distance, similarity):
                match = cluster
                break
        if match is None:
            active.append([event])
        else:
            match.append(event)
    for cluster in active:
        yield cluster


def unique_join(values):
    return ";".join(sorted({str(value) for value in values if value not in {"", None}}))


def catalog_row(cluster, index):
    representative = sorted(
        cluster,
        key=lambda row: (
            row["confidence_tier"] not in {"HIGH_CONFIDENCE", "CALLER_PASS"},
            -row["svlen"],
            row["evidence_id"],
        ),
    )[0]
    methods = {row["discovery_method"] for row in cluster}
    callers = {row["caller"] for row in cluster}
    samples = {row["sample_id"] for row in cluster if row["sample_id"]}
    assemblies = {row["assembly_id"] for row in cluster if row["assembly_id"]}
    qc_tiers = {row.get("assembly_qc_tier", "") for row in cluster if row.get("assembly_qc_tier")}
    high_evidence = any(row["confidence_tier"] in {"HIGH_CONFIDENCE", "CALLER_PASS"} for row in cluster)
    independent_methods = len(methods) >= 2 or len(callers) >= 2
    if high_evidence and (independent_methods or len(samples) >= 2):
        confidence = "HIGH"
    elif high_evidence:
        confidence = "MEDIUM"
    else:
        confidence = "UNCERTAIN"
    poor_only = bool(qc_tiers) and qc_tiers <= {"fragmented_rescue", "not_recommended"}
    validation = "PENDING_HIFI" if poor_only and len(samples) < 2 and not independent_methods else "ASSEMBLY_ONLY_REVIEW"
    key = (
        f"{representative['coordinate_system']}|{representative['chromosome']}|"
        f"{representative['position_0']}|{representative['svtype']}|{representative['svlen']}|{index}"
    )
    event_id = "PSV_" + hashlib.sha1(key.encode()).hexdigest()[:16]
    sequence = representative.get("alternate_sequence", "")
    return {
        "event_id": event_id,
        "coordinate_system": representative["coordinate_system"],
        "chromosome": representative["chromosome"],
        "sex_chromosome_status": unique_join(
            row.get("sex_chromosome_status", "") for row in cluster
        ),
        "position_0": representative["position_0"],
        "end_0": representative["end_0"],
        "svtype": representative["svtype"],
        "svlen": representative["svlen"],
        "reference_allele": "",
        "alternate_allele": sequence,
        "graph_segment": representative.get("graph_segment", ""),
        "segment_offset_0": representative.get("segment_offset_0", ""),
        "carrier_assemblies": unique_join(row["assembly_id"] for row in cluster),
        "carrier_haplotypes": unique_join(
            f"{row['sample_id']}:{row['haplotype']}" for row in cluster if row["sample_id"] and row["haplotype"]
        ),
        "carrier_samples": unique_join(samples),
        "carrier_assembly_count": len(assemblies),
        "independent_sample_count": len(samples),
        "discovery_methods": unique_join(methods),
        "caller_support": unique_join(callers),
        "assembly_qc_tiers": unique_join(qc_tiers),
        "graph_representation_status": "CANDIDATE_MISSING_ALLELE",
        "validation_status": validation,
        "confidence": confidence,
        "evidence_ids": unique_join(row["evidence_id"] for row in cluster),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--graph-candidates", required=True)
    parser.add_argument("--call-manifest", required=True)
    parser.add_argument("--assembly-manifest", required=True)
    parser.add_argument("--catalog-output", required=True)
    parser.add_argument("--evidence-output", required=True)
    parser.add_argument("--min-sv-size", type=int, default=50)
    parser.add_argument("--breakpoint-distance", type=int, default=500)
    parser.add_argument("--length-similarity", type=float, default=0.70)
    parser.add_argument(
        "--exclude-contig-regex",
        default="",
        help=(
            "Full-match regular expression for contigs excluded from the catalog. "
            "This default applies to graph candidates and to call-manifest rows "
            "without their own exclude_contig_regex."
        ),
    )
    parser.add_argument(
        "--temp-dir",
        help="Directory for the disk-backed normalization database; defaults to TMPDIR",
    )
    args = parser.parse_args()
    excluded_contig = (
        re.compile(args.exclude_contig_regex) if args.exclude_contig_regex else None
    )

    manifest_rows = read_tsv(args.assembly_manifest)
    assembly_metadata = {row["assembly_id"]: row for row in manifest_rows}
    sample_metadata = defaultdict(list)
    for row in manifest_rows:
        sample_metadata[row["sample_id"]].append(row)
    def decorate(row):
        metadata = assembly_metadata.get(row["assembly_id"], {})
        if metadata:
            row["graph_member"] = metadata.get("graph_member", "")
            row["assembly_qc_tier"] = metadata.get("rescue_tier", "")
        else:
            sample_rows = sample_metadata.get(row["sample_id"], [])
            tiers = {item.get("rescue_tier", "") for item in sample_rows}
            if tiers and tiers <= {"fragmented_rescue", "not_recommended"}:
                row["assembly_qc_tier"] = (
                    "not_recommended" if "not_recommended" in tiers else "fragmented_rescue"
                )
            elif len(tiers) == 1:
                row["assembly_qc_tier"] = next(iter(tiers))
            else:
                row["assembly_qc_tier"] = "sample_mixed_qc" if tiers else ""
            row["graph_member"] = str(
                any(item.get("graph_member") == "true" for item in sample_rows)
            ).lower() if sample_rows else ""
        row["position_0"] = int_value(row["position_0"])
        row["end_0"] = int_value(row["end_0"])
        row["svlen"] = abs(int_value(row["svlen"]))
        return row

    temp_root = args.temp_dir or os.environ.get("TMPDIR") or str(Path(args.catalog_output).parent)
    Path(temp_root).mkdir(parents=True, exist_ok=True)
    descriptor, database_path = tempfile.mkstemp(prefix="novel-sv-evidence.", suffix=".sqlite", dir=temp_root)
    os.close(descriptor)
    connection = sqlite3.connect(database_path)
    try:
        columns = ", ".join(f'"{field}" TEXT' for field in EVIDENCE_FIELDS)
        connection.execute(f"CREATE TABLE evidence ({columns})")
        placeholders = ",".join("?" for _ in EVIDENCE_FIELDS)
        insert_sql = f"INSERT INTO evidence VALUES ({placeholders})"

        evidence_path = Path(args.evidence_output)
        evidence_path.parent.mkdir(parents=True, exist_ok=True)
        opener = gzip.open if str(evidence_path).endswith(".gz") else open
        inserted = 0
        with opener(evidence_path, "wt", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=EVIDENCE_FIELDS,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="ignore",
            )
            writer.writeheader()

            streams = [
                parse_graph_candidates(
                    args.graph_candidates, args.min_sv_size, excluded_contig
                )
            ]
            streams.extend(
                parse_vcf(source, args.min_sv_size, excluded_contig)
                for source in read_tsv(args.call_manifest)
            )
            for stream in streams:
                for raw_row in stream:
                    row = decorate(raw_row)
                    writer.writerow(row)
                    connection.execute(
                        insert_sql,
                        [str(row.get(field, "")) for field in EVIDENCE_FIELDS],
                    )
                    inserted += 1
                    if inserted % 10000 == 0:
                        connection.commit()
        connection.commit()
        connection.execute(
            "CREATE INDEX evidence_order ON evidence "
            '(coordinate_system, chromosome, CAST(position_0 AS INTEGER), svtype, CAST(svlen AS INTEGER), evidence_id)'
        )
        connection.commit()

        select_fields = ", ".join(f'"{field}"' for field in EVIDENCE_FIELDS)
        cursor = connection.execute(
            f"SELECT {select_fields} FROM evidence ORDER BY "
            "coordinate_system, chromosome, CAST(position_0 AS INTEGER), "
            "svtype, CAST(svlen AS INTEGER), evidence_id"
        )

        def ordered_evidence():
            for values in cursor:
                row = dict(zip(EVIDENCE_FIELDS, values))
                row["position_0"] = int_value(row["position_0"])
                row["end_0"] = int_value(row["end_0"])
                row["svlen"] = abs(int_value(row["svlen"]))
                yield row

        clusters = cluster_sorted_evidence(
            ordered_evidence(), args.breakpoint_distance, args.length_similarity
        )
        catalog = (
            catalog_row(cluster, index)
            for index, cluster in enumerate(clusters, start=1)
        )
        write_tsv(args.catalog_output, catalog, CATALOG_FIELDS)
    finally:
        connection.close()
        Path(database_path).unlink(missing_ok=True)


if __name__ == "__main__":
    main()
