#!/usr/bin/env python3
"""Build read-derived, two-sided local rescue sequences for graph SV alleles.

The script is deliberately separate from graph construction.  It prepares a
reviewable plan, runs one local HiFi assembly at a time, and combines only
passing rescue sequences.  Coordinates in input and output tables are 0-based
and half-open unless a field explicitly says otherwise.
"""

import argparse
import csv
import gzip
import hashlib
import heapq
import json
import re
import shlex
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


PLAN_FIELDS = [
    "rescue_id",
    "plan_status",
    "status_reason",
    "catalog_event_id",
    "coordinate_system",
    "chromosome",
    "position_0",
    "svtype",
    "svlen",
    "catalog_confidence",
    "validation_status",
    "independent_sample_count",
    "linear_supporting_sample_count",
    "discovery_methods",
    "source_event_id",
    "sample_id",
    "assembly_id",
    "haplotype",
    "source_qc_tier",
    "assembly_path",
    "hifi_read_paths",
    "query_name",
    "query_length",
    "query_position_0",
    "query_event_span_bp",
    "seed_start_0",
    "seed_end_0",
    "left_seed_flank_bp",
    "right_seed_flank_bp",
    "source_confidence_tier",
    "source_filter_reasons",
    "source_left_anchor_bp",
    "source_right_anchor_bp",
    "source_identity",
    "source_mapq",
    "source_mapping_pass",
    "graph_segment",
    "segment_offset_0",
]

QC_FIELDS = [
    "rescue_id",
    "status",
    "failure_reasons",
    "catalog_event_id",
    "source_event_id",
    "sample_id",
    "assembly_id",
    "haplotype",
    "svtype",
    "svlen",
    "seed_length_bp",
    "selected_read_count",
    "spanning_read_count",
    "source_allele_spanning_read_count",
    "local_contig_count",
    "selected_local_contig",
    "selected_alignment_strand",
    "selected_alignment_mapq",
    "selected_alignment_identity",
    "selected_target_coverage",
    "selected_event_discordance_bp",
    "rescue_length_bp",
    "rescue_n_count",
    "rescue_sha256",
    "clean_rescue_fasta",
]

SUMMARY_FIELDS = QC_FIELDS + [
    "combined_fasta_status",
    "combined_sequence_name",
]

CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
SAFE_NAME_RE = re.compile(r"[^A-Za-z0-9_.-]+")


def allow_large_tsv_fields():
    """Raise csv's field limit for catalog carrier/evidence fields."""
    limit = sys.maxsize
    while True:
        try:
            csv.field_size_limit(limit)
            return
        except OverflowError:
            limit //= 10


def open_text(path, mode="rt"):
    path = Path(path)
    return gzip.open(path, mode) if str(path).endswith(".gz") else path.open(mode)


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


def split_values(value):
    return [item for item in str(value or "").split(";") if item]


def split_paths(value):
    return [
        item.strip() for item in re.split(r"[;,]", str(value or "")) if item.strip()
    ]


def int_value(value, default=0):
    try:
        return int(float(str(value)))
    except (TypeError, ValueError):
        return default


def float_value(value, default=0.0):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def true_value(value):
    return str(value).lower() in {"1", "true", "yes"}


def safe_name(value):
    name = SAFE_NAME_RE.sub("_", str(value)).strip("._")
    if not name:
        raise ValueError(f"Cannot make a safe name from {value!r}")
    return name


def first(row, *fields, default=""):
    for field in fields:
        value = row.get(field, "")
        if value not in {"", None}:
            return value
    return default


def sequence_sha256(sequence):
    return hashlib.sha256(sequence.upper().encode()).hexdigest()


def read_event_ids(path):
    if not path:
        return set()
    result = set()
    with Path(path).open() as handle:
        for line in handle:
            value = line.strip().split("\t", 1)[0]
            if value and not value.startswith("#") and value != "event_id":
                result.add(value)
    if not result:
        raise ValueError(f"No event IDs found in {path}")
    return result


def slim_catalog_row(row):
    return {
        "event_id": row["event_id"],
        "coordinate_system": row["coordinate_system"],
        "chromosome": row["chromosome"],
        "position_0": int_value(row["position_0"]),
        "svtype": row["svtype"],
        "svlen": abs(int_value(row["svlen"])),
        "confidence": row.get("confidence", ""),
        "validation_status": row.get("validation_status", ""),
        "independent_sample_count": int_value(row.get("independent_sample_count")),
        "linear_supporting_sample_count": int_value(
            row.get("linear_supporting_sample_count")
        ),
        "discovery_methods": row.get("discovery_methods", ""),
        "carrier_assemblies": set(split_values(row.get("carrier_assemblies", ""))),
    }


def catalog_priority(row):
    """Rank recurrent, independently supported, larger events first."""
    return (
        row["independent_sample_count"],
        row["linear_supporting_sample_count"],
        row["svlen"],
    )


def select_catalog_events(args):
    explicit_ids = read_event_ids(args.event_ids)
    required_methods = set(split_values(args.require_methods.replace(",", ";")))
    selected_by_id = {}
    heap = []
    seen_ids = set()

    for raw in iter_tsv(args.catalog):
        event_id = raw.get("event_id", "")
        if explicit_ids:
            if event_id not in explicit_ids:
                continue
        else:
            methods = set(split_values(raw.get("discovery_methods", "")))
            if args.confidence and raw.get("confidence") != args.confidence:
                continue
            if (
                args.validation_status
                and raw.get("validation_status") != args.validation_status
            ):
                continue
            if not required_methods.issubset(methods):
                continue
            if (
                int_value(raw.get("independent_sample_count"))
                < args.min_independent_samples
            ):
                continue

        row = slim_catalog_row(raw)
        seen_ids.add(event_id)
        if explicit_ids:
            selected_by_id[event_id] = row
            continue

        score = catalog_priority(row)
        item = (score, event_id, row)
        if len(heap) < args.max_events:
            heapq.heappush(heap, item)
        elif item[:2] > heap[0][:2]:
            heapq.heapreplace(heap, item)

    if explicit_ids:
        missing = sorted(explicit_ids - seen_ids)
        if missing:
            preview = ", ".join(missing[:10])
            suffix = " ..." if len(missing) > 10 else ""
            raise ValueError(f"Catalog event IDs not found: {preview}{suffix}")
        selected = list(selected_by_id.values())
    else:
        selected = [item[2] for item in heap]

    selected.sort(
        key=lambda row: (
            -row["independent_sample_count"],
            -row["linear_supporting_sample_count"],
            -row["svlen"],
            row["event_id"],
        )
    )
    if not selected:
        raise ValueError(
            "No catalog events passed the planning filters. Supply "
            "--event-ids or relax the confidence, validation, method, or "
            "sample-support filters."
        )
    return selected


def candidate_coordinate(row):
    stable_source = row.get("stable_source", "")
    source_rank = row.get("source_rank", "")
    graph_segment = row.get("graph_segment", "")
    segment_offset = int_value(row.get("segment_offset_0"))
    stable_position = int_value(row.get("stable_position_0"))
    if stable_source and str(source_rank) == "0":
        return (
            "CHM13",
            row.get("chromosome") or stable_source,
            stable_position,
        )
    chromosome = (
        f"{stable_source}|SR{source_rank or 'UNKNOWN'}"
        if stable_source
        else graph_segment
    )
    position = stable_position if stable_source else segment_offset
    return "GRAPH", chromosome, position


def length_compatible(left, right, similarity):
    if left == 0 or right == 0:
        return True
    return min(left, right) / max(left, right) >= similarity


def source_rank(row, event, assembly_metadata):
    metadata = assembly_metadata.get(row.get("assembly_id", ""), {})
    tier_order = {
        "best_rescue": 4,
        "reasonable_rescue": 3,
        "fragmented_rescue": 2,
        "not_recommended": 1,
    }
    position = candidate_coordinate(row)[2]
    size = abs(int_value(row.get("svlen"), int_value(row.get("event_size_bp"))))
    return (
        row.get("confidence_tier") == "HIGH_CONFIDENCE",
        not bool(row.get("filter_reasons")),
        true_value(row.get("primary")),
        min(
            int_value(row.get("left_anchor_bp")),
            int_value(row.get("right_anchor_bp")),
        ),
        int_value(row.get("mapq")),
        float_value(row.get("identity")),
        tier_order.get(metadata.get("rescue_tier", ""), 0),
        -abs(position - event["position_0"]),
        -abs(size - event["svlen"]),
    )


def choose_source_candidates(
    events,
    candidate_path,
    assembly_metadata,
    allowed_tiers,
    breakpoint_distance,
    length_similarity,
):
    event_index = defaultdict(list)
    for event in events:
        event_index[
            (
                event["coordinate_system"],
                event["chromosome"],
                event["svtype"],
            )
        ].append(event)

    best = {}
    for row in iter_tsv(candidate_path):
        if true_value(row.get("graph_member")):
            continue
        assembly_id = row.get("assembly_id", "")
        metadata = assembly_metadata.get(assembly_id, {})
        tier = metadata.get("rescue_tier", "")
        if allowed_tiers is not None and tier not in allowed_tiers:
            continue

        coordinate_system, chromosome, position = candidate_coordinate(row)
        key = (coordinate_system, chromosome, row.get("svtype", ""))
        possible = event_index.get(key, ())
        if not possible:
            continue
        size = abs(int_value(row.get("svlen"), int_value(row.get("event_size_bp"))))
        for event in possible:
            if assembly_id not in event["carrier_assemblies"]:
                continue
            if abs(position - event["position_0"]) > breakpoint_distance:
                continue
            if not length_compatible(size, event["svlen"], length_similarity):
                continue
            rank = source_rank(row, event, assembly_metadata)
            previous = best.get(event["event_id"])
            if (
                previous is None
                or rank > previous[0]
                or (
                    rank == previous[0]
                    and row.get("event_id", "") < previous[1].get("event_id", "")
                )
            ):
                best[event["event_id"]] = (rank, dict(row))
    return {event_id: value[1] for event_id, value in best.items()}


def index_assemblies(path):
    result = {}
    for row in iter_tsv(path):
        assembly_id = row.get("assembly_id", "")
        if not assembly_id:
            raise ValueError(f"Missing assembly_id in {path}")
        if assembly_id in result:
            raise ValueError(f"Duplicate assembly_id in {path}: {assembly_id}")
        result[assembly_id] = row
    return result


def index_reads(path):
    if not path:
        return {}
    result = {}
    for row in iter_tsv(path):
        sample = row.get("sample_id") or row.get("sample")
        if not sample:
            raise ValueError(f"Missing sample_id in {path}")
        if sample in result:
            raise ValueError(f"Duplicate sample_id in {path}: {sample}")
        result[sample] = row
    return result


def read_paths_for_sample(sample, read_metadata, assembly_row):
    row = read_metadata.get(sample, {})
    value = first(
        row,
        "hifi_path",
        "hifi_paths",
        "path_or_accession",
        "hifi_path_or_accession",
        default=first(
            assembly_row,
            "hifi_path_or_accession",
            "hifi_path",
        ),
    )
    return ";".join(split_paths(value))


def query_event_span(row):
    svtype = row.get("svtype", "").upper()
    if svtype in {"INS", "COMPLEX_INDEL"}:
        return max(
            1,
            abs(
                int_value(
                    row.get("event_size_bp"),
                    int_value(row.get("svlen"), 1),
                )
            ),
        )
    return 1


def build_plan(args):
    allow_large_tsv_fields()
    events = select_catalog_events(args)
    assemblies = index_assemblies(args.assembly_manifest)
    reads = index_reads(args.read_manifest)
    allowed_tiers = None
    if args.source_qc_tiers.lower() != "all":
        allowed_tiers = set(split_values(args.source_qc_tiers.replace(",", ";")))
        if not allowed_tiers:
            raise ValueError("--source-qc-tiers resolved to an empty set")

    sources = choose_source_candidates(
        events,
        args.graph_candidates,
        assemblies,
        allowed_tiers,
        args.breakpoint_distance,
        args.length_similarity,
    )

    plan = []
    for event in events:
        source = sources.get(event["event_id"])
        base = {
            "rescue_id": safe_name(event["event_id"]),
            "plan_status": "READY",
            "status_reason": "",
            "catalog_event_id": event["event_id"],
            "coordinate_system": event["coordinate_system"],
            "chromosome": event["chromosome"],
            "position_0": event["position_0"],
            "svtype": event["svtype"],
            "svlen": event["svlen"],
            "catalog_confidence": event["confidence"],
            "validation_status": event["validation_status"],
            "independent_sample_count": event["independent_sample_count"],
            "linear_supporting_sample_count": event["linear_supporting_sample_count"],
            "discovery_methods": event["discovery_methods"],
        }
        if source is None:
            base.update(
                {
                    "plan_status": "NO_MATCHING_SOURCE",
                    "status_reason": (
                        "No compatible carrier graph-residual row was found "
                        "in the requested source QC tiers"
                    ),
                }
            )
            plan.append(base)
            continue

        assembly_id = source["assembly_id"]
        assembly = assemblies.get(assembly_id, {})
        sample = source.get("sample_id") or assembly.get("sample_id", "")
        query_position = int_value(source.get("query_position_0"))
        query_length = int_value(source.get("query_length"))
        event_span = query_event_span(source)
        event_end = min(query_length, query_position + event_span)
        seed_start = max(0, query_position - args.flank_bp)
        seed_end = min(query_length, event_end + args.flank_bp)
        left_flank = query_position - seed_start
        right_flank = seed_end - event_end
        assembly_path = first(
            assembly,
            "cleaned_fasta_path",
            "path",
            "cleaned_path",
        )
        read_paths = read_paths_for_sample(sample, reads, assembly)

        base.update(
            {
                "source_event_id": source.get("event_id", ""),
                "sample_id": sample,
                "assembly_id": assembly_id,
                "haplotype": source.get("haplotype") or assembly.get("haplotype", ""),
                "source_qc_tier": assembly.get("rescue_tier", ""),
                "assembly_path": assembly_path,
                "hifi_read_paths": read_paths,
                "query_name": source.get("query_name", ""),
                "query_length": query_length,
                "query_position_0": query_position,
                "query_event_span_bp": event_span,
                "seed_start_0": seed_start,
                "seed_end_0": seed_end,
                "left_seed_flank_bp": left_flank,
                "right_seed_flank_bp": right_flank,
                "source_confidence_tier": source.get("confidence_tier", ""),
                "source_filter_reasons": source.get("filter_reasons", ""),
                "source_left_anchor_bp": source.get("left_anchor_bp", ""),
                "source_right_anchor_bp": source.get("right_anchor_bp", ""),
                "source_identity": source.get("identity", ""),
                "source_mapq": source.get("mapq", ""),
                "source_mapping_pass": source.get("mapping_pass", ""),
                "graph_segment": source.get("graph_segment", ""),
                "segment_offset_0": source.get("segment_offset_0", ""),
            }
        )
        if not assembly_path:
            base["plan_status"] = "MISSING_ASSEMBLY_PATH"
            base["status_reason"] = "Assembly manifest has no FASTA path"
        elif not read_paths:
            base["plan_status"] = "WAITING_FOR_HIFI"
            base[
                "status_reason"
            ] = "No local HiFi FASTA/FASTQ path is present for the sample"
        elif (
            min(left_flank, right_flank) < args.min_flank_bp
            and not args.allow_contig_end
        ):
            base["plan_status"] = "SOURCE_NEAR_CONTIG_END"
            base["status_reason"] = (
                f"Only {min(left_flank, right_flank)} bp is available on "
                f"the shorter side; require {args.min_flank_bp}"
            )
        plan.append(base)

    write_tsv(args.output, plan, PLAN_FIELDS)
    status_counts = defaultdict(int)
    for row in plan:
        status_counts[row["plan_status"]] += 1
    summary = ", ".join(
        f"{status}={status_counts[status]}" for status in sorted(status_counts)
    )
    print(f"Wrote {len(plan)} rescue plans to {args.output}: {summary}")


def require_executable(name):
    path = shutil.which(name)
    if not path:
        raise RuntimeError(
            f"Required executable is unavailable: {name}. Activate the "
            "sv_rescue conda environment."
        )
    return path


def run_command(command, stdout_path, stderr_path):
    with Path(stdout_path).open("w") as stdout, Path(stderr_path).open("w") as stderr:
        subprocess.run(command, check=True, stdout=stdout, stderr=stderr)


def tool_version(executable):
    for flag in ("--version", "-V"):
        result = subprocess.run(
            [executable, flag],
            capture_output=True,
            text=True,
            check=False,
        )
        value = (result.stdout + "\n" + result.stderr).strip()
        if value:
            return value.splitlines()[0]
    return "version unavailable"


def parse_fasta(path):
    sequences = {}
    name = None
    chunks = []
    with open_text(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    sequences[name] = "".join(chunks)
                name = line[1:].strip().split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        sequences[name] = "".join(chunks)
    return sequences


def write_fasta(path, records, wrap=80):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "wt") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), wrap):
                handle.write(sequence[start : start + wrap] + "\n")


def reverse_complement(sequence):
    table = str.maketrans(
        "ACGTRYMKBDHVNacgtrymkbdhvn",
        "TGCAYRKMVHDBNtgcayrkmvhdbn",
    )
    return sequence.translate(table)[::-1]


def extract_seed(plan, samtools, output, log):
    assembly = Path(plan["assembly_path"])
    if not assembly.is_file():
        raise FileNotFoundError(f"Assembly FASTA is unavailable: {assembly}")
    fai = Path(str(assembly) + ".fai")
    if not fai.is_file():
        raise FileNotFoundError(
            f"Assembly FASTA index is unavailable: {fai}. Create it once "
            f"on the execution server with: samtools faidx {assembly}"
        )
    start = int_value(plan["seed_start_0"])
    end = int_value(plan["seed_end_0"])
    region = f"{plan['query_name']}:{start + 1}-{end}"
    raw = Path(output).with_suffix(".samtools.fa")
    run_command([samtools, "faidx", str(assembly), region], raw, log)
    sequences = parse_fasta(raw)
    if len(sequences) != 1:
        raise RuntimeError(
            f"Expected one sequence from samtools faidx {region}, got "
            f"{len(sequences)}"
        )
    sequence = next(iter(sequences.values()))
    expected = end - start
    if len(sequence) != expected:
        raise RuntimeError(
            f"Seed length mismatch for {region}: {len(sequence)} != {expected}"
        )
    write_fasta(output, [("seed", sequence)])
    raw.unlink(missing_ok=True)
    return sequence


def parse_paf(path):
    rows = []
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed PAF line {line_number} in {path}")
            tags = {}
            for field in fields[12:]:
                parts = field.split(":", 2)
                if len(parts) == 3:
                    tags[parts[0]] = parts[2]
            rows.append(
                {
                    "query_name": fields[0],
                    "query_length": int(fields[1]),
                    "query_start": int(fields[2]),
                    "query_end": int(fields[3]),
                    "strand": fields[4],
                    "target_name": fields[5],
                    "target_length": int(fields[6]),
                    "target_start": int(fields[7]),
                    "target_end": int(fields[8]),
                    "matches": int(fields[9]),
                    "alignment_length": int(fields[10]),
                    "mapq": int(fields[11]),
                    "tags": tags,
                }
            )
    return rows


def read_selection_from_paf(
    path,
    event_start,
    event_end,
    min_mapq,
    min_aligned_bp,
    read_anchor_bp,
    max_event_discordance_bp,
):
    selected = set()
    spanning = set()
    source_allele_spanning = set()
    left_limit = max(0, event_start - read_anchor_bp)
    right_limit = event_end + read_anchor_bp
    for row in parse_paf(path):
        if row["target_name"] != "seed" or row["mapq"] < min_mapq:
            continue
        if row["tags"].get("tp") == "S":
            continue
        target_aligned = row["target_end"] - row["target_start"]
        if target_aligned < min_aligned_bp:
            continue
        selected.add(row["query_name"])
        if row["target_start"] <= left_limit and row["target_end"] >= right_limit:
            spanning.add(row["query_name"])
            discordance = cigar_event_discordance(
                row["tags"].get("cg", ""),
                row["target_start"],
                event_start,
                event_end,
                read_anchor_bp,
            )
            if discordance <= max_event_discordance_bp:
                source_allele_spanning.add(row["query_name"])
    return selected, spanning, source_allele_spanning


def iter_fastx(path):
    """Yield name and sequence from ordinary FASTA or four-line FASTQ."""
    with open_text(path, "rt") as handle:
        first_line = handle.readline()
        if not first_line:
            return
        if first_line.startswith(">"):
            name = first_line[1:].strip().split()[0]
            chunks = []
            for line in handle:
                if line.startswith(">"):
                    yield name, "".join(chunks)
                    name = line[1:].strip().split()[0]
                    chunks = []
                else:
                    chunks.append(line.strip())
            yield name, "".join(chunks)
            return
        if not first_line.startswith("@"):
            raise ValueError(f"Unrecognized FASTA/FASTQ input: {path}")

        header = first_line
        while header:
            name = header[1:].strip().split()[0]
            sequence = handle.readline().strip()
            plus = handle.readline()
            quality = handle.readline().strip()
            if not plus.startswith("+") or len(sequence) != len(quality):
                raise ValueError(
                    f"{path} is not an ordinary four-line FASTQ near {name}. "
                    "Convert it to standard FASTQ before rescue construction."
                )
            yield name, sequence
            header = handle.readline()
            if header and not header.startswith("@"):
                raise ValueError(f"Malformed FASTQ header in {path}: {header!r}")


def extract_selected_reads(read_paths, names, output):
    found = set()
    records = []
    for path in read_paths:
        for name, sequence in iter_fastx(path):
            if name in names and name not in found:
                found.add(name)
                records.append((name, sequence))
    missing = names - found
    if missing:
        preview = ", ".join(sorted(missing)[:5])
        raise RuntimeError(
            f"{len(missing)} mapped read names were not recovered from the "
            f"FASTA/FASTQ files, for example: {preview}"
        )
    write_fasta(output, records)
    return len(records)


def gfa_contigs(gfa_paths):
    records = []
    seen_sequences = set()
    for gfa in gfa_paths:
        prefix = safe_name(gfa.name.replace(".gfa", ""))
        with gfa.open() as handle:
            for line in handle:
                if not line.startswith("S\t"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3 or fields[2] == "*":
                    continue
                sequence = fields[2]
                digest = sequence_sha256(sequence)
                if digest in seen_sequences:
                    continue
                seen_sequences.add(digest)
                records.append((f"{prefix}|{safe_name(fields[1])}", sequence))
    return records


def cigar_event_discordance(cigar, target_start, event_start, event_end, window):
    """Count indel bases in an assembly alignment near the source event."""
    if not cigar:
        return sys.maxsize
    tokens = [(int(length), operation) for length, operation in CIGAR_RE.findall(cigar)]
    if not tokens or "".join(f"{n}{op}" for n, op in tokens) != cigar:
        raise ValueError(f"Malformed PAF CIGAR: {cigar}")
    target_position = target_start
    region_start = max(0, event_start - window)
    region_end = event_end + window
    discordance = 0
    for length, operation in tokens:
        if operation == "I":
            if region_start <= target_position <= region_end:
                discordance += length
        elif operation in {"D", "N"}:
            operation_end = target_position + length
            if operation_end >= region_start and target_position <= region_end:
                discordance += length
        if operation in {"M", "=", "X", "D", "N"}:
            target_position += length
    return discordance


def select_local_contig(
    paf_rows,
    event_start,
    event_end,
    seed_length,
    anchor_bp,
    min_mapq,
    min_identity,
    discordance_window,
):
    left_limit = max(0, event_start - anchor_bp)
    right_limit = min(seed_length, event_end + anchor_bp)
    candidates = []
    for row in paf_rows:
        if row["target_name"] != "seed":
            continue
        identity = (
            row["matches"] / row["alignment_length"] if row["alignment_length"] else 0.0
        )
        anchored = (
            row["target_start"] <= left_limit and row["target_end"] >= right_limit
        )
        target_coverage = (
            (row["target_end"] - row["target_start"]) / seed_length
            if seed_length
            else 0.0
        )
        discordance = cigar_event_discordance(
            row["tags"].get("cg", ""),
            row["target_start"],
            event_start,
            event_end,
            discordance_window,
        )
        row = dict(row)
        row.update(
            {
                "anchored": anchored,
                "identity": identity,
                "target_coverage": target_coverage,
                "event_discordance_bp": discordance,
            }
        )
        if anchored and row["mapq"] >= min_mapq and identity >= min_identity:
            candidates.append(row)
    if not candidates:
        return None
    candidates.sort(
        key=lambda row: (
            row["event_discordance_bp"],
            -row["target_coverage"],
            -row["identity"],
            -row["mapq"],
            -row["alignment_length"],
            row["query_name"],
        )
    )
    return candidates[0]


def load_plan_row(path, rescue_id):
    matches = [row for row in iter_tsv(path) if row.get("rescue_id") == rescue_id]
    if len(matches) != 1:
        raise ValueError(
            f"Expected one row for rescue_id={rescue_id!r} in {path}, "
            f"found {len(matches)}"
        )
    return matches[0]


def write_early_qc_failure(
    result_dir,
    status_path,
    plan,
    rescue_id,
    failure_reasons,
    **metrics,
):
    qc = {field: "" for field in QC_FIELDS}
    qc.update(
        {
            "rescue_id": rescue_id,
            "status": "FAIL",
            "failure_reasons": ";".join(failure_reasons),
            "catalog_event_id": plan["catalog_event_id"],
            "source_event_id": plan["source_event_id"],
            "sample_id": plan["sample_id"],
            "assembly_id": plan["assembly_id"],
            "haplotype": plan["haplotype"],
            "svtype": plan["svtype"],
            "svlen": plan["svlen"],
        }
    )
    qc.update(metrics)
    write_tsv(result_dir / "qc.tsv", [qc], QC_FIELDS)
    status_path.write_text(
        json.dumps(
            {
                "rescue_id": rescue_id,
                "status": "FAIL",
                "failure_reasons": failure_reasons,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    print(f"{rescue_id}: FAIL ({result_dir})")


def run_rescue(args):
    allow_large_tsv_fields()
    plan = load_plan_row(args.plan, args.rescue_id)
    if plan.get("plan_status") != "READY":
        raise ValueError(
            f"Rescue {args.rescue_id} is not READY: "
            f"{plan.get('plan_status')} {plan.get('status_reason', '')}"
        )

    minimap2 = require_executable("minimap2")
    samtools = require_executable("samtools")
    hifiasm = require_executable("hifiasm")

    result_dir = Path(args.output_dir) / safe_name(args.rescue_id)
    if result_dir.exists() and any(result_dir.iterdir()):
        raise FileExistsError(
            f"Result directory is not empty: {result_dir}. Use a new output "
            "directory so an earlier rescue is not overwritten."
        )
    result_dir.mkdir(parents=True, exist_ok=True)
    work = result_dir / "work"
    logs = result_dir / "logs"
    work.mkdir()
    logs.mkdir()
    status_path = result_dir / "run_status.json"

    versions = {
        "minimap2": tool_version(minimap2),
        "samtools": tool_version(samtools),
        "hifiasm": tool_version(hifiasm),
    }
    (result_dir / "tool_versions.json").write_text(
        json.dumps(versions, indent=2, sort_keys=True) + "\n"
    )

    try:
        read_paths = [Path(path) for path in split_paths(plan["hifi_read_paths"])]
        if not read_paths:
            raise ValueError("The plan row has no HiFi read paths")
        missing_reads = [str(path) for path in read_paths if not path.is_file()]
        if missing_reads:
            raise FileNotFoundError(
                "HiFi reads are unavailable: " + ", ".join(missing_reads)
            )

        seed_path = result_dir / "seed.fa"
        seed_sequence = extract_seed(
            plan,
            samtools,
            seed_path,
            logs / "samtools_faidx.log",
        )
        event_start = int_value(plan["query_position_0"]) - int_value(
            plan["seed_start_0"]
        )
        event_end = min(
            len(seed_sequence),
            event_start + int_value(plan["query_event_span_bp"], 1),
        )

        reads_to_seed = result_dir / "reads_to_seed.paf"
        run_command(
            [
                minimap2,
                "-x",
                "map-hifi",
                "-c",
                "--secondary=yes",
                "-N",
                str(args.max_secondary),
                "-t",
                str(args.threads),
                str(seed_path),
                *[str(path) for path in read_paths],
            ],
            reads_to_seed,
            logs / "minimap2_reads_to_seed.log",
        )
        selected_names, spanning_names, source_allele_spanning_names = (
            read_selection_from_paf(
                reads_to_seed,
                event_start,
                event_end,
                args.min_read_mapq,
                args.min_read_aligned_bp,
                args.read_anchor_bp,
                args.max_read_event_discordance_bp,
            )
        )
        selected_reads = result_dir / "selected_reads.fa.gz"
        selected_count = extract_selected_reads(
            read_paths,
            selected_names,
            selected_reads,
        )
        read_failure_reasons = []
        if selected_count < args.min_selected_reads:
            read_failure_reasons.append(
                f"SELECTED_READS<{args.min_selected_reads}"
            )
        if len(source_allele_spanning_names) < args.min_source_allele_reads:
            read_failure_reasons.append(
                f"SOURCE_ALLELE_READS<{args.min_source_allele_reads}"
            )
        if read_failure_reasons:
            write_early_qc_failure(
                result_dir,
                status_path,
                plan,
                args.rescue_id,
                read_failure_reasons,
                seed_length_bp=len(seed_sequence),
                selected_read_count=selected_count,
                spanning_read_count=len(spanning_names),
                source_allele_spanning_read_count=len(
                    source_allele_spanning_names
                ),
            )
            return

        prefix = work / "local"
        hifiasm_command = [
            hifiasm,
            "-o",
            str(prefix),
            "-t",
            str(args.threads),
            "-f0",
            *shlex.split(args.hifiasm_extra),
            str(selected_reads),
        ]
        run_command(
            hifiasm_command,
            logs / "hifiasm.stdout.log",
            logs / "hifiasm.stderr.log",
        )
        gfa_paths = sorted(work.glob("local*.p_ctg.gfa"))
        if not gfa_paths:
            raise RuntimeError("hifiasm produced no primary/haplotype p_ctg GFA")
        local_records = gfa_contigs(gfa_paths)
        if not local_records:
            raise RuntimeError("hifiasm GFA files contain no sequence records")
        local_contigs = result_dir / "local_contigs.fa"
        write_fasta(local_contigs, local_records)

        contigs_to_seed = result_dir / "local_contigs_to_seed.paf"
        run_command(
            [
                minimap2,
                "-x",
                "asm5",
                "-c",
                "--secondary=yes",
                "-N",
                str(args.max_secondary),
                "-t",
                str(args.threads),
                str(seed_path),
                str(local_contigs),
            ],
            contigs_to_seed,
            logs / "minimap2_contigs_to_seed.log",
        )
        selected = select_local_contig(
            parse_paf(contigs_to_seed),
            event_start,
            event_end,
            len(seed_sequence),
            args.contig_anchor_bp,
            args.min_contig_mapq,
            args.min_contig_identity,
            args.discordance_window_bp,
        )
        if selected is None:
            write_early_qc_failure(
                result_dir,
                status_path,
                plan,
                args.rescue_id,
                ["NO_TWO_SIDED_LOCAL_CONTIG"],
                seed_length_bp=len(seed_sequence),
                selected_read_count=selected_count,
                spanning_read_count=len(spanning_names),
                source_allele_spanning_read_count=len(
                    source_allele_spanning_names
                ),
                local_contig_count=len(local_records),
            )
            return

        local_sequences = dict(local_records)
        sequence = local_sequences[selected["query_name"]][
            selected["query_start"] : selected["query_end"]
        ]
        if selected["strand"] == "-":
            sequence = reverse_complement(sequence)
        rescue_name = (
            f"{safe_name(plan['catalog_event_id'])}|{safe_name(plan['sample_id'])}|"
            f"{safe_name(plan['haplotype'])}|{safe_name(plan['source_event_id'])}"
        )
        candidate_rescue = result_dir / "candidate_rescue.fa"
        write_fasta(candidate_rescue, [(rescue_name, sequence)])

        failure_reasons = []
        discordance = selected["event_discordance_bp"]
        if discordance > args.max_event_discordance_bp:
            failure_reasons.append(f"EVENT_DISCORDANCE>{args.max_event_discordance_bp}")
        n_count = sequence.upper().count("N")
        if n_count:
            failure_reasons.append("RESCUE_CONTAINS_N")
        if len(sequence) < args.min_rescue_length_bp:
            failure_reasons.append(f"RESCUE_LENGTH<{args.min_rescue_length_bp}")

        status = "PASS" if not failure_reasons else "FAIL"
        clean_path = ""
        if status == "PASS":
            clean_rescue = result_dir / "clean_rescue.fa"
            write_fasta(clean_rescue, [(rescue_name, sequence)])
            clean_path = str(clean_rescue.resolve())

        qc = {
            "rescue_id": args.rescue_id,
            "status": status,
            "failure_reasons": ";".join(failure_reasons),
            "catalog_event_id": plan["catalog_event_id"],
            "source_event_id": plan["source_event_id"],
            "sample_id": plan["sample_id"],
            "assembly_id": plan["assembly_id"],
            "haplotype": plan["haplotype"],
            "svtype": plan["svtype"],
            "svlen": plan["svlen"],
            "seed_length_bp": len(seed_sequence),
            "selected_read_count": selected_count,
            "spanning_read_count": len(spanning_names),
            "source_allele_spanning_read_count": len(
                source_allele_spanning_names
            ),
            "local_contig_count": len(local_records),
            "selected_local_contig": selected["query_name"],
            "selected_alignment_strand": selected["strand"],
            "selected_alignment_mapq": selected["mapq"],
            "selected_alignment_identity": f"{selected['identity']:.8f}",
            "selected_target_coverage": f"{selected['target_coverage']:.8f}",
            "selected_event_discordance_bp": discordance,
            "rescue_length_bp": len(sequence),
            "rescue_n_count": n_count,
            "rescue_sha256": sequence_sha256(sequence),
            "clean_rescue_fasta": clean_path,
        }
        write_tsv(result_dir / "qc.tsv", [qc], QC_FIELDS)
        status_path.write_text(
            json.dumps(
                {"rescue_id": args.rescue_id, "status": status},
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        print(f"{args.rescue_id}: {status} ({result_dir})")
    except Exception as error:
        status_path.write_text(
            json.dumps(
                {
                    "rescue_id": args.rescue_id,
                    "status": "ERROR",
                    "error": str(error),
                },
                indent=2,
                sort_keys=True,
            )
            + "\n"
        )
        raise


def summarize_rescues(args):
    allow_large_tsv_fields()
    plan = {row["rescue_id"]: row for row in iter_tsv(args.plan)}
    summary = []
    fasta_records = []
    graph_rows = []
    sequence_to_name = {}
    order = 0

    for rescue_id in sorted(plan):
        result_dir = Path(args.results_dir) / rescue_id
        qc_path = result_dir / "qc.tsv"
        if not qc_path.is_file():
            row = {
                "rescue_id": rescue_id,
                "status": "NOT_RUN",
                "failure_reasons": "No qc.tsv was found",
                "catalog_event_id": plan[rescue_id].get("catalog_event_id", ""),
                "combined_fasta_status": "NOT_INCLUDED",
            }
            summary.append(row)
            continue
        qc_rows = read_tsv(qc_path)
        if len(qc_rows) != 1:
            raise ValueError(f"Expected one QC row in {qc_path}")
        row = dict(qc_rows[0])
        row["combined_fasta_status"] = "NOT_INCLUDED"
        row["combined_sequence_name"] = ""
        if row.get("status") == "PASS":
            fasta = Path(row["clean_rescue_fasta"])
            records = parse_fasta(fasta)
            if len(records) != 1:
                raise ValueError(f"Expected one rescue sequence in {fasta}")
            original_name, sequence = next(iter(records.items()))
            digest = sequence_sha256(sequence)
            if digest in sequence_to_name:
                row["combined_fasta_status"] = "DUPLICATE_SEQUENCE"
                row["combined_sequence_name"] = sequence_to_name[digest]
            else:
                order += 1
                name = f"rescue_{order:05d}|{original_name}"
                sequence_to_name[digest] = name
                fasta_records.append((name, sequence))
                row["combined_fasta_status"] = "INCLUDED"
                row["combined_sequence_name"] = name
                graph_rows.append(
                    {
                        "order": order,
                        "role": "validated_local_rescue",
                        "rescue_id": rescue_id,
                        "catalog_event_id": row["catalog_event_id"],
                        "source_event_id": row["source_event_id"],
                        "sample_id": row["sample_id"],
                        "haplotype": row["haplotype"],
                        "svtype": row["svtype"],
                        "svlen": row["svlen"],
                        "sequence_name": name,
                        "sequence_sha256": digest,
                        "source_fasta": str(fasta.resolve()),
                    }
                )
        summary.append(row)

    write_tsv(args.summary, summary, SUMMARY_FIELDS)
    write_fasta(args.fasta, fasta_records)
    graph_fields = [
        "order",
        "role",
        "rescue_id",
        "catalog_event_id",
        "source_event_id",
        "sample_id",
        "haplotype",
        "svtype",
        "svlen",
        "sequence_name",
        "sequence_sha256",
        "source_fasta",
    ]
    write_tsv(args.graph_inputs, graph_rows, graph_fields)
    print(
        f"Included {len(fasta_records)} unique PASS rescues in {args.fasta}; "
        f"summary: {args.summary}"
    )


def add_plan_parser(subparsers):
    parser = subparsers.add_parser(
        "plan",
        help="Select catalog loci and choose one source haplotype/event per locus",
    )
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--graph-candidates", required=True)
    parser.add_argument("--assembly-manifest", required=True)
    parser.add_argument(
        "--read-manifest",
        help=(
            "Optional sample-level TSV with sample_id and hifi_path. If "
            "omitted, read paths are taken from the assembly manifest."
        ),
    )
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--event-ids",
        help=(
            "Optional one-ID-per-line file of explicit catalog PSV IDs. "
            "Explicit IDs bypass catalog confidence/method filters."
        ),
    )
    parser.add_argument(
        "--confidence",
        default="UNCERTAIN",
        help="Default batch filter; ignored with --event-ids",
    )
    parser.add_argument(
        "--validation-status",
        default="PENDING_HIFI",
        help="Default batch filter; ignored with --event-ids",
    )
    parser.add_argument(
        "--require-methods",
        default="graph_residual;assembly_dipcall;assembly_svim_asm",
        help="Semicolon/comma-separated required methods; ignored with --event-ids",
    )
    parser.add_argument(
        "--min-independent-samples",
        type=int,
        default=2,
        help="Default batch filter; ignored with --event-ids",
    )
    parser.add_argument(
        "--max-events",
        type=int,
        default=100,
        help="Maximum default-filtered events; ignored with --event-ids",
    )
    parser.add_argument(
        "--source-qc-tiers",
        default="fragmented_rescue;not_recommended",
        help="Allowed source assembly tiers, or 'all'",
    )
    parser.add_argument("--breakpoint-distance", type=int, default=500)
    parser.add_argument("--length-similarity", type=float, default=0.70)
    parser.add_argument("--flank-bp", type=int, default=50_000)
    parser.add_argument("--min-flank-bp", type=int, default=20_000)
    parser.add_argument(
        "--allow-contig-end",
        action="store_true",
        help="Allow source events with less than --min-flank-bp on one side",
    )
    parser.set_defaults(func=build_plan)


def add_run_parser(subparsers):
    parser = subparsers.add_parser(
        "run",
        help="Recruit HiFi reads and locally reassemble one planned rescue",
    )
    parser.add_argument("--plan", required=True)
    parser.add_argument("--rescue-id", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--max-secondary", type=int, default=20)
    parser.add_argument("--min-read-mapq", type=int, default=10)
    parser.add_argument("--min-read-aligned-bp", type=int, default=2_000)
    parser.add_argument("--read-anchor-bp", type=int, default=1_000)
    parser.add_argument(
        "--max-read-event-discordance-bp",
        type=int,
        default=30,
        help=(
            "Maximum indel bases near the event for a spanning read to count "
            "as support for the source allele"
        ),
    )
    parser.add_argument("--min-selected-reads", type=int, default=8)
    parser.add_argument("--min-source-allele-reads", type=int, default=3)
    parser.add_argument("--contig-anchor-bp", type=int, default=10_000)
    parser.add_argument("--min-contig-mapq", type=int, default=20)
    parser.add_argument("--min-contig-identity", type=float, default=0.99)
    parser.add_argument("--discordance-window-bp", type=int, default=500)
    parser.add_argument("--max-event-discordance-bp", type=int, default=20)
    parser.add_argument("--min-rescue-length-bp", type=int, default=20_000)
    parser.add_argument(
        "--hifiasm-extra",
        default="",
        help="Additional hifiasm arguments, shell-style quoted",
    )
    parser.set_defaults(func=run_rescue)


def add_summarize_parser(subparsers):
    parser = subparsers.add_parser(
        "summarize",
        help="Combine unique PASS rescues and write a graph-input audit table",
    )
    parser.add_argument("--plan", required=True)
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--graph-inputs", required=True)
    parser.set_defaults(func=summarize_rescues)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    add_plan_parser(subparsers)
    add_run_parser(subparsers)
    add_summarize_parser(subparsers)
    return parser.parse_args()


def main():
    args = parse_args()
    if getattr(args, "max_events", 1) < 1:
        raise ValueError("--max-events must be at least 1")
    if getattr(args, "threads", 1) < 1:
        raise ValueError("--threads must be at least 1")
    args.func(args)


if __name__ == "__main__":
    main()
