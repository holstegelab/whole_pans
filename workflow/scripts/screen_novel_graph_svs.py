#!/usr/bin/env python3
"""Screen complete assemblies for SV-sized differences from a frozen rGFA.

The command is deliberately scheduler-agnostic.  ``tasks`` creates stable,
1-based batches and ``run`` processes one batch.  Snakemake (or a cluster array)
provides the parallelism.
"""

import argparse
import csv
import gzip
import hashlib
import os
import re
import shlex
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


ASSEMBLY_RE = re.compile(
    r"^(?P<sample>.+?)\.hifi\.hifiasm\.bp\.(?P<haplotype>hap[12])\.p_ctg(?:\..*)?$"
)
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
PATH_RE = re.compile(r"([><])([^><]+)")

ALIGNMENT_FIELDS = [
    "assembly_id",
    "graph_member",
    "mapping_pass",
    "query_name",
    "query_length",
    "query_start",
    "query_end",
    "strand",
    "graph_path",
    "path_length",
    "path_start",
    "path_end",
    "matches",
    "block_length",
    "identity",
    "mapq",
    "primary",
    "cigar",
    "alignment_status",
]

EVENT_FIELDS = [
    "event_id",
    "assembly_id",
    "sample_id",
    "haplotype",
    "graph_member",
    "mapping_pass",
    "query_name",
    "query_length",
    "query_position_0",
    "strand",
    "graph_path",
    "path_position_0",
    "graph_segment",
    "segment_offset_0",
    "stable_source",
    "stable_position_0",
    "source_rank",
    "chromosome",
    "svtype",
    "svlen",
    "event_size_bp",
    "cigar_operation",
    "cluster_member_count",
    "left_anchor_bp",
    "right_anchor_bp",
    "alignment_length",
    "identity",
    "mapq",
    "primary",
    "confidence_tier",
    "filter_reasons",
]

COMPLEX_FIELDS = [
    "assembly_id",
    "sample_id",
    "haplotype",
    "graph_member",
    "mapping_pass",
    "query_name",
    "classification",
    "alignment_count",
    "query_breakpoints",
    "graph_paths",
    "strands",
    "mapqs",
    "notes",
]

SUMMARY_FIELDS = [
    "assembly_id",
    "sample_id",
    "haplotype",
    "path",
    "graph_member",
    "rescue_tier",
    "total_contigs",
    "total_bp",
    "primary_alignment_count",
    "secondary_alignment_count",
    "primary_aligned_bp",
    "primary_callable_bp",
    "primary_aligned_fraction",
    "primary_callable_fraction",
    "sensitivity_pass_run",
    "sensitivity_callable_fraction",
    "raw_event_count",
    "high_confidence_event_count",
    "review_event_count",
    "subthreshold_event_count",
    "split_or_complex_contig_count",
    "screen_status",
]

TASK_FIELDS = [
    "task_id",
    "task_index",
    "assembly_index",
    "assembly_id",
    "sample_id",
    "haplotype",
    "path",
    "graph_member",
    "rescue_tier",
    "graph_context",
]

TASK_OUTPUT_FIELDS = ["assembly_id", "kind", "relative_path", "size_bytes"]


def open_text(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def write_tsv(path, rows, fields, compress=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if compress is None:
        compress = str(path).endswith(".gz")
    opener = gzip.open if compress else open
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


def read_tsv(path):
    with open_text(path, "rt") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


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


def parse_sample_haplotype(assembly_id):
    match = ASSEMBLY_RE.match(normalize_assembly_id(assembly_id))
    if match:
        return match.group("sample"), match.group("haplotype")
    return normalize_assembly_id(assembly_id), "unknown"


def manifest_path(row):
    for field in ("cleaned_fasta_path", "path", "cleaned_path"):
        if row.get(field):
            return row[field]
    raise ValueError("Manifest needs one of: cleaned_fasta_path, path, cleaned_path")


def graph_members(path):
    members = set()
    if not path:
        return members
    for row in read_tsv(path):
        if row.get("role", "assembly") != "assembly":
            continue
        value = row.get("assembly_id") or row.get("path")
        if value:
            members.add(normalize_assembly_id(value))
    return members


def parse_tag_fields(fields):
    tags = {}
    for field in fields:
        parts = field.split(":", 2)
        if len(parts) == 3:
            tags[parts[0]] = parts[2]
    return tags


def index_graph(gfa, output):
    rows = []
    seen = set()
    with open_text(gfa, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.startswith("S\t"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"Malformed S record at {gfa}:{line_number}")
            segment = fields[1]
            if segment in seen:
                raise ValueError(f"Duplicate graph segment: {segment}")
            seen.add(segment)
            tags = parse_tag_fields(fields[3:])
            sequence = fields[2]
            if sequence == "*":
                if "LN" not in tags:
                    raise ValueError(f"Segment {segment} has '*' sequence and no LN tag")
                length = int(tags["LN"])
            else:
                length = len(sequence)
                if "LN" in tags and int(tags["LN"]) != length:
                    raise ValueError(f"LN disagrees with sequence length for {segment}")
            rows.append(
                {
                    "segment_id": segment,
                    "length": length,
                    "sn": tags.get("SN", ""),
                    "so": tags.get("SO", ""),
                    "sr": tags.get("SR", ""),
                }
            )
    if not rows:
        raise ValueError(f"No graph segments found in {gfa}")
    write_tsv(output, rows, ["segment_id", "length", "sn", "so", "sr"])


def priority_key(row):
    if str(row.get("graph_context", "")) == "restore_missing_mate_for_graph_sample":
        tier = 0
    elif str(row.get("graph_member", "")).lower() in {"1", "true", "yes"}:
        tier = 1
    else:
        tier = {
            "best_rescue": 2,
            "reasonable_rescue": 3,
            "fragmented_rescue": 4,
            "not_recommended": 5,
        }.get(row.get("rescue_tier", ""), 6)
    return tier, row.get("sample_id", row.get("sample", "")), row["assembly_id"]


def build_tasks(manifest, graph_assemblies, output, batch_size=1):
    if batch_size < 1:
        raise ValueError("--batch-size must be at least 1")
    members = graph_members(graph_assemblies)
    rows = []
    seen = set()
    for source in read_tsv(manifest):
        assembly_id = normalize_assembly_id(source["assembly_id"])
        if assembly_id in seen:
            raise ValueError(f"Duplicate assembly ID in manifest: {assembly_id}")
        seen.add(assembly_id)
        sample, haplotype = parse_sample_haplotype(assembly_id)
        row = dict(source)
        row.update(
            {
                "assembly_id": assembly_id,
                "path": manifest_path(source),
                "sample_id": source.get("sample_id", source.get("sample", sample)),
                "haplotype": source.get("haplotype", haplotype),
                "graph_member": "true" if assembly_id in members or str(source.get("graph_member", "")).lower() in {"1", "true", "yes"} else "false",
                "rescue_tier": source.get("rescue_tier", "graph_member" if assembly_id in members else "unclassified"),
            }
        )
        rows.append(row)
    rows.sort(key=priority_key)
    task_rows = []
    for index, row in enumerate(rows, start=1):
        task_index = (index - 1) // batch_size + 1
        task_rows.append(
            {
                "task_id": f"{task_index:04d}",
                "task_index": task_index,
                "assembly_index": index,
                "assembly_id": row["assembly_id"],
                "sample_id": row["sample_id"],
                "haplotype": row["haplotype"],
                "path": row["path"],
                "graph_member": row["graph_member"],
                "rescue_tier": row["rescue_tier"],
                "graph_context": row.get("graph_context", ""),
            }
        )
    write_tsv(output, task_rows, TASK_FIELDS, compress=False)


def load_segment_index(path):
    result = {}
    for row in read_tsv(path):
        result[row["segment_id"]] = {
            "length": int(row["length"]),
            "sn": row.get("sn", ""),
            "so": int(row["so"]) if row.get("so", "") not in {"", "NA"} else None,
            "sr": int(row["sr"]) if row.get("sr", "") not in {"", "NA"} else None,
        }
    return result


def parse_graph_path(path):
    if path.startswith(">") or path.startswith("<"):
        parsed = PATH_RE.findall(path)
        if "".join(direction + name for direction, name in parsed) != path:
            return []
        return parsed
    return [(">", path)]


def normalize_chromosome(stable_source):
    value = stable_source
    for prefix in ("CHM13.", "GRCh38."):
        if value.startswith(prefix):
            value = value[len(prefix) :]
    return value


def locate_graph_position(path, path_position, segments):
    traversal = parse_graph_path(path)
    if not traversal:
        return {"segment": "", "offset": "", "sn": "", "so": "", "sr": "", "stable_position": ""}
    remaining = max(0, int(path_position))
    for direction, name in traversal:
        segment = segments.get(name)
        if not segment:
            return {"segment": "", "offset": "", "sn": "", "so": "", "sr": "", "stable_position": ""}
        length = segment["length"]
        if remaining < length or (remaining == length and name == traversal[-1][1]):
            local = min(remaining, max(0, length - 1))
            offset = local if direction == ">" else max(0, length - local - 1)
            stable = segment["so"] + offset if segment["so"] is not None else ""
            return {
                "segment": name,
                "offset": offset,
                "sn": segment["sn"],
                "so": segment["so"] if segment["so"] is not None else "",
                "sr": segment["sr"] if segment["sr"] is not None else "",
                "stable_position": stable,
            }
        remaining -= length
    return {"segment": "", "offset": "", "sn": "", "so": "", "sr": "", "stable_position": ""}


def parse_cigar(cigar):
    tokens = [(int(length), operation) for length, operation in CIGAR_RE.findall(cigar)]
    if not tokens or "".join(f"{length}{op}" for length, op in tokens) != cigar:
        raise ValueError(f"Malformed CIGAR: {cigar}")
    return tokens


def aligned_anchors(tokens, token_indices):
    first = min(token_indices)
    last = max(token_indices)
    aligned_ops = {"M", "=", "X"}
    left = sum(length for length, op in tokens[:first] if op in aligned_ops)
    right = sum(length for length, op in tokens[last + 1 :] if op in aligned_ops)
    return left, right


def extract_cigar_events(cigar, query_start=0, target_start=0, raw_min_size=30, cluster_gap=50):
    """Return atomic indels and qualifying adjacent-indel clusters.

    Coordinates are zero-based offsets in the query and GAF path traversal.
    A cluster is emitted in addition to its atomic sensitivity records so no
    30--49 bp evidence is lost before local normalization.
    """
    tokens = parse_cigar(cigar)
    query_position = int(query_start)
    target_position = int(target_start)
    all_indels = []
    for token_index, (length, operation) in enumerate(tokens):
        if operation in {"I", "D"}:
            all_indels.append(
                {
                    "kind": "atomic",
                    "operation": operation,
                    "length": length,
                    "query_position": query_position,
                    "target_position": target_position,
                    "token_indices": [token_index],
                    "members": 1,
                }
            )
        if operation in {"M", "=", "X", "I", "S"}:
            query_position += length
        if operation in {"M", "=", "X", "D", "N"}:
            target_position += length

    atoms = [atom for atom in all_indels if atom["length"] >= raw_min_size]
    clusters = []
    current = []
    for atom in all_indels:
        if not current:
            current = [atom]
            continue
        previous = current[-1]
        token_gap = tokens[previous["token_indices"][0] + 1 : atom["token_indices"][0]]
        query_gap = sum(length for length, op in token_gap if op in {"M", "=", "X", "I", "S"})
        target_gap = sum(length for length, op in token_gap if op in {"M", "=", "X", "D", "N"})
        if max(query_gap, target_gap) <= cluster_gap:
            current.append(atom)
        else:
            if len(current) > 1:
                clusters.append(current)
            current = [atom]
    if len(current) > 1:
        clusters.append(current)

    events = list(atoms)
    for cluster in clusters:
        if sum(atom["length"] for atom in cluster) < 50:
            continue
        insertion_bp = sum(atom["length"] for atom in cluster if atom["operation"] == "I")
        deletion_bp = sum(atom["length"] for atom in cluster if atom["operation"] == "D")
        operation = "I" if deletion_bp == 0 else "D" if insertion_bp == 0 else "ID"
        events.append(
            {
                "kind": "cluster",
                "operation": operation,
                "length": insertion_bp + deletion_bp,
                "svlen": insertion_bp - deletion_bp,
                "query_position": cluster[0]["query_position"],
                "target_position": cluster[0]["target_position"],
                "token_indices": [atom["token_indices"][0] for atom in cluster],
                "members": len(cluster),
            }
        )

    for event in events:
        event["left_anchor"], event["right_anchor"] = aligned_anchors(
            tokens, event["token_indices"]
        )
        if "svlen" not in event:
            event["svlen"] = event["length"] if event["operation"] == "I" else -event["length"]
    return events


def parse_gaf(path, mapping_pass):
    records = []
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise ValueError(f"Malformed GAF line at {path}:{line_number}")
            tags = parse_tag_fields(fields[12:])
            block_length = int(fields[10])
            record = {
                "mapping_pass": mapping_pass,
                "query_name": fields[0],
                "query_length": int(fields[1]),
                "query_start": int(fields[2]),
                "query_end": int(fields[3]),
                "strand": fields[4],
                "graph_path": fields[5],
                "path_length": int(fields[6]),
                "path_start": int(fields[7]),
                "path_end": int(fields[8]),
                "matches": int(fields[9]),
                "block_length": block_length,
                "identity": int(fields[9]) / block_length if block_length else 0.0,
                "mapq": int(fields[11]),
                "primary": tags.get("tp", "P") == "P",
                "cigar": tags.get("cg", ""),
            }
            record["alignment_status"] = "OK" if record["cigar"] else "NO_CIGAR"
            records.append(record)
    return records


def merge_intervals(intervals):
    if not intervals:
        return []
    merged = []
    start, end = sorted(intervals)[0]
    for next_start, next_end in sorted(intervals)[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            merged.append((start, end))
            start, end = next_start, next_end
    merged.append((start, end))
    return merged


def interval_bp(intervals):
    return sum(end - start for start, end in merge_intervals(intervals))


def fasta_lengths(path):
    lengths = {}
    current = None
    length = 0
    with open_text(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if current is not None:
                    lengths[current] = length
                current = line[1:].strip().split()[0]
                if not current:
                    raise ValueError(f"Empty FASTA identifier at {path}:{line_number}")
                if current in lengths:
                    raise ValueError(f"Duplicate FASTA identifier in {path}: {current}")
                length = 0
            elif current is not None:
                length += len(line.strip())
            elif line.strip():
                raise ValueError(f"Sequence before first FASTA header at {path}:{line_number}")
    if current is not None:
        lengths[current] = length
    if not lengths:
        raise ValueError(f"No sequences found in FASTA: {path}")
    return lengths


def callable_metrics(records, total_bp, min_alignment, min_identity, min_mapq):
    aligned = defaultdict(list)
    callable_intervals = defaultdict(list)
    for record in records:
        if not record["primary"]:
            continue
        interval = (record["query_start"], record["query_end"])
        aligned[record["query_name"]].append(interval)
        if (
            record["block_length"] >= min_alignment
            and record["identity"] >= min_identity
            and record["mapq"] >= min_mapq
        ):
            callable_intervals[record["query_name"]].append(interval)
    aligned_bp = sum(interval_bp(values) for values in aligned.values())
    callable_bp = sum(interval_bp(values) for values in callable_intervals.values())
    denominator = max(1, total_bp)
    return aligned_bp, callable_bp, aligned_bp / denominator, callable_bp / denominator


def classify_split_alignments(records, task):
    by_query = defaultdict(list)
    for record in records:
        if record["primary"]:
            by_query[(record["mapping_pass"], record["query_name"])].append(record)
    rows = []
    for (mapping_pass, query_name), alignments in sorted(by_query.items()):
        if len(alignments) < 2:
            continue
        ordered = sorted(alignments, key=lambda row: (row["query_start"], row["query_end"]))
        strands = {row["strand"] for row in ordered}
        path_heads = []
        for row in ordered:
            traversal = parse_graph_path(row["graph_path"])
            path_heads.append(traversal[0][1] if traversal else row["graph_path"])
        classification = "SPLIT_ALIGNMENT"
        notes = []
        orientation_change = len(strands) > 1
        distant_join = len(set(path_heads)) > 1
        if orientation_change:
            classification = "POTENTIAL_INVERSION_OR_COMPLEX"
            notes.append("primary alignments have inconsistent orientation")
        if distant_join:
            classification = "POTENTIAL_DISTANT_JOIN_OR_TRANSLOCATION"
            notes.append("primary alignments begin on different graph segments")
        if orientation_change and distant_join:
            classification = "POTENTIAL_INVERSION_OR_TRANSLOCATION"
        rows.append(
            {
                "assembly_id": task["assembly_id"],
                "sample_id": task["sample_id"],
                "haplotype": task["haplotype"],
                "graph_member": task["graph_member"],
                "mapping_pass": mapping_pass,
                "query_name": query_name,
                "classification": classification,
                "alignment_count": len(ordered),
                "query_breakpoints": ";".join(
                    f"{row['query_start']}-{row['query_end']}" for row in ordered
                ),
                "graph_paths": ";".join(row["graph_path"] for row in ordered),
                "strands": ";".join(row["strand"] for row in ordered),
                "mapqs": ";".join(str(row["mapq"]) for row in ordered),
                "notes": "; ".join(notes),
            }
        )
    return rows


def events_from_alignments(records, task, segments, thresholds):
    rows = []
    for alignment_index, record in enumerate(records, start=1):
        if not record["cigar"]:
            continue
        events = extract_cigar_events(
            record["cigar"],
            query_start=record["query_start"],
            target_start=record["path_start"],
            raw_min_size=thresholds["raw_min_size"],
            cluster_gap=thresholds["cluster_gap"],
        )
        for event_index, event in enumerate(events, start=1):
            operation = event["operation"]
            svtype = "INS" if operation == "I" else "DEL" if operation == "D" else "COMPLEX_INDEL"
            reasons = []
            if event["length"] < thresholds["min_sv_size"]:
                reasons.append("BELOW_50_BP_PENDING_NORMALIZATION")
            if not record["primary"]:
                reasons.append("SECONDARY_ALIGNMENT")
            if record["block_length"] < thresholds["min_alignment"]:
                reasons.append("SHORT_ALIGNMENT")
            if record["identity"] < thresholds["min_identity"]:
                reasons.append("LOW_IDENTITY")
            if record["mapq"] < thresholds["min_mapq"]:
                reasons.append("LOW_MAPQ")
            if event["left_anchor"] < thresholds["min_anchor"]:
                reasons.append("WEAK_LEFT_ANCHOR")
            if event["right_anchor"] < thresholds["min_anchor"]:
                reasons.append("WEAK_RIGHT_ANCHOR")
            if event["left_anchor"] == 0 or event["right_anchor"] == 0:
                reasons.append("TERMINAL_EVENT_WITHOUT_SECOND_BREAKPOINT")
            if event["length"] < thresholds["min_sv_size"]:
                tier = "SENSITIVITY_SUBTHRESHOLD"
            elif reasons:
                tier = "REVIEW"
            else:
                tier = "HIGH_CONFIDENCE"

            location = locate_graph_position(
                record["graph_path"], event["target_position"], segments
            )
            event_key = (
                f"{task['assembly_id']}|{record['mapping_pass']}|{record['query_name']}|"
                f"{alignment_index}|{event_index}|{event['query_position']}|{operation}|{event['length']}"
            )
            event_id = "GSV_" + hashlib.sha1(event_key.encode()).hexdigest()[:16]
            rows.append(
                {
                    "event_id": event_id,
                    "assembly_id": task["assembly_id"],
                    "sample_id": task["sample_id"],
                    "haplotype": task["haplotype"],
                    "graph_member": task["graph_member"],
                    "mapping_pass": record["mapping_pass"],
                    "query_name": record["query_name"],
                    "query_length": record["query_length"],
                    "query_position_0": event["query_position"],
                    "strand": record["strand"],
                    "graph_path": record["graph_path"],
                    "path_position_0": event["target_position"],
                    "graph_segment": location["segment"],
                    "segment_offset_0": location["offset"],
                    "stable_source": location["sn"],
                    "stable_position_0": location["stable_position"],
                    "source_rank": location["sr"],
                    "chromosome": normalize_chromosome(location["sn"]),
                    "svtype": svtype,
                    "svlen": event["svlen"],
                    "event_size_bp": event["length"],
                    "cigar_operation": operation,
                    "cluster_member_count": event["members"],
                    "left_anchor_bp": event["left_anchor"],
                    "right_anchor_bp": event["right_anchor"],
                    "alignment_length": record["block_length"],
                    "identity": f"{record['identity']:.6f}",
                    "mapq": record["mapq"],
                    "primary": str(record["primary"]).lower(),
                    "confidence_tier": tier,
                    "filter_reasons": ";".join(dict.fromkeys(reasons)),
                }
            )
    return rows


def run_minigraph(graph, assembly, output, log, threads, min_chain_score, secondary, extra):
    command = [
        "minigraph",
        "-cxasm",
        f"-l{min_chain_score}",
        "-N",
        str(secondary),
        "-t",
        str(threads),
    ] + shlex.split(extra) + [str(graph), str(assembly)]
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(log).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(output, "wb") as stdout, open(log, "ab") as stderr:
        stderr.write(("COMMAND: " + shlex.join(command) + "\n").encode())
        subprocess.run(command, stdout=stdout, stderr=stderr, check=True)


def analyze_one(task, args, segments):
    assembly_id = task["assembly_id"]
    output_dir = Path(args.output_dir)
    primary_gaf = output_dir / "alignments" / f"{assembly_id}.gaf.gz"
    sensitivity_gaf = output_dir / "alignments" / f"{assembly_id}.sensitivity.gaf.gz"
    log = output_dir / "logs" / f"{assembly_id}.minigraph.log"
    if log.exists():
        log.unlink()
    run_minigraph(
        args.graph,
        task["path"],
        primary_gaf,
        log,
        args.threads,
        args.primary_min_chain_score,
        args.primary_secondary,
        args.minigraph_extra,
    )
    primary = parse_gaf(primary_gaf, "primary")
    lengths = fasta_lengths(task["path"])
    total_bp = sum(lengths.values())
    primary_metrics = callable_metrics(
        primary,
        total_bp,
        args.min_alignment,
        args.min_identity,
        args.min_mapq,
    )

    sensitivity = []
    if primary_metrics[3] < args.rescue_callable_fraction:
        run_minigraph(
            args.graph,
            task["path"],
            sensitivity_gaf,
            log,
            args.threads,
            args.sensitivity_min_chain_score,
            args.sensitivity_secondary,
            args.sensitivity_minigraph_extra,
        )
        sensitivity = parse_gaf(sensitivity_gaf, "sensitivity")
    elif sensitivity_gaf.exists():
        sensitivity_gaf.unlink()

    records = primary + sensitivity
    for record in records:
        record["assembly_id"] = assembly_id
        record["graph_member"] = task["graph_member"]
    thresholds = {
        "raw_min_size": args.raw_min_size,
        "min_sv_size": args.min_sv_size,
        "cluster_gap": args.adjacent_indel_gap,
        "min_alignment": args.min_alignment,
        "min_identity": args.min_identity,
        "min_mapq": args.min_mapq,
        "min_anchor": args.min_anchor,
    }
    events = events_from_alignments(records, task, segments, thresholds)
    complex_rows = classify_split_alignments(records, task)
    sensitivity_metrics = callable_metrics(
        sensitivity,
        total_bp,
        args.min_alignment,
        args.min_identity,
        args.min_mapq,
    ) if sensitivity else (0, 0, 0.0, 0.0)

    primary_count = sum(record["primary"] for record in primary)
    secondary_count = len(primary) - primary_count
    high_count = sum(row["confidence_tier"] == "HIGH_CONFIDENCE" for row in events)
    review_count = sum(row["confidence_tier"] == "REVIEW" for row in events)
    subthreshold_count = sum(row["confidence_tier"] == "SENSITIVITY_SUBTHRESHOLD" for row in events)
    if not records or max(primary_metrics[3], sensitivity_metrics[3]) == 0:
        status = "UNASSESSABLE_NO_CALLABLE_ALIGNMENT"
    elif complex_rows:
        status = "REVIEW_SPLIT_OR_COMPLEX_ALIGNMENT"
    elif high_count:
        status = "CANDIDATE_GRAPH_MISSING_SV"
    elif review_count or subthreshold_count:
        status = "REVIEW_RESIDUAL_EVENT"
    else:
        status = "SCREENED_NO_HIGH_CONFIDENCE_RESIDUAL"

    summary = {
        "assembly_id": assembly_id,
        "sample_id": task["sample_id"],
        "haplotype": task["haplotype"],
        "path": task["path"],
        "graph_member": task["graph_member"],
        "rescue_tier": task.get("rescue_tier", ""),
        "total_contigs": len(lengths),
        "total_bp": total_bp,
        "primary_alignment_count": primary_count,
        "secondary_alignment_count": secondary_count,
        "primary_aligned_bp": primary_metrics[0],
        "primary_callable_bp": primary_metrics[1],
        "primary_aligned_fraction": f"{primary_metrics[2]:.8f}",
        "primary_callable_fraction": f"{primary_metrics[3]:.8f}",
        "sensitivity_pass_run": str(bool(sensitivity)).lower(),
        "sensitivity_callable_fraction": f"{sensitivity_metrics[3]:.8f}",
        "raw_event_count": len(events),
        "high_confidence_event_count": high_count,
        "review_event_count": review_count,
        "subthreshold_event_count": subthreshold_count,
        "split_or_complex_contig_count": len({row["query_name"] for row in complex_rows}),
        "screen_status": status,
    }
    write_tsv(output_dir / "alignment_tables" / f"{assembly_id}.alignments.tsv.gz", records, ALIGNMENT_FIELDS)
    write_tsv(output_dir / "candidates" / f"{assembly_id}.residual_svs.tsv.gz", events, EVENT_FIELDS)
    write_tsv(output_dir / "complex" / f"{assembly_id}.complex.tsv.gz", complex_rows, COMPLEX_FIELDS)
    write_tsv(output_dir / "summaries" / f"{assembly_id}.summary.tsv", [summary], SUMMARY_FIELDS, compress=False)


def run_task(args):
    if args.tasks:
        tasks = read_tsv(args.tasks)
    else:
        temporary = Path(args.output_dir) / "tasks.tsv"
        build_tasks(args.manifest, args.graph_assemblies, temporary, args.batch_size)
        tasks = read_tsv(temporary)
    task_id = args.task_id or (f"{args.task_index:04d}" if args.task_index is not None else None)
    selected = [row for row in tasks if row.get("task_id", f"{int(row['task_index']):04d}") == task_id]
    if not selected:
        available = sorted(
            {row.get("task_id", f"{int(row['task_index']):04d}") for row in tasks}
        )
        raise ValueError(
            f"Task {task_id} does not exist; available range is "
            f"{available[:1]}..{available[-1:]}"
        )
    segments = load_segment_index(args.segment_index)
    final_dir = Path(args.output_dir)
    final_dir.parent.mkdir(parents=True, exist_ok=True)
    if final_dir.exists():
        raise FileExistsError(
            f"Task output directory already exists: {final_dir}. "
            "Remove this one task output before a manual rerun; Snakemake manages it automatically."
        )
    temporary_dir = Path(
        tempfile.mkdtemp(prefix=f".{final_dir.name}.", dir=final_dir.parent)
    )
    original_output_dir = args.output_dir
    try:
        args.output_dir = str(temporary_dir)
        for task in selected:
            analyze_one(task, args, segments)
        output_rows = []
        for task in selected:
            assembly_id = task["assembly_id"]
            expected = {
                "primary_gaf": Path("alignments") / f"{assembly_id}.gaf.gz",
                "alignment_table": Path("alignment_tables") / f"{assembly_id}.alignments.tsv.gz",
                "candidates": Path("candidates") / f"{assembly_id}.residual_svs.tsv.gz",
                "complex": Path("complex") / f"{assembly_id}.complex.tsv.gz",
                "summary": Path("summaries") / f"{assembly_id}.summary.tsv",
                "minigraph_log": Path("logs") / f"{assembly_id}.minigraph.log",
            }
            for kind, relative in expected.items():
                path = temporary_dir / relative
                if not path.is_file() or path.stat().st_size == 0:
                    raise RuntimeError(f"Missing or empty {kind} output for {assembly_id}: {path}")
                output_rows.append(
                    {
                        "assembly_id": assembly_id,
                        "kind": kind,
                        "relative_path": str(relative),
                        "size_bytes": path.stat().st_size,
                    }
                )
        write_tsv(
            temporary_dir / "task_outputs.tsv",
            output_rows,
            TASK_OUTPUT_FIELDS,
            compress=False,
        )
        (temporary_dir / ".complete").write_text(
            "\n".join(row["assembly_id"] for row in selected) + "\n"
        )
        temporary_dir.replace(final_dir)
    except Exception:
        # Keep failed task data for diagnosis but outside the declared output.
        failed_dir = temporary_dir.with_name(temporary_dir.name + ".failed")
        if temporary_dir.exists():
            temporary_dir.replace(failed_dir)
        raise
    finally:
        args.output_dir = original_output_dir
    if args.completion_marker:
        marker = Path(args.completion_marker)
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.write_text("\n".join(row["assembly_id"] for row in selected) + "\n")


def concatenate_tables(paths, output, fields):
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if str(output).endswith(".gz") else open
    with opener(output, "wt", newline="") as output_handle:
        writer = csv.DictWriter(
            output_handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for path in paths:
            with open_text(path, "rt") as input_handle:
                reader = csv.DictReader(input_handle, delimiter="\t")
                if reader.fieldnames != fields:
                    raise ValueError(
                        f"Unexpected columns in {path}: {reader.fieldnames}; expected {fields}"
                    )
                writer.writerows(reader)


def summarize(args):
    output_dir = Path(args.output_dir)
    manifest_rows = read_tsv(args.manifest)
    assembly_ids = [normalize_assembly_id(row["assembly_id"]) for row in manifest_rows]
    task_roots = {}
    if args.task_dirs:
        if not args.tasks:
            raise ValueError("--tasks is required with --task-dir")
        for path in map(Path, args.task_dirs):
            task_id = path.name
            if task_id in task_roots:
                raise ValueError(f"Duplicate task directory for {task_id}: {path}")
            task_roots[task_id] = path
        task_rows = read_tsv(args.tasks)
        expected_task_ids = {
            row.get("task_id", f"{int(row['task_index']):04d}") for row in task_rows
        }
        if set(task_roots) != expected_task_ids:
            missing_tasks = sorted(expected_task_ids - set(task_roots))
            extra_tasks = sorted(set(task_roots) - expected_task_ids)
            raise ValueError(
                f"Task directory set does not match tasks.tsv; "
                f"missing={missing_tasks[:20]}, extra={extra_tasks[:20]}"
            )
        assembly_roots = {
            normalize_assembly_id(row["assembly_id"]): task_roots[
                row.get("task_id", f"{int(row['task_index']):04d}")
            ]
            for row in task_rows
        }
        if set(assembly_roots) != set(assembly_ids):
            raise ValueError("Assembly IDs in tasks.tsv do not match the discovery manifest")
    else:
        # Backward-compatible manual mode for pre-checkpoint output layouts.
        assembly_roots = {assembly_id: output_dir for assembly_id in assembly_ids}
    missing = []
    summary_paths = []
    candidate_paths = []
    complex_paths = []
    alignment_paths = []
    for assembly_id in assembly_ids:
        task_root = assembly_roots[assembly_id]
        expected = {
            "summary": task_root / "summaries" / f"{assembly_id}.summary.tsv",
            "candidate": task_root / "candidates" / f"{assembly_id}.residual_svs.tsv.gz",
            "complex": task_root / "complex" / f"{assembly_id}.complex.tsv.gz",
            "alignment": task_root / "alignment_tables" / f"{assembly_id}.alignments.tsv.gz",
        }
        for kind, path in expected.items():
            if not path.is_file() or path.stat().st_size == 0:
                missing.append(f"{assembly_id}:{kind}:{path}")
        summary_paths.append(expected["summary"])
        candidate_paths.append(expected["candidate"])
        complex_paths.append(expected["complex"])
        alignment_paths.append(expected["alignment"])
    if missing:
        preview = "\n".join(missing[:20])
        raise FileNotFoundError(f"Missing per-assembly screen outputs ({len(missing)}):\n{preview}")
    summary_dir = output_dir / "summary"
    concatenate_tables(summary_paths, summary_dir / "all_assembly_novel_sv_summary.tsv", SUMMARY_FIELDS)
    concatenate_tables(candidate_paths, summary_dir / "all_residual_sv_candidates.tsv.gz", EVENT_FIELDS)
    concatenate_tables(complex_paths, summary_dir / "all_complex_alignments.tsv.gz", COMPLEX_FIELDS)
    concatenate_tables(alignment_paths, summary_dir / "all_contig_alignments.tsv.gz", ALIGNMENT_FIELDS)


def add_screen_thresholds(parser):
    parser.add_argument("--min-sv-size", type=int, default=50)
    parser.add_argument("--raw-min-size", type=int, default=30)
    parser.add_argument("--adjacent-indel-gap", type=int, default=50)
    parser.add_argument("--min-alignment", type=int, default=5000)
    parser.add_argument("--min-anchor", type=int, default=2000)
    parser.add_argument("--min-mapq", type=int, default=5)
    parser.add_argument("--min-identity", type=float, default=0.90)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    index_parser = subparsers.add_parser("index-graph", help="Index rGFA segment coordinates")
    index_parser.add_argument("--gfa", required=True)
    index_parser.add_argument("--output", required=True)

    tasks_parser = subparsers.add_parser("tasks", help="Create deterministic assembly batches")
    tasks_parser.add_argument("--manifest", required=True)
    tasks_parser.add_argument("--graph-assemblies", required=True)
    tasks_parser.add_argument("--output", required=True)
    tasks_parser.add_argument("--batch-size", type=int, default=1)

    run_parser = subparsers.add_parser("run", help="Map and analyze one 1-based task batch")
    run_parser.add_argument("--graph", required=True)
    run_parser.add_argument("--segment-index", required=True)
    run_parser.add_argument("--manifest", required=True)
    run_parser.add_argument("--graph-assemblies", required=True)
    run_parser.add_argument("--tasks")
    run_parser.add_argument("--output-dir", required=True)
    task_group = run_parser.add_mutually_exclusive_group(required=True)
    task_group.add_argument("--task-id")
    task_group.add_argument("--task-index", type=int)
    run_parser.add_argument("--batch-size", type=int, default=1)
    run_parser.add_argument("--threads", type=int, default=16)
    run_parser.add_argument("--primary-min-chain-score", type=int, default=5000)
    run_parser.add_argument("--primary-secondary", type=int, default=5)
    run_parser.add_argument("--minigraph-extra", default="")
    run_parser.add_argument("--rescue-callable-fraction", type=float, default=0.85)
    run_parser.add_argument("--sensitivity-min-chain-score", type=int, default=1000)
    run_parser.add_argument("--sensitivity-secondary", type=int, default=20)
    run_parser.add_argument("--sensitivity-minigraph-extra", default="")
    run_parser.add_argument("--completion-marker")
    add_screen_thresholds(run_parser)

    summary_parser = subparsers.add_parser("summarize", help="Combine all assembly outputs")
    summary_parser.add_argument("--manifest", required=True)
    summary_parser.add_argument("--tasks")
    summary_parser.add_argument("--task-dir", dest="task_dirs", action="append")
    summary_parser.add_argument("--output-dir", required=True)
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    if args.command == "index-graph":
        index_graph(args.gfa, args.output)
    elif args.command == "tasks":
        build_tasks(args.manifest, args.graph_assemblies, args.output, args.batch_size)
    elif args.command == "run":
        run_task(args)
    elif args.command == "summarize":
        summarize(args)
    else:
        raise AssertionError(args.command)


if __name__ == "__main__":
    main()
