#!/usr/bin/env python3
"""Rewrite rGFA source tags in resumable byte-range jobs."""

import argparse
import csv
import gzip
import json
import os
import sys
import tempfile
from collections import Counter
from pathlib import Path

from minigraph_cactus_common import atomic_bytes, canonical, file_sha256, read_json


def source_map(path):
    names = {}
    targets = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"gfa_source_name", "pansn_source_name", "seqfile_name"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("seed contig map lacks source-name columns")
        for row in reader:
            source, target = row["gfa_source_name"], row["pansn_source_name"]
            if not source or not target or source.encode() in names or target in targets:
                raise ValueError(f"duplicate/empty rGFA source mapping: {source} -> {target}")
            names[source.encode()] = (target.encode(), row["seqfile_name"])
            targets.add(target)
    if not names:
        raise ValueError("seed contig map is empty")
    return names


def checked_source(path, preflight):
    source = Path(path)
    if source.suffix == ".gz":
        raise ValueError("byte-range seed jobs require an uncompressed source GFA")
    if str(source) != preflight["seed_gfa"]["path"]:
        raise ValueError("seed GFA differs from preflight")
    stat = source.stat()
    if (stat.st_size, stat.st_mtime_ns) != (preflight["seed_gfa"]["size"],
                                         preflight["seed_gfa"]["mtime_ns"]):
        raise ValueError("seed GFA changed since preflight")
    return source, stat


def rewrite_line(line, names, counts, source_events, rank0):
    if line.startswith(b"S\t"):
        fields = line.rstrip(b"\r\n").split(b"\t")
        tagged = [i for i, field in enumerate(fields) if field.startswith(b"SN:Z:")]
        if len(tagged) > 1:
            raise ValueError("rGFA segment has multiple SN tags")
        if tagged:
            source_name = fields[tagged[0]][5:]
            mapped = names.get(source_name)
            if mapped is None:
                # Minigraph can prefix an existing PanSN header with its seqfile label.
                label, separator, suffix = source_name.partition(b"#")
                sample, dot, hap = label.rpartition(b".")
                original = sample + b"#" + suffix
                mapped = names.get(original)
                if (not separator or not dot or hap not in {b"1", b"2"}
                        or not suffix.startswith(hap + b"#") or mapped is None
                        or mapped[1].encode() != label):
                    raise ValueError(
                        f"rGFA source absent from FASTA manifest: {source_name[:180]!r}")
                counts["seqfile_prefixed_segments"] += 1
            target, event = mapped
            fields[tagged[0]] = b"SN:Z:" + target
            source_events[event] += 1
            if b"SR:i:0" in fields:
                rank0[event] += 1
            line = b"\t".join(fields) + b"\n"
        counts["segments"] += 1
        if len(fields) > 2 and fields[2] != b"*":
            counts["segment_bp"] += len(fields[2])
    elif line.startswith(b"L\t"):
        counts["links"] += 1
    elif line.startswith((b"P\t", b"W\t")):
        counts["paths"] += 1
    return line


def normalize_chunk(args):
    preflight = read_json(args.preflight_json)
    source, before = checked_source(args.source_gfa, preflight)
    destination = Path(args.output_gfa)
    if source == destination:
        raise ValueError("seed graph output must differ from the original GFA")
    if not 0 <= args.chunk_index < args.chunks:
        raise ValueError("chunk index is outside the configured range")
    names = source_map(args.contig_map)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    counts, source_events, rank0 = Counter(), Counter(), Counter()
    try:
        with source.open("rb") as input_file, os.fdopen(fd, "wb") as output:
            boundary = before.st_size * args.chunk_index // args.chunks
            end_boundary = before.st_size * (args.chunk_index + 1) // args.chunks
            if boundary:
                input_file.seek(boundary - 1)
                if input_file.read(1) != b"\n":
                    input_file.readline()
            start = input_file.tell()
            while input_file.tell() < end_boundary:
                line = input_file.readline()
                if not line:
                    break
                output.write(rewrite_line(line, names, counts, source_events, rank0))
            end = input_file.tell()
        after = source.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError("source GFA changed during normalization")
        os.replace(temporary, destination)
        report = {"chunk_index": args.chunk_index, "chunk_count": args.chunks,
                  "input_start": start, "input_end": end, "source_size": before.st_size,
                  "output_gfa": str(destination), "output_size": destination.stat().st_size,
                  "counts": dict(counts),
                  "source_segments_by_event": dict(sorted(source_events.items())),
                  "rank0_segments_by_event": dict(sorted(rank0.items()))}
        atomic_bytes(args.report, canonical(report) + b"\n")
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def finalize(args):
    preflight = read_json(args.preflight_json)
    source, source_stat = checked_source(args.source_gfa, preflight)
    destination = Path(args.output_gfa)
    reports = [read_json(path) for path in args.chunk_reports]
    counts, source_events, rank0 = Counter(), Counter(), Counter()
    end = 0
    for index, row in enumerate(reports):
        if (row["chunk_index"] != index or row["chunk_count"] != len(reports)
                or row["input_start"] != end or row["source_size"] != source_stat.st_size):
            raise ValueError(f"seed chunk {index:02d} is missing, out of order or from another source")
        chunk = Path(row["output_gfa"])
        if not chunk.is_file() or chunk.stat().st_size != row["output_size"]:
            raise ValueError(f"seed chunk output changed: {chunk}")
        end = row["input_end"]
        counts.update(row["counts"])
        source_events.update(row["source_segments_by_event"])
        rank0.update(row["rank0_segments_by_event"])
    if not reports or end != source_stat.st_size:
        raise ValueError("seed chunks do not cover the source GFA")
    if rank0["CHM13"] == 0:
        raise ValueError("source GFA has no rank-0 CHM13 segments")
    if counts["segments"] == 0 or counts["links"] == 0:
        raise ValueError("source GFA lacks segments or links")
    if not destination.is_file() or destination.stat().st_size == 0:
        raise ValueError("normalized seed GFA is missing or empty")
    report = {"source_gfa": str(source), "source_sha256": preflight["seed_gfa"]["sha256"],
              "normalized_gfa": str(destination), "normalized_sha256": file_sha256(destination),
              "counts": dict(counts), "source_segments_by_event": dict(sorted(source_events.items())),
              "rank0_segments_by_event": dict(sorted(rank0.items())),
              "normalization": "SN tags only; GFA segment IDs, sequences, links and paths retained"}
    atomic_bytes(args.report, canonical(report) + b"\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    chunk = commands.add_parser("chunk")
    final = commands.add_parser("finalize")
    for command in (chunk, final):
        for key in ("source-gfa", "output-gfa", "preflight-json", "report"):
            command.add_argument(f"--{key}", required=True)
    chunk.add_argument("--contig-map", required=True)
    chunk.add_argument("--chunk-index", type=int, required=True)
    chunk.add_argument("--chunks", type=int, required=True)
    final.add_argument("--chunk-reports", nargs="+", required=True)
    args = parser.parse_args()
    if args.command == "chunk":
        normalize_chunk(args)
    else:
        finalize(args)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, KeyError, UnicodeError, json.JSONDecodeError) as exc:
        sys.exit(f"Minigraph-Cactus seed normalization failed: {exc}")
