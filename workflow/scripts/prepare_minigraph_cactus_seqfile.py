#!/usr/bin/env python3
"""Validate selected FASTAs on the server and write Cactus input metadata."""

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import re
import sys
import tempfile
from pathlib import Path

from minigraph_cactus_common import atomic_bytes, canonical, read_json, sha256


ASSEMBLY = re.compile(r"^(.+)\.hifi\.hifiasm\.bp\.hap([12])\.p_ctg\.clean\.fa(?:\.gz)?$")
SUBRANGE = re.compile(r":[0-9]+-[0-9]+$")
BASES = b"ACGTNRYKMSWBDHVacgtnrykmswbdhv"


def required_path(raw):
    path = Path(raw)
    if not path.is_absolute() or any(char.isspace() for char in raw):
        raise ValueError(f"FASTA path must be absolute and contain no whitespace: {raw}")
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing or empty FASTA: {path}")
    return path


def effective_name(header):
    name = header.rsplit("#", 1)[-1].split(" ", 1)[0].split("\t", 1)[0]
    match = SUBRANGE.search(name)
    if match:
        return name[:match.start()].replace(":", "_") + match.group()
    return name.replace(":", "_")


def scan_fasta(path, label, contig_writer):
    before = path.stat()
    digest = hashlib.sha256()
    seen = set()
    contigs = []
    total_bases = 0
    original = None
    effective = None
    sequence_length = 0
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rb") as handle:
        for line_no, line in enumerate(handle, 1):
            digest.update(line)
            if line.startswith(b">"):
                if original is not None:
                    if sequence_length == 0:
                        raise ValueError(f"empty contig before line {line_no}: {path}")
                    contig_writer.writerow([label, original, effective, sequence_length])
                    total_bases += sequence_length
                try:
                    original = line[1:].rstrip(b"\r\n").decode("utf-8")
                except UnicodeDecodeError as exc:
                    raise ValueError(f"non-UTF-8 FASTA header at {path}:{line_no}") from exc
                if not original or any(c in original for c in "\r\n"):
                    raise ValueError(f"empty or malformed FASTA header at {path}:{line_no}")
                effective = effective_name(original)
                if not effective or effective in seen:
                    raise ValueError(f"empty or duplicate Cactus contig ID {effective!r}: {path}:{line_no}")
                seen.add(effective)
                contigs.append(effective)
                sequence_length = 0
            else:
                bases = line.rstrip(b"\r\n")
                if not seen or not bases or bases.translate(None, BASES):
                    raise ValueError(f"malformed FASTA sequence at {path}:{line_no}")
                sequence_length += len(bases)
    if not seen or sequence_length == 0:
        raise ValueError(f"FASTA has no nonempty contigs: {path}")
    contig_writer.writerow([label, original, effective, sequence_length])
    total_bases += sequence_length
    after = path.stat()
    if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
        raise ValueError(f"FASTA changed while scanning: {path}")
    return {"size": after.st_size, "mtime_ns": after.st_mtime_ns,
            "sha256": digest.hexdigest(), "contig_count": len(seen),
            "total_bases": total_bases}, contigs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("included-list", "chm13", "grch38", "seqfile", "input-manifest",
                 "contig-map", "preflight-json"):
        parser.add_argument(f"--{name}", required=True)
    args = parser.parse_args()
    included = Path(args.included_list)
    if not included.is_absolute() or not included.is_file():
        raise ValueError("included list must be an existing absolute server path")
    list_before = included.stat()
    list_bytes = included.read_bytes()
    list_after = included.stat()
    if (list_before.st_size, list_before.st_mtime_ns) != (list_after.st_size, list_after.st_mtime_ns):
        raise ValueError("included list changed while reading")
    selected = []
    labels = {"CHM13", "GRCh38"}
    paths = set()
    for line_no, raw in enumerate(list_bytes.decode("utf-8").splitlines(), 1):
        if not raw.strip() or raw.lstrip().startswith("#"):
            continue
        if raw != raw.strip():
            raise ValueError(f"whitespace around FASTA path at inclusion-list line {line_no}")
        path = required_path(raw)
        match = ASSEMBLY.fullmatch(path.name)
        if not match or "." in match.group(1) or any(c.isspace() for c in match.group(1)):
            raise ValueError(f"unrecognized assembly name at inclusion-list line {line_no}: {path.name}")
        name = f"{match.group(1)}.{match.group(2)}"
        if name in labels or raw in paths:
            raise ValueError(f"duplicate sample label or FASTA path at line {line_no}: {raw}")
        labels.add(name)
        paths.add(raw)
        selected.append((name, match.group(1), match.group(2), path, line_no))
    if not selected:
        raise ValueError("inclusion list has no assembly FASTAs")
    references = [("CHM13", "CHM13", "", required_path(args.chm13), 0),
                  ("GRCh38", "GRCh38", "", required_path(args.grch38), 0)]
    for _, _, _, path, _ in references:
        if str(path) in paths:
            raise ValueError(f"reference FASTA also appears in inclusion list: {path}")
        paths.add(str(path))
    rows = references + selected
    seqfile = "".join(f"{name} {path}\n" for name, _, _, path, _ in rows).encode()
    inputs = []
    reference_contigs = {}
    metadata_dir = Path(args.preflight_json).parent
    metadata_dir.mkdir(parents=True, exist_ok=True)
    fd, staged_map = tempfile.mkstemp(prefix=".contig-map.", dir=metadata_dir)
    os.close(fd)
    try:
        with open(staged_map, "wb") as raw_map:
            with gzip.GzipFile(fileobj=raw_map, mode="wb", filename="", mtime=0) as compressed:
                with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as text_map:
                    writer = csv.writer(text_map, delimiter="\t")
                    writer.writerow(["seqfile_name", "original_header", "cactus_contig_id", "length_bp"])
                    for name, sample, haplotype, path, source_line in rows:
                        stats, contigs = scan_fasta(path, name, writer)
                        inputs.append({"name": name, "sample": sample, "haplotype": haplotype,
                                       "path": str(path), "source_line": source_line, **stats})
                        if name in {"CHM13", "GRCh38"}:
                            reference_contigs[name] = contigs
        preflight = {
            "preflight_version": 2,
            "included_list": {"path": str(included), "size": list_after.st_size,
                              "mtime_ns": list_after.st_mtime_ns, "sha256": sha256(list_bytes)},
            "inputs": inputs, "seqfile_sha256": sha256(seqfile),
            "haplotype_count": len(selected), "reference_contigs": reference_contigs,
            "total_input_bases": sum(row["total_bases"] for row in inputs),
        }
        preflight["identity_sha256"] = sha256(canonical(preflight))
        previous_path = Path(args.preflight_json)
        locked_path = metadata_dir / "preflight.locked.json"
        state_paths = [metadata_dir / "state.json", metadata_dir / "success.json",
                       *metadata_dir.glob("*.state.json"), *metadata_dir.glob("*.success.json")]
        stage_markers = [*metadata_dir.glob("*.done"),
                         *metadata_dir.glob("*.cactus_finished"),
                         metadata_dir / "settings.locked.json"]
        preflight_hash = sha256(canonical(preflight))
        if locked_path.exists():
            if read_json(locked_path) != preflight:
                raise ValueError("inputs changed for an occupied build ID; choose a new build ID")
        elif any(path.exists() for path in [*state_paths, *stage_markers]):
            raise ValueError("build state exists without its locked preflight; choose a new build ID")
        for state_path in state_paths:
            if state_path.exists() and read_json(state_path).get("preflight_sha256") != preflight_hash:
                raise ValueError("preflight differs from existing build state; choose a new build ID")
        manifest = io.StringIO()
        writer = csv.writer(manifest, delimiter="\t", lineterminator="\n")
        writer.writerow(["seqfile_name", "sample", "haplotype", "source_path",
                         "source_list_line", "size", "mtime_ns", "sha256_uncompressed",
                         "contig_count", "total_bases"])
        for row in inputs:
            writer.writerow([row["name"], row["sample"], row["haplotype"], row["path"],
                             row["source_line"], row["size"], row["mtime_ns"], row["sha256"],
                             row["contig_count"], row["total_bases"]])
        atomic_bytes(args.seqfile, seqfile)
        atomic_bytes(args.input_manifest, manifest.getvalue().encode())
        os.replace(staged_map, args.contig_map)
        atomic_bytes(previous_path, canonical(preflight) + b"\n")
    finally:
        if os.path.exists(staged_map):
            os.unlink(staged_map)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, EOFError, gzip.BadGzipFile) as exc:
        sys.exit(f"Minigraph-Cactus preflight failed: {exc}")
